classdef AudioCD
    % AudioCD class: simulation of the CIRC encoding/decoding and interpolation of the CD digital audio standard
    %
    % Author: Johannes Van Wonterghem, Jan 2017
    %
    % See the static test() method for an example usage of this class
    
    properties
        Fs; % Sample rate of the audio
        
        configuration; % % 0: no CIRC; 1: CIRC as described in standard; 2: Concatenated RS, no interleaving; 3: Single 32,24 RS
        max_interpolation; % The maximum number of interpolated audio samples
        
        enc2; % Outer RS code of the CIRC
        dec2;
        gpoly2;
        
        enc1; % Inner RS code of the CIRC
        dec1;
        gpoly1;
        
        gpoly_8_parity; % Single RS code for configuration 3
        enc_8_parity;
        dec_8_parity;
                
        cd_bits; % Bits written to disk (before EFM)
        
        scaled_quantized_padded_original; % Reference to compare the output of readCD to
    end
    
    methods

        function obj=AudioCD(Fs,configuration,max_interpolation)  
            % Constructor of the AudioCD class
			% INPUT:
            % -Fs: The sample rate of the audio
            % -configuration: 0: no CIRC; 1: CIRC as described in standard; 2: Concatenated RS, no interleaving; 3: Single 32,24 RS
			% -max_interpolation: The maximum number of interpolated audio samples
            % OUTPUT:
            % -obj: the AudioCD object

            obj.Fs = Fs;
            obj.max_interpolation = max_interpolation;
            
            % Initialize the RS encoders and decoders
            primpoly = [1,0,0,0,1,1,1,0,1]; % Primitive polynomial from the standard
            if (configuration == 1) || (configuration == 2)
                obj.gpoly2 = rsgenpoly(255,251,bi2de(fliplr(primpoly)),0);
                obj.enc2 = comm.RSEncoder(255,251,obj.gpoly2,24,'PrimitivePolynomialSource','Property','PrimitivePolynomial',primpoly);
                obj.dec2 = comm.RSDecoder(255,251,obj.gpoly2,24,'PrimitivePolynomialSource','Property','PrimitivePolynomial',primpoly,'ErasuresInputPort',true);

                obj.gpoly1 = rsgenpoly(255,251,bi2de(fliplr(primpoly)),0);
                obj.enc1 = comm.RSEncoder(255,251,obj.gpoly1,28,'PrimitivePolynomialSource','Property','PrimitivePolynomial',primpoly);
                obj.dec1 = comm.RSDecoder(255,251,obj.gpoly1,28,'PrimitivePolynomialSource','Property','PrimitivePolynomial',primpoly,'ErasuresInputPort',true);
            elseif configuration == 3
                obj.gpoly_8_parity = rsgenpoly(255,247,bi2de(fliplr(primpoly)),0);
                obj.enc_8_parity = comm.RSEncoder(255,247,obj.gpoly_8_parity,24,'PrimitivePolynomialSource','Property','PrimitivePolynomial',primpoly);
                obj.dec_8_parity = comm.RSDecoder(255,247,obj.gpoly_8_parity,24,'PrimitivePolynomialSource','Property','PrimitivePolynomial',primpoly);           
            end
            obj.configuration = configuration;
        end
        
        function obj=writeCd(obj,audiofile)
            % Write an audiofile to the CD
            % INPUT:
            % -obj: the AudioCD object
            % -audiofile: [Nsamples x 2] matrix containing the left and right audio track as samples of the double datatype
            % OUTPUT:
            % -obj: the updated AudioCD object            

            assert(size(audiofile,2) == 2);
            
            xscaled = audiofile / max(max(abs(audiofile))); % normalize to -1:1
            x = uencode(xscaled,16); % convert to 16 bit signed values -> THIS IS UNSIGNED!!! Correct way would be x_sign = uencode(xscaled,16,1,'signed');

            xlr16 = reshape(x',[],1); % serialize left and right audio channel -> 2 column L|R to 1 column L,R,L2,R2,R3,L3...
            xlr8 = typecast(xlr16,'uint8'); % split into 8 bit words

            xlr8_padded = [xlr8 ; zeros(24-(rem(numel(xlr8)-1,24)+1),1)]; % pad with zeros to fill an integer number of frames
            n_frames = numel(xlr8_padded)/24; % every frame contains 24 8 bit words
            
            ylr16 = typecast(uint8(xlr8_padded),'uint16');
            y = reshape(ylr16,2,[])';
            obj.scaled_quantized_padded_original = udecode(y,16); % Reference to compare the output of readCD to
            
            switch obj.configuration
                case 0 % no CIRC
                    xlrb = de2bi(xlr8,8);
                case 1 % CIRC as described in standard
                    [delay_interleaved,n_frames] = obj.CIRC_enc_delay_interleave(xlr8_padded,n_frames);
                    [C2_encoded,n_frames] = obj.CIRC_enc_C2(delay_interleaved,n_frames);
                    [delay_unequal,n_frames] = obj.CIRC_enc_delay_unequal(C2_encoded,n_frames);
                    [C1_encoded,n_frames] = obj.CIRC_enc_C1(delay_unequal,n_frames);
                    [delay_inv,n_frames] = obj.CIRC_enc_delay_inv(C1_encoded,n_frames);
                                        
                    xlrb = de2bi(delay_inv,8);
                case 2 % Concatenated RS, no interleaving
                    [C2_encoded,n_frames] = obj.CIRC_enc_C2(xlr8_padded,n_frames);
                    [C1_encoded,n_frames] = obj.CIRC_enc_C1(C2_encoded,n_frames);
                    
                    xlrb = de2bi(C1_encoded,8);
                case 3 % Single 32,24 RS
                    [encoded,n_frames] = obj.C3_enc_8_parity(xlr8_padded,n_frames);
                    
                    xlrb = de2bi(encoded,8);
                otherwise
                    error('Invalid configuration selected');
            end
            
            xlrbserial = reshape(xlrb',[],1);
            
            obj.cd_bits = logical(xlrbserial);            
        end
        
        function obj=bitErrorsCd(obj,p)
            % Add uniform bit errors to cd
            % INPUT:
            % -obj: the AudioCD object
            % -p: the bit error probability, i.e., an obj.cd_bits bit is flipped with probability p
            % OUTPUT:
            % -obj: the updated AudioCD object            

            noise = rand(size(obj.cd_bits))<p;
            
            obj.cd_bits = xor(obj.cd_bits,noise);
        end
        
        function obj=scratchCd(obj,length_scratch,location_scratch)
            % Add a scratch to the cd
            % INPUT:
            % -obj: the AudioCD object
            % -length_scratch: the length of the scratch (in number of bit)
            % -location_scratch: the location of the scratch (in bit offset from start of obj.cd_bits)
            % OUTPUT:
            % -obj: the updated AudioCD object            
            
            obj.cd_bits(location_scratch:min(location_scratch+length_scratch-1,numel(obj.cd_bits))) = 0;
        end
        
        function obj=scratchCd_v2(obj,length_scratch,location_scratch)
            % Add a scratch to the cd
            % INPUT:
            % -obj: the AudioCD object
            % -length_scratch: the length of the scratch (in number of bit)
            % -location_scratch: the location of the scratch (in bit offset from start of obj.cd_bits)
            % OUTPUT:
            % -obj: the updated AudioCD object with the bits inversed at the specified location  
            
            o = obj.cd_bits(location_scratch:min(location_scratch+length_scratch-1,numel(obj.cd_bits)));
            
            inv = ~o;
            
            obj.cd_bits(location_scratch:min(location_scratch+length_scratch-1,numel(obj.cd_bits))) = inv;
        end
        
        function [audio_out,interpolation_flags]=readCd(obj)
            % Read an audiofile from the CD
            % INPUT:
            % -obj: the AudioCD object, containing the data in obj.cd_bits
            % OUTPUT:
            % -audio_out: [Nsamples x 2] matrix containing the left and right audio track as samples of the double datatype  
            % -interpolation_flags: [Nsamples x 2] matrix containing a 0 where no erasure was flagged, a 1 where an erasure was interpolated and a -1
            % where interpolation failed

            
            ylrb = reshape(obj.cd_bits,8,[])';
            ylr8 = bi2de(ylrb);
            
            switch obj.configuration
                case 0 % no CIRC
                    ylr16 = typecast(uint8(ylr8),'uint16');
                    y = reshape(ylr16,2,[])';
                    audio_out = udecode(y,16);
                    interpolation_flags = zeros(size(audio_out));
                case 1
                    n_frames = numel(ylr8)/32;
                    assert(n_frames*32 == numel(ylr8));
                    
                    [delay_inv,n_frames] = obj.CIRC_dec_delay_inv(ylr8,n_frames);
                    [C1_decoded,erasure_flags,n_frames] = obj.CIRC_dec_C1(delay_inv,n_frames);
                    [delay_unequal,erasure_flags,n_frames] = obj.CIRC_dec_delay_unequal(C1_decoded,erasure_flags,n_frames);
                    [C2_decoded,erasure_flags,n_frames] = obj.CIRC_dec_C2(delay_unequal,erasure_flags,n_frames);
                    [deinterleave_delay,erasure_flags,n_frames] = obj.CIRC_dec_deinterleave_delay(C2_decoded,erasure_flags,n_frames);
                                        
                    ylr16 = typecast(uint8(deinterleave_delay),'uint16');
                    y = reshape(ylr16,2,[])';
                    
                    erasure_flags = reshape(erasure_flags,2,[]);
                    erasure_flags = or(erasure_flags(1,:),erasure_flags(2,:))';
                    erasure_flags = reshape(erasure_flags,2,[])';
                    
                    % Linear Interpolation
                    interpolation_failed = zeros(size(erasure_flags),'logical');
                    [y(:,1),interpolation_failed(:,1)] = obj.interpolator(y(:,1),erasure_flags(:,1)); % Left
                    [y(:,2),interpolation_failed(:,2)] = obj.interpolator(y(:,2),erasure_flags(:,2)); % Right
                    
                    audio_out = udecode(y,16);
                    interpolation_flags = zeros(size(audio_out));
                    interpolation_flags(erasure_flags) = 1;
                    interpolation_flags(interpolation_failed) = -1;
                    
                case 2 % Concatenated RS, no interleaving
                    n_frames = numel(ylr8)/32;
                    assert(n_frames*32 == numel(ylr8));
                    
                    [C1_decoded,erasure_flags,n_frames] = obj.CIRC_dec_C1(ylr8,n_frames);
                    erasure_flags_t = erasure_flags;
                    [C2_decoded,erasure_flags,n_frames] = obj.CIRC_dec_C2(C1_decoded,erasure_flags,n_frames);
                    
                    if(numel(erasure_flags)  ~= numel(C2_decoded))
                        disp('Something wrong!');
                    end
                    
                    ylr16 = typecast(uint8(C2_decoded),'uint16');
                    y = reshape(ylr16,2,[])';
                    
                    
                    erasure_flags = reshape(erasure_flags,2,[]);
                    erasure_flags = or(erasure_flags(1,:),erasure_flags(2,:))';
                    erasure_flags = reshape(erasure_flags,2,[])';
                    
                    % Linear Interpolation
                    interpolation_failed = zeros(size(erasure_flags),'logical');
                    [y(:,1),interpolation_failed(:,1)] = obj.interpolator(y(:,1),erasure_flags(:,1)); % Left
                    [y(:,2),interpolation_failed(:,2)] = obj.interpolator(y(:,2),erasure_flags(:,2)); % Right
                    
                    audio_out = udecode(y,16);
                    interpolation_flags = zeros(size(audio_out));
                    interpolation_flags(erasure_flags) = 1;
                    interpolation_flags(interpolation_failed) = -1;
                case 3 % Single 32,24 RS
                    n_frames = numel(ylr8)/32;
                    assert(n_frames*32 == numel(ylr8));
                    
                    [decoded,erasure_flags,n_frames] = obj.C3_dec_8_parity(ylr8,n_frames);
                    
                    ylr16 = typecast(uint8(decoded),'uint16');
                    y = reshape(ylr16,2,[])';
                    
                    erasure_flags = reshape(erasure_flags,2,[]);
                    erasure_flags = or(erasure_flags(1,:),erasure_flags(2,:))';
                    erasure_flags = reshape(erasure_flags,2,[])';
                    
                    % Linear Interpolation
                    interpolation_failed = zeros(size(erasure_flags),'logical');
                    [y(:,1),interpolation_failed(:,1)] = obj.interpolator(y(:,1),erasure_flags(:,1)); % Left
                    [y(:,2),interpolation_failed(:,2)] = obj.interpolator(y(:,2),erasure_flags(:,2)); % Right
                    
                    audio_out = udecode(y,16);
                    interpolation_flags = zeros(size(audio_out));
                    interpolation_flags(erasure_flags) = 1;
                    interpolation_flags(interpolation_failed) = -1;
                otherwise
                    error('Invalid configuration selected');
            end
            
        end
    
        function [output,n_frames] = CIRC_enc_delay_interleave(obj,input,n_frames)
            % CIRC Encoder: Delay of 2 frames + interleaving sequence
            % INPUT:
            % -obj: the AudioCD object
            % -input: the input to this block of the CIRC encoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC encoder  
            % -n_frames: the length of the output expressed in frames (changed from input because of delay!)
            
            %% Create variables to store sequences of the frame symbols
            w12_A = input(1:24:end);%L
            w12_B = input(2:24:end);
            
            w12_1_A = input(3:24:end);%R
            w12_1_B = input(4:24:end);
            
            w12_2_A = input(5:24:end);%L1
            w12_2_B = input(6:24:end);
                     
            w12_3_A = input(7:24:end);%R1
            w12_3_B = input(8:24:end);
                      
            w12_4_A = input(9:24:end);%L2
            w12_4_B = input(10:24:end);
                  
            w12_5_A = input(11:24:end);%R2
            w12_5_B = input(12:24:end);
                    
            w12_6_A = input(13:24:end);%L3
            w12_6_B = input(14:24:end);
                     
            w12_7_A = input(15:24:end);%R3
            w12_7_B = input(16:24:end);
                      
            w12_8_A = input(17:24:end);%L4
            w12_8_B = input(18:24:end);
                      
            w12_9_A = input(19:24:end);%R4
            w12_9_B = input(20:24:end);
            
            w12_10_A = input(21:24:end);%L5
            w12_10_B = input(22:24:end);
            
            w12_11_A = input(23:24:end);%R5
            w12_11_B = input(24:24:end);
            
            %% Delay of 2 frames for even samples
            w12_A = [zeros(2,1); w12_A];%L
            w12_B = [zeros(2,1); w12_B];
            
            w12_1_A = [zeros(2,1); w12_1_A];%R
            w12_1_B = [zeros(2,1); w12_1_B];
            
            w12_4_A = [zeros(2,1); w12_4_A];%L2
            w12_4_B = [zeros(2,1); w12_4_B];
                                             
            w12_5_A = [zeros(2,1); w12_5_A];%R2
            w12_5_B = [zeros(2,1); w12_5_B];
                      
            w12_8_A = [zeros(2,1); w12_8_A];%L4
            w12_8_B = [zeros(2,1); w12_8_B];
                                              
            w12_9_A = [zeros(2,1); w12_9_A];%R4
            w12_9_B = [zeros(2,1); w12_9_B];
            
            %% pad non delayed (odd) symbols with zeros to make them equal length
            % with the delayed symbols
            w12_2_A = [w12_2_A; zeros(2,1)];%L1
            w12_2_B = [w12_2_B; zeros(2,1)];
                     
            w12_3_A = [w12_3_A; zeros(2,1)];%R1
            w12_3_B = [w12_3_B; zeros(2,1)];
                      
            w12_6_A = [w12_6_A; zeros(2,1)];%L3
            w12_6_B = [w12_6_B; zeros(2,1)];
                     
            w12_7_A = [w12_7_A; zeros(2,1)];%R3
            w12_7_B = [w12_7_B; zeros(2,1)];
                      
            w12_10_A = [w12_10_A; zeros(2,1)];%L5
            w12_10_B = [w12_10_B; zeros(2,1)];
            
            w12_11_A = [w12_11_A; zeros(2,1)];%R5
            w12_11_B = [w12_11_B; zeros(2,1)];
            
            %% update amount of frames
            n_frames = n_frames + 2; %2 aditional frames because of delay
            
            %% Combine and interleave symbol sequences
            output = zeros([length(w12_A)*24, 1]);
            
            output(1:24:end) = w12_A;%L
            output(2:24:end) = w12_B;
            
            output(3:24:end) = w12_4_A;%L2
            output(4:24:end) = w12_4_B;
            
            output(5:24:end) = w12_8_A;%L4
            output(6:24:end) = w12_8_B;
            
            output(7:24:end) = w12_1_A;%R
            output(8:24:end) = w12_1_B;
			
			output(9:24:end) = w12_5_A;%R2
            output(10:24:end) = w12_5_B;
			
			output(11:24:end) = w12_9_A;%R4
            output(12:24:end) = w12_9_B;
			
			output(13:24:end) = w12_2_A;%L1
            output(14:24:end) = w12_2_B;
			
			output(15:24:end) = w12_6_A;%L3
            output(16:24:end) = w12_6_B;
			
			output(17:24:end) = w12_10_A;%L5
            output(18:24:end) = w12_10_B;
			
			output(19:24:end) = w12_3_A;%R1
            output(20:24:end) = w12_3_B;
			
			output(21:24:end) = w12_7_A;%R3
            output(22:24:end) = w12_7_B;
			
			output(23:24:end) = w12_11_A;%R5
            output(24:24:end) = w12_11_B;
            
            
          
        end
        
        function [output,n_frames] = CIRC_enc_C2(obj,input,n_frames)
            % CIRC Encoder: Generation of 4 parity symbols (C2)
            % INPUT:
            % -obj: the AudioCD object
            % -input: the input to this block of the CIRC encoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC encoder  
            % -n_frames: the length of the output expressed in frames
            
            output = zeros(n_frames*28,1,'uint8');
            for i = 1:n_frames
                temp = step(obj.enc2,input((i-1)*24+1:i*24));
                
                % place Q parity symbols of C2 in the midle
                outputframe = zeros(1, 28)';
                outputframe(1:12) = temp(1:12);
                outputframe(13:16) = temp(25:28);
                outputframe(17:28) = temp(13:24);
                output((i-1)*28+1:i*28) = outputframe;
            end
        end
        
        function [output,n_frames] = CIRC_enc_delay_unequal(obj,input,n_frames)
            % CIRC Encoder: Delay lines of unequal length
            % INPUT:
            % -obj: the AudioCD object
            % -input: the input to this block of the CIRC encoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC encoder  
            % -n_frames: the length of the output expressed in frames (changed from input because of delay!)
            
            %% get sequence per symbol            
            symbols = zeros(n_frames,28);
            for i = 1:28
                symbols(:,i) = input(i:28:end);
            end
            
            %% delay frames with unequal length with minimal of 4 frames except for the first symbol.
            n_frames = n_frames+4*27;
            symbols_delayed = zeros(n_frames,28);
            for i = 1:28
                symbols_delayed(:,i) = [zeros((i-1)*4,1); symbols(:,i); zeros((28-i)*4,1)];
            end
            
            %% Combine sequences
            output = zeros(n_frames*28,1);
            for i = 1:28
                output(i:28:end) = symbols_delayed(:,i);
            end

        end
        
        function [output,n_frames] = CIRC_enc_C1(obj,input,n_frames)
            % CIRC Encoder: Generation of 4 parity symbols (C1)
            % INPUT:
            % -obj: the AudioCD object
            % -input: the input to this block of the CIRC encoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC encoder  
            % -n_frames: the length of the output expressed in frames
            
            output = zeros(n_frames*32,1,'uint8');
            for i = 1:n_frames
                output((i-1)*32+1:i*32) = step(obj.enc1,input((i-1)*28+1:i*28));
            end
            

        end
        
        function [output,n_frames] = CIRC_enc_delay_inv(obj,input,n_frames)
            % CIRC Encoder: Delay of 1 frame + inversions
            % INPUT:
            % -obj: the AudioCD object
            % -input: the input to this block of the CIRC encoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC encoder  
            % -n_frames: the length of the output expressed in frames (changed from input because of delay!)
            
            %% get sequence per symbol            
            symbols = zeros(n_frames,32);
            for i = 1:32
                symbols(:,i) = input(i:32:end);
            end
            
            %% delay even symbols with 1 frame 
            n_frames = n_frames+1;
            symbols_delayed = zeros(n_frames,32);
            for i = 1:2:32
                symbols_delayed(:,i) = [0; symbols(:,i)];
            end
            for i = 2:2:32
                symbols_delayed(:,i) = [symbols(:,i); 0];
            end
            
            %% Do inversions on parity symbols and combine sequences
            output = zeros(n_frames*32,1,'uint8');
            for i = 1:32
                if (((i > 12) && (i < 17)) || (i > 28))
                    output(i:32:end) = bitcmp(uint8(symbols_delayed(:,i)));
                else
                    output(i:32:end) = symbols_delayed(:,i);
                end
            end
        end
        
        function [output,n_frames] = CIRC_dec_delay_inv(obj,input,n_frames)
            % CIRC Decoder: Delay of 1 frame + inversions
            % INPUT:
            % -obj: the AudioCD object
            % -input: the input to this block of the CIRC decoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC decoder  
            % -n_frames: the length of the output expressed in frames (changed from input because of delay!)
            
            % Note: remove empty frames such that obj.CIRC_dec_delay_inv(obj.CIRC_enc_delay_inv(x)) == x!
            
            %% get sequence per symbol            
            symbols = zeros(n_frames,32);
            for i = 1:32
                symbols(:,i) = input(i:32:end);
            end
            
            %% delay odd symbols with 1 frame 
            n_frames = n_frames+1;
            symbols_delayed = zeros(n_frames,32);
            for i = 2:2:32
                delay_symbol = 0;
                if (((i > 12) && (i < 17)) || (i > 28))% because of inverted parity need to add 1's instead of zeros for delay
                    delay_symbol = 255;
                end
                symbols_delayed(:,i) = [delay_symbol; symbols(:,i)];
            end
            for i = 1:2:32
                delay_symbol = 0;
                if (((i > 12) && (i < 17)) || (i > 28))% because of inverted parity need to add 1's instead of zeros for delay
                    delay_symbol = 255;
                end
                symbols_delayed(:,i) = [symbols(:,i); delay_symbol];
            end
            
            %% Do inversions on parity symbols
            for i = 1:32
                if (((i > 12) && (i < 17)) || (i > 28))
                    symbols_delayed(:,i) = bitcmp(uint8(symbols_delayed(:,i)));
                end
            end
            
            %% check if first and last frame is empty (it should be because of delay at end of encoder part and this delay)
            force_removal = 1;
            row_start = 0;
            row_end = 0;
            if (force_removal == 1)
                row_start = 1;
                row_end = 1;
                n_frames = n_frames - 2;
            else
                if (symbols_delayed(1,:) == 0)
                    row_start = row_start + 1;
                    n_frames = n_frames -1;
                end
                if (symbols_delayed(end,:) == 0)
                    row_end = row_end + 1;
                    n_frames = n_frames -1;
                end
            end
            
            %% Combine sequences with respect to empty first and last frame
            output = zeros(n_frames*32,1,'uint8');
            for i = 1:32
                output(i:32:end) = symbols_delayed((1+row_start):(end-row_end),i);
            end

        end
        
        
        function [output,erasure_flags_out,n_frames] = CIRC_dec_C1(obj,input,n_frames)
            % CIRC Decoder: C1 decoder
            % INPUT:
            % -obj: the AudioCD object
            % -input: the input to this block of the CIRC decoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC decoder  
            % -erasure_flags_out: the erasure flags at the output of this block, follow the decoding algorithm from the assignment
            % -n_frames: the length of the output expressed in frames
            
            % Decode according to decoding strategy algorithm 1
            output = zeros(n_frames*28,1,'uint8');
            erasure_flags_out = zeros(n_frames*28,1,'logical');
            for i = 1:n_frames
                output_dec = 0;
                ERR = 0;
                try
                [output_dec,ERR] = step(obj.dec1,input((i-1)*32+1:i*32),zeros(32,1));
                catch ME
                    if isequal(ME.message,'The sum of the number of erasures and punctures cannot exceed twice the error correcting capability of the code.')
                        ERR = -1;
                    else
                        rethrow(ME);
                    end
                end
                
                if ((ERR == 0) || (ERR == 1)) % modify symbol accordingly
                    output((i-1)*28+1:i*28) = output_dec;
%                 elseif ERR > 1 % assign erasure flags to all symbols and do not correct
%                     for e = 1:28
%                         erasure_flags_out((i-1)*28+1+(e-1)) = logical(sum(bitxor(output_dec(e),input((i-1)*32+1+(e-1)))));
%                     end
%                     output((i-1)*28+1:i*28) = input((i-1)*32+1:(i-1)*32+28);
                else % uncorectable error so no correction
                    output((i-1)*28+1:i*28) = input((i-1)*32+1:(i-1)*32+28);
                    erasure_flags_out((i-1)*28+1:i*28) = 1;
                end                
            end
        end

        
        function [output,erasure_flags_out,n_frames] = CIRC_dec_delay_unequal(obj,input,erasure_flags_in,n_frames)
            % CIRC Decoder: Delay lines of unequal length
            % INPUT:
            % -obj: the AudioCD object
            % -input: the input to this block of the CIRC decoder
            % -erasure_flags_in: the erasure flags at the input of this block of the CIRC decoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC decoder  
            % -erasure_flags_out: the erasure flags at the output of this block (erasure flags follow same interleaving as the data)
            % -n_frames: the length of the output expressed in frames
            
            % Note: remove empty frames such that obj.CIRC_dec_delay_unequal(obj.CIRC_enc_delay_unequal(x)) == x!
            
            %% get sequence per symbol
            % do same for erasures
            symbols = zeros(n_frames,28);
            erasures = zeros(n_frames,28,'logical');
            for i = 1:28
                symbols(:,i) = input(i:28:end);
                erasures(:,i) = erasure_flags_in(i:28:end);
            end
            
            %% delay frames with unequal length with minimal of 4 frames except for the last symbol.
            %% do same for erasure_flags
            n_frames = n_frames+4*27;
            symbols_delayed = zeros(n_frames,28);
            erasures_delayed = zeros(n_frames,28,'logical');
            for i = 1:28
                symbols_delayed(:,29-i) = [zeros((i-1)*4,1); symbols(:,29-i); zeros((28-i)*4,1)];
                erasures_delayed(:,29-i) = [zeros((i-1)*4,1); erasures(:,29-i); zeros((28-i)*4,1)];
            end
            
            %% check for empty frames (it should be because of delays at end of encoder part and this delays->27*4 frames front and back)
            force_removal = 1;
            row_start = 0;
            row_end = 0;
            if (force_removal == 1)
                row_start = 27*4;
                row_end = 27*4;
                n_frames = n_frames - 2*27*4;
            else
                for i = 1:(27*4)
                    if (symbols_delayed(i,:) == 0)
                        row_start = row_start + 1;
                        n_frames = n_frames -1;
                    end
                    if (symbols_delayed(end-(i-1),:) == 0)
                        row_end = row_end + 1;
                        n_frames = n_frames -1;
                    end
                end
            end            
            
            %% Combine sequences with respect to empty frames
            % do same for erasures
            output = zeros(n_frames*28,1,'uint8');
            erasure_flags_out = zeros(n_frames*28,1,'logical');
            for i = 1:28
                output(i:28:end) = symbols_delayed((1+row_start):(end-row_end),i);
                erasure_flags_out(i:28:end) = erasures_delayed((1+row_start):(end-row_end),i);
            end
            

        end

        
        function [output,erasure_flags_out,n_frames] = CIRC_dec_C2(obj,input,erasure_flags_in,n_frames)
            % CIRC Decoder: C2 decoder
            % INPUT:
            % -obj: the AudioCD object
            % -input: the input to this block of the CIRC decoder
            % -erasure_flags_in: the erasure flags at the input of this block of the CIRC decoder            
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC decoder  
            % -erasure_flags_out: the erasure flags at the output of this block, follow the decoding algorithm from the assignment
            % -n_frames: the length of the output expressed in frames
            
            output = zeros(n_frames*24,1,'uint8');
            erasure_flags_out = zeros(n_frames*24,1,'logical');
            
            for i = 1:n_frames
                % place Q parity symbols of C2 from the midle to the back
                % again, same for erasure flags
                temp = input((i-1)*28+1:i*28);
                input_frame = zeros(28,1);
                input_frame(1:12) = temp(1:12);
                input_frame(13:24) = temp(17:28);
                input_frame(25:28) = temp(13:16);
                temp_erasure = erasure_flags_in((i-1)*28+1:i*28);
                input_flags = zeros(28,1,'logical');
                input_flags(1:12) = temp_erasure(1:12);
                input_flags(13:24) = temp_erasure(17:28);
                input_flags(25:28) = temp_erasure(13:16);
                
                % Decode acording to decoding strategy of algorithm 1
                output_dec = 0;
                ERR = 0;
                try
                [output_dec,ERR] = step(obj.dec2,input_frame,input_flags);
                catch ME
                    if isequal(ME.message,'The sum of the number of erasures and punctures cannot exceed twice the error correcting capability of the code.')
                        ERR = -1;
                    else
                        rethrow(ME);
                    end
                end
                
                if ((ERR == 0) || (ERR == 1)) % modify symbol accordingly
                    output((i-1)*24+1:i*24) = output_dec;
                else
                    if (sum(input_flags)) > 2 % copy erasure flags from C1 and do not correct
						erasure_flags_out((i-1)*24+1:i*24) = input_flags(1:24);
						output((i-1)*24+1:i*24) = input_frame(1:24);
                    elseif ((sum(input_flags) == 2) && (ERR > -1)) % modify symbol accordingly when C1 has only 2 flags and error correcting successful
                        output((i-1)*24+1:i*24) = output_dec;
                    else% uncorectable error so no correction
                        erasure_flags_out((i-1)*24+1:i*24) = 1;
                        output((i-1)*24+1:i*24) = input_frame(1:24);
%                     else % assign erasure flags to all symbols and do not correct
%                         output((i-1)*24+1:i*24) = input_frame(1:24);                    
%                         for e = 1:24
%                             erasure_flags_out((i-1)*24+1+(e-1)) = logical(sum(bitxor(output_dec(e),input_frame(e))));
%                         end
                    end
                end
            end        
        end
        
        
        function [output,erasure_flags_out,n_frames] = CIRC_dec_deinterleave_delay(obj,input,erasure_flags_in,n_frames)
            % CIRC Decoder: De-interleaving sequence + delay of 2 frames
            % INPUT:
            % -obj: the AudioCD object
            % -input: the input to this block of the CIRC decoder
            % -erasure_flags_in: the erasure flags at the input of this block of the CIRC decoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC decoder  
            % -erasure_flags_out: the erasure flags at the output of this block (eraure flags follow same interleaving as the data)
            % -n_frames: the length of the output expressed in frames
            
            % Note: remove empty frames such that obj.CIRC_dec_deinterleave_delay(obj.CIRC_enc_delay_interleave(x)) == x!
            
            %% deinterleave
            % do same for erasures
            w12_A = input(1:24:end);%L
            w12_B = input(2:24:end);
			w12_A_e = erasure_flags_in(1:24:end);
			w12_B_e = erasure_flags_in(2:24:end);
            
            w12_4_A = input(3:24:end);%L2
            w12_4_B = input(4:24:end);
			w12_4_A_e = erasure_flags_in(3:24:end);
            w12_4_B_e = erasure_flags_in(4:24:end);  
			
            w12_8_A = input(5:24:end);%L4
            w12_8_B = input(6:24:end);
			w12_8_A_e = erasure_flags_in(5:24:end);
			w12_8_B_e = erasure_flags_in(6:24:end);
                      
            w12_1_A = input(7:24:end);%R
            w12_1_B = input(8:24:end);
			w12_1_A_e = erasure_flags_in(7:24:end);
			w12_1_B_e = erasure_flags_in(8:24:end);
			
			w12_5_A = input(9:24:end);%R2
            w12_5_B = input(10:24:end);
			w12_5_A_e = erasure_flags_in(9:24:end);
            w12_5_B_e = erasure_flags_in(10:24:end);
			
			
			w12_9_A = input(11:24:end);%R4
            w12_9_B = input(12:24:end);
			w12_9_A_e = erasure_flags_in(11:24:end);
            w12_9_B_e = erasure_flags_in(12:24:end);
			
			          
			w12_2_A = input(13:24:end);%L1
            w12_2_B = input(14:24:end);
			w12_2_A_e = erasure_flags_in(13:24:end);
			w12_2_B_e = erasure_flags_in(14:24:end);
			          
			w12_6_A = input(15:24:end);%L3
            w12_6_B = input(16:24:end);
			w12_6_A_e = erasure_flags_in(15:24:end);
			w12_6_B_e = erasure_flags_in(16:24:end);
			
			w12_10_A = input(17:24:end);%L5
            w12_10_B = input(18:24:end);
			w12_10_A_e = erasure_flags_in(17:24:end);
			w12_10_B_e = erasure_flags_in(18:24:end);
			
			w12_3_A = input(19:24:end);%R1
            w12_3_B = input(20:24:end);
			w12_3_A_e = erasure_flags_in(19:24:end);
			w12_3_B_e = erasure_flags_in(20:24:end);          
			
			w12_7_A = input(21:24:end);%R3
            w12_7_B = input(22:24:end);
			w12_7_A_e = erasure_flags_in(21:24:end);
			w12_7_B_e = erasure_flags_in(22:24:end);
			
			w12_11_A = input(23:24:end);%R5
            w12_11_B = input(24:24:end);
			w12_11_A_e = erasure_flags_in(23:24:end);
			w12_11_B_e = erasure_flags_in(24:24:end);
            
            deinterleaved = zeros(n_frames*24,1,'uint8');
			e_deinterleaved = zeros(n_frames*24,1,'logical');
            
            deinterleaved(1:24:end) = w12_A;%L
            deinterleaved(2:24:end) = w12_B;
			e_deinterleaved(1:24:end) = w12_A_e;
			e_deinterleaved(2:24:end) = w12_B_e;
            
            deinterleaved(3:24:end) = w12_1_A;%R
            deinterleaved(4:24:end) = w12_1_B;
			e_deinterleaved(3:24:end) = w12_1_A_e;
			e_deinterleaved(4:24:end) = w12_1_B_e;
            
            deinterleaved(5:24:end) = w12_2_A;%L1
            deinterleaved(6:24:end) = w12_2_B;
			e_deinterleaved(5:24:end) = w12_2_A_e;
            e_deinterleaved(6:24:end) = w12_2_B_e;
            
			deinterleaved(7:24:end) = w12_3_A;%R1
            deinterleaved(8:24:end) = w12_3_B;
			e_deinterleaved(7:24:end) = w12_3_A_e;
			e_deinterleaved(8:24:end) = w12_3_B_e;
			
			deinterleaved(9:24:end) = w12_4_A;%L2
            deinterleaved(10:24:end) = w12_4_B;
			e_deinterleaved(9:24:end) = w12_4_A_e;
			e_deinterleaved(10:24:end) = w12_4_B_e;
			
			deinterleaved(11:24:end) = w12_5_A;%R2
            deinterleaved(12:24:end) = w12_5_B;
			e_deinterleaved(11:24:end) = w12_5_A_e;
			e_deinterleaved(12:24:end) = w12_5_B_e;
			
			deinterleaved(13:24:end) = w12_6_A;%L3
            deinterleaved(14:24:end) = w12_6_B;
			e_deinterleaved(13:24:end) = w12_6_A_e;
			e_deinterleaved(14:24:end) = w12_6_B_e;
			
			deinterleaved(15:24:end) = w12_7_A;%R3
            deinterleaved(16:24:end) = w12_7_B;
			e_deinterleaved(15:24:end) = w12_7_A_e;
			e_deinterleaved(16:24:end) = w12_7_B_e;
			
			deinterleaved(17:24:end) = w12_8_A;%L4
            deinterleaved(18:24:end) = w12_8_B;
			e_deinterleaved(17:24:end) = w12_8_A_e;
			e_deinterleaved(18:24:end) = w12_8_B_e;
			
			deinterleaved(19:24:end) = w12_9_A;%R4
            deinterleaved(20:24:end) = w12_9_B;
			e_deinterleaved(19:24:end) = w12_9_A_e;
			e_deinterleaved(20:24:end) = w12_9_B_e;
			
			deinterleaved(21:24:end) = w12_10_A;%L5
            deinterleaved(22:24:end) = w12_10_B;
			e_deinterleaved(21:24:end) = w12_10_A_e;
			e_deinterleaved(22:24:end) = w12_10_B_e;
			
			deinterleaved(23:24:end) = w12_11_A;%R5
            deinterleaved(24:24:end) = w12_11_B;
			e_deinterleaved(23:24:end) = w12_11_A_e;
            e_deinterleaved(24:24:end) = w12_11_B_e;
            
            %% get sequence per symbol
            % do same for erasures
            symbols = zeros(n_frames,24);
            erasures = zeros(n_frames,24,'logical');
            for i = 1:24
                symbols(:,i) = deinterleaved(i:24:end);
                erasures(:,i) = e_deinterleaved(i:24:end);
            end
            
            %% delay odd samples/words with 2 frames (L1;R1...) 
            n_frames = n_frames+2;
            deinterleaved_delayed = zeros(n_frames,24);
            erasures_deinter_delay = zeros(n_frames,24,'logical');
            for i = 1:24
                if (((i > 4) && (i < 9)) || ((i > 12) && (i < 17)) || (i > 20))
                    deinterleaved_delayed(:,i) = [0; 0; symbols(:,i)];
                    erasures_deinter_delay(:,i) = [0; 0; erasures(:,i)];
                else
                    deinterleaved_delayed(:,i) = [symbols(:,i); 0; 0];
                    erasures_deinter_delay(:,i) = [erasures(:,i); 0; 0];
                end
            end
            
            %% check for empty frames (it should be because of delays at end of encoder part and this delays->27*4 frames front and back)
            force_removal = 1;
            row_start = 0;
            row_end = 0;
            if (force_removal == 1)
                row_start = 2;
                row_end = 2;
                n_frames = n_frames - 2*2;
            else
                for i = 1:2
                    if (deinterleaved_delayed(i,:) == 0)
                        row_start = row_start + 1;
                        n_frames = n_frames -1;
                    end
                    if (deinterleaved_delayed(end-(i-1),:) == 0)
                        row_end = row_end + 1;
                        n_frames = n_frames -1;
                    end
                end
            end            
            
            %% Combine sequences with respect to empty frames
            % do same for erasures
            output = zeros(n_frames*24,1,'uint8');
            erasure_flags_out = zeros(n_frames*24,1,'logical');
            for i = 1:24
                output(i:24:end) = deinterleaved_delayed((1+row_start):(end-row_end),i);
                erasure_flags_out(i:24:end) = erasures_deinter_delay((1+row_start):(end-row_end),i);
            end
            
        end
        
        function [output,n_frames] = C3_enc_8_parity(obj,input,n_frames)
            % Configuration 3: Generation of 8 parity symbols
            % INPUT:
            % -obj: the AudioCD object
            % -input: the input to this block
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block
            % -n_frames: the length of the output expressed in frames
            
            output = zeros(n_frames*32,1,'uint8');
            for i = 1:n_frames
                output((i-1)*32+1:i*32) = step(obj.enc_8_parity,input((i-1)*24+1:i*24));
            end
        end
        
        function [output,erasure_flags_out,n_frames] = C3_dec_8_parity(obj,input,n_frames)
            % Configuration 3: Decoder
            % INPUT:
            % -obj: the AudioCD object
            % -input: the input to this block
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block
            % -erasure_flags_out: the erasure flags at the output of this block
            % -n_frames: the length of the output expressed in frames

            output = zeros(n_frames*24,1,'uint8');
            erasure_flags_out = zeros(n_frames*24,1,'logical');
            for i = 1:n_frames
                [output_dec,ERR] = step(obj.dec_8_parity,input((i-1)*32+1:i*32));
                if ERR == -1
                    output((i-1)*24+1:i*24) = output_dec;
                    erasure_flags_out((i-1)*24+1:i*24) = 1;
                else
                    output((i-1)*24+1:i*24) = output_dec;
                end
                
            end
        end
        
        function [output,interpolation_failed] = interpolator(obj,input,erasure_flags_in)
            % Interpolation: Linear interpolation
            % INPUT:
            % -obj: the AudioCD object
            % -input: the input to this block
            % -erasure_flags_in: the erasure flags at the input of this block
            % OUTPUT:
            % -output: linear interpolation of the input where there are no more than obj.max_interpolation consecutive erasures
            % -interpolation_failed: equal to one at the samples where interpolation failed

			% Make sure first and last sample are defined. If erasure: set to 0 (corresponds to the unsigned integer 2^15)
            if erasure_flags_in(1) ~= 0
                erasure_flags_in(1) = 0;
                input(1) = 2^15;
            end
            if erasure_flags_in(end) ~= 0
                erasure_flags_in(end) = 0;
                input(end) = 2^15;
            end
            
            erasure_burst = zeros(size(erasure_flags_in(:))); % Number of consecutive erasures
            ii = strfind([0,erasure_flags_in(:)'],[0 1]);
            erasure_burst(ii) = strfind([erasure_flags_in(:)',0],[1 0]) - ii + 1;
            
            output = input;
            interpolation_failed = erasure_flags_in;
            for i = find( (erasure_burst>0) & (erasure_burst<= obj.max_interpolation))'
                output(i:i+erasure_burst(i)-1) = uint16(round(double(output(i-1))+(1:erasure_burst(i))*(double(output(i+erasure_burst(i)))-double(output(i-1)))/(erasure_burst(i)+1)));
                interpolation_failed(i:i+erasure_burst(i)-1) = 0;
            end


        end
        
    end
    
    methods(Static)
        function test()
            % Test the Matlab code of this class
            
%             audio_file = audioread('Hallelujah.wav');
%             Fs = 44.1e3;
            
            audio_file = audioread('Hallelujah_22050.wav');
            Fs = 22.05e3;


            cd = AudioCD(Fs,1,8);
            tic
            cd = cd.writeCd(audio_file);
            toc
            
%             cd = cd.scratchCd_v2(3640,36065); %29th symbol of a frame
            
%             T_scratch = 600000; % Scratch at a diameter of approx. 66 mm
%             for i = 1:floor(numel(cd.cd_bits)/T_scratch)
%                 cd = cd.scratchCd(3000,30000+(i-1)*T_scratch);
%             end
            cd = cd.bitErrorsCd(0.005);
            tic
            [out,interpolation_flags] = cd.readCd();
            toc

            fprintf('Number samples with erasure flags: %d\n',sum(sum(interpolation_flags~=0)));
            fprintf('Number samples with failed interpolations: %d\n',sum(sum(interpolation_flags==-1)));
            fprintf('Number undetected errors: %d\n',sum(sum(out(interpolation_flags==0) ~= cd.scaled_quantized_padded_original(interpolation_flags==0))));

            sound(out,Fs);     
            audiowrite('Hallelujah_22050_output.wav',out,Fs);
        end
        
        % Method for question 5.3.5, allows to choose the length of the
        % scratch, the error probability and the configuration
        function [num_samples, erasure_flags_count, failed_interpolations_count, undetected_error_count] = test_conf(configuration, scratch_length, error_probability, play_sound, write_file, show_log)
            % Test the Matlab code of this class
            
%             audio_file = audioread('Hallelujah.wav');
%             Fs = 44.1e3;

            if (show_log)
                fprintf("=====\nConfiguration: %d\nScratch length: %d\nError probability: %.3f\n", configuration, scratch_length, error_probability);
            end
            
            audio_file = audioread('Hallelujah_22050.wav');
            Fs = 22.05e3;

            cd = AudioCD(Fs,configuration,8);
            tic
            cd = cd.writeCd(audio_file);
            toc
            
            T_scratch = 600000; % Scratch at a diameter of approx. 66 mm
            for i = 1:floor(numel(cd.cd_bits)/T_scratch)
                cd = cd.scratchCd(scratch_length,30000+(i-1)*T_scratch);
            end
            cd = cd.bitErrorsCd(error_probability);
            tic
            [out,interpolation_flags] = cd.readCd();
            toc
            
            num_samples = 2 * size(interpolation_flags, 1);
            erasure_flags_count = sum(sum(interpolation_flags ~= 0));
            failed_interpolations_count = sum(sum(interpolation_flags == -1));
            undetected_error_count = sum(sum(out(interpolation_flags == 0) ~= cd.scaled_quantized_padded_original(interpolation_flags == 0)));

            if (show_log)
                fprintf('Number samples with erasure flags: %d\n',erasure_flags_count);
                fprintf('Number samples with failed interpolations: %d\n',failed_interpolations_count);
                fprintf('Number undetected errors: %d\n',undetected_error_count);
            end

            filename = sprintf("conf_%d-scratchlength_%d-error_probability_%.3f.wav", configuration, scratch_length, error_probability);
            
            if (play_sound)
                sound(out,Fs);
            end
            
            if (write_file)
                audiowrite(filename,out,Fs);
            end
        end
        
        % generates all audiofiles for question 5.3.5, part 1
        function testSuite1()
            diary("command_output_suite_1.log")
            diary on
            erasure_probs = zeros(3, 4);
            undetected_errors = zeros(3, 4);
            for conf = 0:3
                index = 1;
                for scratch_length = [100, 3000, 10000]
                   [ns, ef, fi, ue] = AudioCD.test_conf(conf, scratch_length, 0, false, true, true);
                   erasure_probs(index, conf+1) = ef;
                   undetected_errors(index, conf+1) = ue;
                   index = index + 1;
                end
            end
            
            figure
            bar(erasure_probs);
            title("Number of samples marked as erasure");
            xlabel("Scrath length (bit)");
            ylabel("Number of samples marked as erasure");
            ax = gca;
            set(gca, 'XTickLabel', {'100', '3000', '10000'});
            set(gcf, 'position', [10, 10, 800, 400]);
            ax.YLim = [0, 15000];
            legend(["Configuration 0", "Configuration 1", "Configuration 2", "Configuration 3"], "location", "northwest");
            saveas(gcf, "erasure_flags_bar.png");
            savefig(gcf, "erasure_flags_bar.fig");
            
            figure
            bar(undetected_errors);
            title("Number of undetected errors");
            xlabel("Scrath length (bit)");
            ylabel("Number of undetected errors");
            ax = gca;
            set(gca, 'XTickLabel', {'100', '3000', '10000'});
            set(gcf, 'position', [10, 10, 800, 400]);
            ax.YLim = [0, 7500];
            legend(["Configuration 0", "Configuration 1", "Configuration 2", "Configuration 3"], "location", "northwest");
            saveas(gcf, "failed_interpolation_bar.png");
            savefig(gcf, "failed_interpolation_bar.fig");
            
            diary off
        end
        
        % plots data for question 5.3.5 part 2
        function testSuite2()
            diary("command_output_suite_2.log")
            diary on
            colors = ['r', 'g', 'b'];
                
            x = (logspace(-1-log10(2), -3, 10));
            y1 = cell(1, length(colors));
            y2 = cell(1, length(colors));
            
            for conf = 1:3
                i = 0;
                erasure_probs = zeros(size(logspace(-1-log10(2), -3, 10), 2), 1);
                failed_interpolation_probs = zeros(size(logspace(-1-log10(2), -3, 10), 2), 1);
                
                for error_prob = logspace(-1-log10(2), -3, 10)
                    i = i + 1;
                    [ns, ef, fi, ue] = AudioCD.test_conf(conf, 0, error_prob, false, false, true);
                    erasure_probs(i) = ef / ns;
                    failed_interpolation_probs(i) = fi / ns;
                end
                
                y1{conf} = erasure_probs;
                
                figure
                semilogx(x, erasure_probs);
                title(sprintf("Percentage of samples marked as erasure (configuration %d)", conf));
                xlabel("Error probability in source");
                ylabel("Probability that a sample is marked as erasure");
                set(gcf, 'position', [10, 10, 600, 400]);
                saveas(gcf, sprintf("erasure_flags_probs_conf_%d.png", conf));
                savefig(gcf, sprintf("erasure_flags_probs_conf_%d.fig", conf));
                
                y2{conf} = failed_interpolation_probs;
                
                figure
                semilogx(x, failed_interpolation_probs);
                title(sprintf("Percentage of samples that failed interpolation (configuration %d)", conf));
                xlabel("Error probability in source");
                ylabel("Probability of a failed interpolation");
                set(gcf, 'position', [10, 10, 600, 400]);
                saveas(gcf, sprintf("failed_interpolations_conf_%d.png", conf));
                savefig(gcf, sprintf("failed_interpolations_conf_%d.fig", conf));
            end
            
            figure
            for conf = 1:3
                hold on
                semilogx(x, y1{conf}, 'color', colors(conf), 'LineWidth', 2);
                grid on
            end
            
            title("Percentage of samples marked as erasure");
            xlabel("Error probability in source");
            ylabel("Probability that a sample is marked as erasure");
            set(gcf, 'position', [10, 10, 600, 400]);
            ax = gca;
            ax.YLim = [0.0, 1.0];
            set(gca, 'XScale', 'log');
            legend(["Configuration 1", "Configuration 2", "Configuration 3"], "location", "southeast");
            saveas(gcf, "erasure_flags_probs.png");
            savefig(gcf, "erasure_flags_probs.fig");
            
            figure
            for conf = 1:3
                hold on
                semilogx(x, y2{conf}, 'color', colors(conf), 'LineWidth', 2);
                grid on
            end
            
            title("Percentage of samples that failed interpolation");
            xlabel("Error probability in source");
            ylabel("Probability of a failed interpolation");
            set(gcf, 'position', [10, 10, 600, 400]);
            ax = gca;
            ax.YLim = [0.0, 1.0];
            set(gca, 'XScale', 'log');
            legend(["Configuration 1", "Configuration 2", "Configuration 3"], "location", "southeast");
            saveas(gcf, "failed_interpolations.png");
            savefig(gcf, "failed_interpolations.fig");
            
            diary off
        end
    end
    
end
