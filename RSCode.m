classdef RSCode
	% RSCode class: allows encoding and decoding using a Reed-Solomon code
	%
	% Author: Johannes Van Wonterghem, Jan 2017
	%
	% See the static test() method for an example usage of this class
	
	properties
		m; % GF(2^m) field
		
		n; % Code length
		k; % Information length
		t; % Error correction capability
		l; % Shortened information length (-> shortened code length = l+n-k)
		
		m0; % m0 of the Reed-Solomon code, determines first root of generator
		
		g; % rowvector containing the GF(2^m) generator polynomial coefficients in order of descending powers
	end
	
	methods
	
		function obj=RSCode(m,t,l,m0)
			% Constructor of the RSCode class
			% INPUT:
			% -m: Defines the GF(2^m) field
			% -t: Error correction capability
			% -l: Shortened information length (-> shortened code length = l+n-k)
			% -m0: m0 of the Reed-Solomon code, determines first root of generator
			% OUTPUT:
			% -obj: the RSCode object
			obj.m = m;
		  
			obj.t = t;
			obj.l = l;
			
			obj.n = 2^m-1;
			obj.k = obj.n-2*t;
		   
			obj.m0 = m0;
			
			obj.g = RSCode.makeGenerator(m,t,m0);

		end
		
		function code = encode(obj,msg)
			% Systematically encode information words using the Reed-Solomon code
			% INPUT:
			% -obj: the RSCode object, defines the properties of the code
			% -msg: every row contains a length obj.l information word consisting of GF(2^m) elements
			% OUTPUT:
			% -code: every row contains a GF(2^m) codeword corresponding to systematic Reed-Solomon coding of the corresponding information word			
			
            % we assert we have the right size of message comming in
			assert(size(msg,2) == obj.l);
			
            % this method of encoding is based on the syllabus page 145
            % we could just multiply the codeword by the generator
            % polynomial, we have to divide anyway to be sure its a
            % codeword, and this way we have a simpler convolution on paper
            
            %page 145 in the syllabus
            %a = msg -> vermenigvuldigen met x^(n-k)
            %b = rest van msg delen door g(x)
            %encoded = a + b
            
            % we define a polynomial with one term: x^(n-k)
            extender = [gf(1, obj.m) gf(zeros(1, obj.n-obj.k), obj.m)];
          
            %define code as the output matrix where every row is a codeword
            code = gf(zeros(size(msg,1), obj.l+2*obj.t),8);

            for i=1:size(msg,1)
                % we multiply the extender with our message, we can do this
                % with a convolution. This just raises the degree of our
                % message, which is inefficient to do this with a matlab
                % convolution but its what you would do on paper, so this reads
                % more clearly
                extended_message = conv(msg(i,:) , extender);
                
                % To ensure we result in a codeword, we divide by the
                % generator polynomial
                [~, remainder] = deconv(extended_message, obj.g);
            
                % We add the extended message to the remainder as explained on
                % page 145
                temp = remainder + extended_message;
                code(i,:) = temp;
            end
            
            
						
        end
        
        
		
		function [decoded, nERR] = decode(obj,code)
			% Decode Reed-Solomon codes
			% INPUT:
			% -obj: the RSCode object, defines the properties of the code
			% -code:  every row contains a GF(2^m) codeword of length obj.l+2*obj.t
			% OUTPUT:
			% -decoded: every row contains a GF(2^m) information word corresponding to decoding of the corresponding Reed-Solomon codeword
			% -nERR: column vector containing the number of corrected symbols for every codeword, -1 if error correction failed
			
            decoded = gf(zeros(size(code,1),obj.l), obj.m);
            nERR = zeros(size(code, 1),1);
            
            code_size = size(code,1);
            
            for code_index = 1:code_size
                %disp("starting code word: " + code_index);
                
                received = code(code_index,:);

                assert(size(received,2) == obj.l+2*obj.t);

                % The primitive element is 0000 0010, which is 2 in binary,
                % this element will be used as fourier kernel
                pe = gf(2, obj.m);

                % Calulating the syndrome polynomial by calculating ffFT(code),
                % this can be easily done with S_j-1 = r(alpha^j)

                    % Creating an array filled with the roots of our generator pol.
                    % Notice that we need to provide the option of m_0 being
                    % either 0 or 1
                roots = pe .^ flip(obj.m0:2*obj.t-1+obj.m0);
                    % The "polyval" fucntion evaluates the polynomial "received" in
                    % the point "alpha_power"
                    % Equivalent python: S = [polyval(received, alpha_power) for alpha_power in roots]
                S = arrayfun(@(alpha_power) polyval(received, alpha_power), roots);
                
                % If S is a zero-vector, our received vector is a codeword,
                % in that case we skip the reconstruction
                if any(S)
                    %disp("we found a mistake in the received codeword")
                    
                    % Solving the key equations to extract the error locator
                    % polynomial by applying the extended Euclidian algorithm
                    
                        % We can't efficiently use the "linear equations"-approach because
                        % in regular binary-codes, if it's wrong, it can only be one other thing.
                        % Now here, we have a non-binary code, which means we must not only find 
                        % the error position, but also the error values.
                        % This results in twice as much unknowns, which makes "linear
                        % equations"-approach very inefficient.


                            % Initialise as in (6.71)
                    lambda_old = gf(0,obj.m);
                    lambda_cur = gf(1,obj.m);

                    omega_old = [gf(1, obj.m) gf(zeros(1,2*obj.t), obj.m)];
                    omega_cur = S;

                            % Run Euclids algorithm as in (6.72) and (6.73)
                            % Stop when deg(OMEGA(z)) < t, determine "nu"
                    nu = 0;
                    while ~(length(omega_cur) <= (obj.t))
                        nu = nu + 1;

                        [q_cur, ~] = deconv(omega_old, omega_cur);


                        % Isn't this just the remainder from the deconv above?
                        omega_new = obj.polysub(omega_old, conv(q_cur,omega_cur), obj);
                        lambda_new = obj.polysub(lambda_old, conv(q_cur,lambda_cur), obj);

                        omega_old = omega_cur;
                        omega_cur = omega_new;

                        lambda_old = lambda_cur;
                        lambda_cur = lambda_new;

                    end

                        % Set LAMBDA(z) and OMEGA(z)
                        
                    LAMBDA = lambda_cur;
                    
                    % We now have the error locator polynomial LAMBDA
                    % but to use the key_equation it is favorable that
                    % the constant term is 1.
                    % To this end we mutiply the LAMBDA(z) polynomial
                    % with a certain power of our primitive element, we
                    % can do this without the need of compensating for
                    % this because LAMBDA(z) = 0 :D
                    % We can find this certain power of alpha because
                    % we know that alpha^obj.n = 1
                    % e.g.with alpha^0 = alpha^15 = 1:
                    %   z*alpha^3 + alpha^5 = 0
                    %   (z*alpha^3 + alpha^5)*alpha^10 = 0
                    %   z*alpha^16 + 1 = 0
                    lambda_0 = LAMBDA(size(LAMBDA,2));
                    if lambda_0 ~= 1
                       certain_exponent = obj.n - log(lambda_0);
                    else
                       certain_exponent = 0;
                    end
                    
                    LAMBDA = conv(LAMBDA, pe^certain_exponent);

                    % Use the key equation (6.61) with our syndrome and our 
                    % error locator polynomial, this should result in E 
                    E = obj.use_key_equation(LAMBDA,S,obj);

                    % Inverse Fourier transform of E, invffFT(E) = e
                    % this is defined in (6.14), we also need to flip the
                    % error polynomial because the positioning is not what
                    % we would expect on paper.
                    
                    roots_inv = pe .^ (-1*(0:obj.n));
                    e = (1/mod(obj.n,2))*arrayfun(@(alpha_power) polyval(E, alpha_power), roots_inv);
                    e = flip(e(1:obj.l+2*obj.t));
                       
                    % We reconstruct the entire codeword, including the
                    % parity symbols
                    code_reconstructed = obj.polysub(received, e, obj);
                    
                    % Now it's possible that there were more errors than we
                    % could ever correct, we now check if our reconstructed
                    % codeword is actually a codeword by dividing it by our
                    % generator polynomial and expecting a remainder of 
                    % zero, if this is not the case, we return -1 as nERR
                    % because the reconstruction failed.
                    [~, r] = deconv(code_reconstructed, obj.g);
                    if (any(r))
                        nu = -1;
                    end
                    
                else
                    % disp("we did not find a mistake in the received codeword")
                    code_reconstructed = received;
                    nu = 0;
                end
                
                % As stated in the syllabus, nu is the amount of errors we
                % (tried to)correct
                nERR(code_index,1) = nu;
                
                % Now we need to actually decode the message, since we have
                % a systematic code, this is fairly easy
                % The first obj.l indices is the codeword, the remaining
                % elements are the parity symbols
                info_reconstructed = code_reconstructed(1, 1:obj.l);
                decoded(code_index,:) = info_reconstructed;
            
            end
		
        end
		
    end
        
methods(Static)
    
        %function for 5.2 question 1, checking minimal polynomial for alpha
        %(which should give the primitive polynomial)
        function question_1()
            i = 2;
            %alpha is represented as 2 in decimal, in gf of power 8
            alpha = gf(2,8);
            %the first factor of the minimal polynomial is (x-alpha)
            poly = [gf(1,8) -alpha];
            %multiply the polynomial by (x-alpha^2)*(x-alpha^4),... while
            %the power is smaller than 2^8 (so smaller than 256)
            while i < 256
                poly = conv(poly, [gf(1,8) -(alpha^i)]);
                i = 2*i;
            end
            %display the polynomial (which indeed is equal to the primitive
            %polynomial)
            disp(poly);

        end
        
            
        
        %function for 5.2 question 2, constructing the generator polyn.
        function question_2()
            rs = RSCode(8,2,10,0);
            disp(rs.g)
            
        end
        
        %function for 5.2 question 3, constructing the check polyn.
        %In the 5th exercise on polynomial codes, the exercise mentions to 
        %define the check polynomial as (x^n - 1)/g(x).
        %I'm assuming I can define h(x) = 
        %   (x-alpha^4)(x-alpha^5)(x-alpha^6)...(x-alpha^254) 
        %as the check polynomial
        function question_3()
            
            m = 8;
            m0 = 0;
            
            %The primitive element is 0000 0010, which is 2 in binary 
            pe = gf(2, m);
            
            %The first factor of the check polynomial is (x - alpha^4),
            check = [gf(1,m) pe^m0];
            
            %we now iterativly add the other factors
            for i = 5:254
                check = conv(check, [gf(1, m) pe^i]);
            end
			
            disp("first 5 elements");
            disp(check(1:5));
            disp("last 5 elements");
            disp(check(end-5+1:end));
            
        end
            
		function generator = makeGenerator(m, t, m0)
			% Generate the Reed-Solomon generator polynomial with error correcting capability t over GF(2^m)
			% INPUT:
			% -m: positive integer
			% -t: error correction capability of the Reed-Solomon code, positive integer > 1
			% -m0: determines first root of generator polynomial, positve integer >= 0
			% OUTPUT:
			% -generator: rowvector containing the GF(2^m) generator polynomial coefficients in order of descending powers
			
            
            %Lets put a quick introduction to the "gf() function" here:
            
            %gf(3, 4) represents the third element of the galois field with
            %size 2^4, so it will output the value "3" but it denotes the
            %element 11 in binary, which represents the polynomial 1 + x
            
            %An example of addition on these objects:
            %gf(1,4) + gf(3,4) = "gf(2,4)"
            %0001 + 0011 = 0010 in GF(2)
            
            %An example of multiplication can be done with the
            %gf(3, 4)*gf(3, 4) = gf(5,4)
            %0011 * 0011 = 0101 in GF(2)
            %(x+1)*(x+1) = x^2 + 1 in GF(2)
            
            %An example of a polynomial with elements of the field GF(256)
            %The position in the array denotes the power of x, a
            %convolution of coefficient of two polynomials results in the
            %coefficients of the product of these two polynomials
            %
            %conv([gf(3, 4) gf(1, 4)], [gf(3, 4) gf(2, 4) gf(0, 4)]) 
            %      = [gf(5, 4) gf(5,4) gf(2,4) gf(0,4)]
            %
            %(3*x + 1)*(3*x^2 + 2*x) = 
            %   (gf(3, 4)*gf(3, 4))*x^3 + ((gf(3,4)*gf(2,4) + gf(1,4)*gf(3,4))) 
            %       + (gf(1,4)*gf(2,4))*x 
            %                        =
            %   (gf(5,4)*x^3 + gf(5,4)*x^2 + gf(2,4)*x + gf(0,4))
            
            %The primitive element is 0000 0010, which is 2 in binary 
            pe = gf(2, m);
            
            %The formula for the generator polyn. is given in 6.47(p154)
            %The position in the array denotes the power of x
            
            %The first factor of the generator polynomial is (x - alpha^m0),
            %So this line means generator = gf(1, 4)*x + pe^m0
            gen = [gf(1,m) pe^m0];
            
            %we now iterativly add the other factors
            for i = m0+1:(m0+2*t-1)
                gen = conv(gen, [gf(1, m) pe^i]);
            end
			
            generator = gen;
        end
        
        %function to test the encoding
        function test_encoding()
            m0 = 1; 
            % Construct the RSCode object and a message
			rs = RSCode(8,5,10,m0);
            msg = gf(randi([0,2^8-1],5,10),8);
            
            % information length and code length                                     
            %(used in the matlab rsenc function)
            code_length = rs.l + 2*rs.t ;
            information_length = rs.l ;

            % encode using the build-in matlab function
            enc1 = rsenc(msg,code_length, information_length);
            
            % encode using our own encode function
            enc2 = rs.encode(msg);      
            
            % if they are equal, display "SUCCEEDED"
            if(eq(enc1, enc2))
                disp('SUCCEEDED');
            end
        end
        
        function c = polyadd( a, b, obj)
            % addaptation of https://stackoverflow.com/questions/10516083/summing-two-polynomials
            assert( isrow(a))
            assert( isrow(b))

            maxL = max( [length(a), length(b)]);
            a = [ gf(zeros( 1, maxL - length(a)), obj.m), a];
            b = [ gf(zeros( 1, maxL - length(b)), obj.m), b];
            
            c = a + b ;
            
            index = find(c ~= 0, 1, 'first');
            c = c(index : end);
            return;
        end
        
        function c = polysub( a, b, obj)
            % addaptation of https://stackoverflow.com/questions/10516083/summing-two-polynomials
            assert( isrow(a))
            assert( isrow(b))

            maxL = max( [length(a), length(b)]);
            a = [ gf(zeros( 1, maxL - length(a)), obj.m), a];
            b = [ gf(zeros( 1, maxL - length(b)), obj.m), b];
            
            c = a - b ;
            
            index = find(c ~= 0, 1, 'first');
            c = c(index : end);
            return;
        end
        
        function E = use_key_equation(LAMBDA,S,obj)
            
            % S is ordered from large degree to small, I do it flipped so I
            % flip it (I flip E back at the end so everything is correct)
            S_flipped = flip(S);
            
            % create the error vector
            E_flipped = gf(zeros(1,(obj.n)), obj.m);
            
            %we know that E(j) = S(j-1), so j (or i) has to start at 2 (so
            %j-1 is 1 and thus first index in matlab)
            
            for i=obj.m0 + 1:size(S,2)+obj.m0
                E_flipped(i) = S_flipped(i-obj.m0);
            end

            % (so E(1)) has to be calculated a bit different, should be
            % same principle but I'm not sure about this
            x = size(LAMBDA,2) - 1;

            % this is which Ej's we still look for with formula 6.61
            % these are the ones going from where we didn't have syndromes
            % anymore (so size S +2 )
            codeword_length = obj.n;
            % recursevely find these remaining E(j)
            for i=(size(S,2)+ obj.m0 + 1):codeword_length
                E_cur = 0 ;
                % use formule 6.61
                for j=1:(x)
                    E_cur = E_cur + LAMBDA(j)*E_flipped(i - (1 + x-j));
                end
                % add this to the E (the flipped one)
                E_flipped(i) = E_cur;
            end	
            
            for i=obj.m0:-1:1
                E_cur = 0;
                for j=2:x+1
                    E_cur = E_cur + LAMBDA(j)*E_flipped(i+j-1);
                end
                E_flipped(i) = E_cur/LAMBDA(1);
            end
            % on correct GF, is this correct to do this way?
            E_flipped = gf(E_flipped, obj.m);
            
            % flip back since we need the one where left is largets degree
            E = flip(E_flipped);
            return
        end

         
		function test_decoding()
			% Test the Matlab code of this class
			
            % Also test with other values of m0, it's 0 in the CD specs!
			m0 = 0; 
			
            % Construct the RSCode object
			rs = RSCode(8,5,10,m0);
			
            % Generate a fixed message 1,2,3,....
            msg = 1:10;
			msg = gf(msg,8);
            
            disp("msg=")
            disp(msg)
            
			code = rs.encode(msg); % Encode this message
            
            disp("code.x=")
            disp(code.x)
			
            % Introduce errors
			code(1,8) = 0;
            
            disp("code.x with errors=")
            disp(code.x)
			
			[decoded, ~] = rs.decode(code); % Decode
            
            disp("decoded=")
            disp(decoded)
            
			assert(all(all(decoded(1,:) == msg(1,:))))
			
        end
        
        function test()
            % Test the Matlab code of this class
            
            m0 = 1; % Also test with other values of m0!
            
            rs = RSCode(8,5,10,m0); % Construct the RSCode object
            
            msg = gf(randi([0,2^8-1],5,10),8); % Generate a random message of 5 information words
            
            code = rs.encode(msg); % Encode this message
            
            % Introduce errors
            code(2,[3 18]) = code(2,[5 18])+1;
            code(3,8) = 0;
            code(4,[4 2 19 20 6]) = gf(randi([0,2^8-1],1,5),8);
            code(5,[4 2 19 20 6 13]) = gf(randi([0,2^8-1],1,6),8);
            
            
            [decoded,nERR] = rs.decode(code); % Decode
            
            nERR
            
            assert(all(all(decoded(1:4,:) == msg(1:4,:))))
            
        end
		
	end
	
end