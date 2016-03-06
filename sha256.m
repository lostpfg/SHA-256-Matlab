function out = sha256( msg )
    % Initial Hash Values(8 constant 32-bit words).(§5.3.3)
    default_hash = [
                    '6a09e667';
                    'bb67ae85';
                    '3c6ef372';
                    'a54ff53a';
                    '510e527f';
                    '9b05688c';
                    '1f83d9ab';
                    '5be0cd19'
                    ];
    % Constant value array (64 constant 32-bit words) to be used for the iteration t of the hash computation.(§4.2.2)
    K = [
        '428a2f98'; '71374491'; 'b5c0fbcf'; 'e9b5dba5'; 
        '3956c25b'; '59f111f1'; '923f82a4'; 'ab1c5ed5';
        'd807aa98'; '12835b01'; '243185be'; '550c7dc3'; 
        '72be5d74'; '80deb1fe'; '9bdc06a7'; 'c19bf174';
        'e49b69c1'; 'efbe4786'; '0fc19dc6'; '240ca1cc';
        '2de92c6f'; '4a7484aa'; '5cb0a9dc'; '76f988da';
        '983e5152'; 'a831c66d'; 'b00327c8'; 'bf597fc7';
        'c6e00bf3'; 'd5a79147'; '06ca6351'; '14292967';
        '27b70a85'; '2e1b2138'; '4d2c6dfc'; '53380d13'; 
        '650a7354'; '766a0abb'; '81c2c92e'; '92722c85';
        'a2bfe8a1'; 'a81a664b'; 'c24b8b70'; 'c76c51a3'; 
        'd192e819'; 'd6990624'; 'f40e3585'; '106aa070';
        '19a4c116'; '1e376c08'; '2748774c'; '34b0bcb5'; 
        '391c0cb3'; '4ed8aa4a'; '5b9cca4f'; '682e6ff3';
        '748f82ee'; '78a5636f'; '84c87814'; '8cc70208'; 
        '90befffa'; 'a4506ceb'; 'bef9a3f7'; 'c67178f2'
    ];
    % First padd the input message to be a multiple of 512(bits).(§5)
    [padded_msg,padded_len] = padder( msg );
    % Split padded message to N (512-bit) blocks.(§6)
    [M,total_blocks] = split2block( padded_msg,padded_len );
    W = zeros( 64, 32 );
    % Main SHA-256 compuation process.(§6.2.2)
    for i = 1:total_blocks % For every block M(i).
        if i == 1 % Load initial hash values at first iteration.
            h0 = hexToBinaryVector( default_hash( 1 , : ) , 32 );
            h1 = hexToBinaryVector( default_hash( 2 , : ) , 32 );
            h2 = hexToBinaryVector( default_hash( 3 , : ) , 32 ); 
            h3 = hexToBinaryVector( default_hash( 4 , : ) , 32 );
            h4 = hexToBinaryVector( default_hash( 5 , : ) , 32 );
            h5 = hexToBinaryVector( default_hash( 6 , : ) , 32 );
            h6 = hexToBinaryVector( default_hash( 7 , : ) , 32 );
            h7 = hexToBinaryVector( default_hash( 8 , : ) , 32 );
        end
        % Step 1 - Prepare the message schedule.
        for j = 1:64
            if j >= 1 && j <= 16
                W( j, 1:32 ) = M( i, 32*(j-1)+1:j*32 );
            else
                W( j, 1:32 ) = mod32add( sigma1( W( j-2, : ) ), W( j-7, : ) , sigma0( W( j-15, : ) ), W( j-16, : ) );
            end
        end
        % Step 2 - Initialize the eight working variables, a, b, c, d, e, f, g,
        % and h, with the (i-1)st hash value.
        a = h0;
        b = h1;
        c = h2;
        d = h3;
        e = h4;
        f = h5;
        g = h6;
        h = h7;
        % For t=0 to 63.
        for t = 1:64
            T1 = mod32add( h, capSigma1(e),ch( e, f, g ),hexToBinaryVector( K( t,: ) , 32 ), W( t,:  ) );
            T2 = mod32add( capSigma0(a), maj( a, b ,c )  );
            h = g;
            g = f;
            f = e;
            e = mod32add( d, T1 );
            d = c;
            c = b;
            b = a;
            a = mod32add( T1, T2 );
        end
        % Step 4 - Compute the ith intermediate hash value H(i).
        h0 =  mod32add( a, h0 );
        h1 =  mod32add( b, h1 );
        h2 =  mod32add( c, h2 );
        h3 =  mod32add( d, h3 );
        h4 =  mod32add( e, h4 );
        h5 =  mod32add( f, h5 );
        h6 =  mod32add( g, h6 );
        h7 =  mod32add( h, h7 );    
    end
    % Final Step - After hash process the resulting 256-bit message digest
    % of the message, M, is:
    out = calcout( h0, h1, h2, h3, h4, h5 ,h6, h7 );
end

function [out,len] = padder( msg )
    % Function padder : Padds the input message.(§5.1.1)
    padded = []; % Initialize output.
    l = length(msg)*8; % Length of the input message in dec.
    for i = 1:length(msg) % First append message body.
        padded = strcat(padded,dec2bin(msg(i),8));
    end
    padded( end + 1 ) = '1'; % Append bit '1' at the end of message body.
    % Calculate number of zeros to be added at the padded message.
    k = mod( 447 - l , 512 );
    padded( end + 1 : end + k ) = '0';  % Append k bits '0' at the end of message body.
    % Append the length of the input message (in 64-bits).
    padded( end + 1: end + 64 ) = reshape( dec2bin( l, 64 ), 1, [] );
    out = logical(padded(:)'-'0'); % Convery to logical array.
    len = length( padded ); % Return also length of the padded message.
end

function [out,total_blocks] = split2block( padded_msg,padded_len )
    % Function split2block : Splits the padded message to N 512-bit blocks M(N).(§5.2.1)
    total_blocks = padded_len/512; % Calculate total number of blocks.
    out = zeros( total_blocks, 512 );
    for i = 1:total_blocks % Split per 512 bits (Big-Endianess).
        out( i, 1:512 ) = padded_msg( (i-1)*512 + 1:i*512 );
    end
end

function out = fix2mod( x )
    % Function fix2mod : Converts the input logical word to binary.
    out = num2str( x );
    out(isspace(out)) = '';
    out = bin2dec(out);
end

function out = mod32add( varargin )
    % Function mod32add : Performs addition modulo 32.(§3.2.1)
    narginchk(1, inf)
    out = 0; % initialise return arguments
    nums = zeros( 1, length(varargin) ); % Set array of input variables.
    for i=1:length(varargin) % Calculate addition
        nums(i) = fix2mod(varargin{i});
        out = out + nums(i);
    end
    % Perform modulo 32 operation.
    out1 = dec2bin( mod(( out ),2^32),32 ) ;
    out = logical( out1(:)'-'0' ); % Cast output to logical array.
end

function out = rotr( word, pos )
    % Function rotr : Performs ROTR (Cirular right shift) operation.(§3.2.1)
    out = zeros( 1, length( word ) );
    out( pos + 1 : end ) = word( 1 : end - pos + 0 ); 
    out( 1 : pos ) = word( end - pos + 1 : end ); 
end

function out = shr( word, pos )
    % Function rotr : Performs SHR (Right shift) operation.(§3.2.1)
    out = zeros( 1, length( word ) );
    out( 1 + pos : end ) = word( 1 : end - pos ); 
end

function out = maj( x, y, z )
    % Function maj : Performs MAJ operation.(§4.1.2)
    out = bitxor( bitxor( x & y, x & z ) , y & z );
end

function out = ch( x, y, z )
    % Function maj : Performs CH operation.(§4.1.2)
    out = bitxor( x & y ,~x & z );
end

function out = capSigma0( x )
    % Function Ó0 : Performs Ó0 operation.(§4.1.2)
    out = bitxor( bitxor( rotr( x, 2 ), rotr( x, 13 ) ) , rotr( x, 22 ) );
end

function out = capSigma1( word )
    % Function Ó1 : Performs Ó0 operation.(§4.1.2)
    out = bitxor( bitxor( rotr( word, 6 ), rotr( word, 11 ) ) , rotr( word, 25 ) );
end

function out = sigma0( word )
    % Function ó0 : Performs Ó0 operation.(§4.1.2)
    out = bitxor( bitxor( rotr( word, 7 ), rotr( word, 18 ) ) , shr( word, 3 ) );
end

function out = sigma1( word )
    % Function ó1 : Performs Ó0 operation.(§4.1.2)
    out = bitxor( rotr( word, 17 ), rotr( word, 19 ) );
    out = bitxor( out, shr( word, 10 ) );
end

function out = calcout( h0, h1, h2, h3, h4, h5 ,h6, h7 )
    % Function calcout : Merges the output string (256-bit).
    t1 = binaryVectorToHex(h0);
    t2 = binaryVectorToHex(h1);
    t3 = binaryVectorToHex(h2);
    t4 = binaryVectorToHex(h3);
    t5 = binaryVectorToHex(h4);
    t6 = binaryVectorToHex(h5);
    t7 = binaryVectorToHex(h6);
    t8 = binaryVectorToHex(h7);
    out = strcat( t1, t2, t3, t4, t5, t6, t7, t8 );
end
