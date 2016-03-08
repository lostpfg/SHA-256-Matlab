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
    H = zeros( 8, 32 );
    % Main SHA-256 compuation process.(§6.2.2)
    for j = 1:8 % Load initial hash values at first iteration.
        H(j,:) = hexToBinaryVector( default_hash( j , : ) , 32 );
    end
    for i = 1:total_blocks % For every block M(i).
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
        a = H(1,:);
        b = H(2,:);
        c = H(3,:);
        d = H(4,:);
        e = H(5,:);
        f = H(6,:);
        g = H(7,:);
        h = H(8,:);
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
        H(1,:) =  mod32add( a, H(1,:)  );
        H(2,:) =  mod32add( b, H(2,:)  );
        H(3,:) =  mod32add( c, H(3,:)  );
        H(4,:) =  mod32add( d, H(4,:)  );
        H(5,:) =  mod32add( e, H(5,:)  );
        H(6,:) =  mod32add( f, H(6,:)  );
        H(7,:) =  mod32add( g, H(7,:)  );
        H(8,:) =  mod32add( h, H(8,:)  );    

    end
    % Final Step - After hash process the resulting 256-bit message digest
    % of the message, M, is:
    out = binaryVectorToHex( horzcat( H(1,:), H(2,:), H(3,:), H(4,:), H(5,:), H(6,:) ,H(7,:), H(8,:)  ) );
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
    out = 0; % initialise return arguments
    for i = 1:length( varargin ) % Calculate addition
        out = out + fix2mod(varargin{i});
    end
    % Perform modulo 32 operation.
    out = dec2bin( mod( ( out ), 2^32 ),32 ) ;
    out = logical( out(:)'-'0' ); % Cast output to logical array.
end

function out = rotr( word, pos )
    % Function rotr : Performs ROTR (Cirular right shift) operation.(§3.2.1)
    out = zeros( 1, length( word ) );
    out( pos + 1:end ) = word( 1:end - pos + 0 ); 
    out( 1:pos ) = word( end - pos + 1 : end ); 
end

function out = shr( word, pos )
    % Function rotr : Performs SHR (Right shift) operation.(§3.2.1)
    out = zeros( 1, length( word ) );
    out( 1 + pos:end ) = word( 1:end - pos ); 
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
    out = bitxor( bitxor( rotr( word, 17 ), rotr( word, 19 ) ), shr( word, 10 ) );
end
