function [ ml ] = maximalLength10()
a = ones( 1, 11 );
a(1) = 0;
ml(1) = a(1);
for l = 2 : 2^10
    b = xor( a(8), a(11) )
    for i = 11 : -1 : 2
        a( i ) = a( i - 1 );
    end
    a(1) = b;
    ml( l ) = a(1);
    a
end
end
