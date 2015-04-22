function [ ml ] = maximalLength10(r, N)
index = r - 1;
l_table = {[1 2] [2 3] [3 4] [3 5] [5 6] [6 7] [2 3 4 8] [5 9] [7 10] [9 11] [2 10 11 12]};
a = ones( 1, r );  
for l = 1 : 2^r-1
    if(r == 8 || r == 12)
        b = xor( a(l_table{index}(1)), a(l_table{index}(2)) );
        b = xor(b, a(l_table{index}(3)));
        b = xor(b, a(l_table{index}(4)));
    else
        b = xor( a(l_table{index}(1)), a(l_table{index}(2)) );
    end
    for i = r : -1 : 2
        a( i ) = a( i - 1 );
    end
    a(1) = b;
    ml( l ) = a(1);
end
b = [];
for i = 1:length(ml)
    if ml(i) == 0
        b(i) = -1;
    else
        b(i) = 1;
    end
end
length(ml);
sum(b);
1/mean( b);

ml = [ml(1:N) ml];
length(ml)   ;
end
