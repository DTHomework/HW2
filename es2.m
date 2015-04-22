clear all
close all
N = 1;

r = 3; % ordine ML

ml =maximalLength10(r, N);
[d1 hreal] = h( ml ); % 4 colonne, una per fase

%lavoriamo con le fasi singolarmante 
N = 1;
L = 2^r-1;
for n = 1 : 4 %per ogni fase
    for i = 1 : L
        I( i, : ) = ml( (N-1)+i : -1 : i );
    end
    o = d1( N-1+1 : N-1+L , n );
    
    PHI = I'*I;
    theta = I'*o;
    hstim(n) = inv(PHI)*I'*o;
    
end
error = abs(hreal- hstim)
    

