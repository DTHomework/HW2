clear all
close all
N = 1;

r = 4; % ordine ML
for t = 1 : 1000
    ml =maximalLength10(r, N);
    [d1 hreal sigmawsq] = h( ml ); % 4 colonne, una per fase

    %lavoriamo con le fasi singolarmante 
    N = 1;
    L = 2^r-1;
    for n = 1 : 4 %per ogni fase
        for i = 1 : L
            I( i, : ) = ml( (N-1)+i : -1 : i );
        end
        o = d1( N-1+1 : N-1+L , n );

        PHI{n} = I'*I;
        theta = I'*o;
        hstim(n) = inv(PHI{n})*I'*o;

    end
    %error = abs(hreal- hstim)
    a{t} = hreal- hstim;
    a{t} = abs(a{t}).^2;
    errorReal{t} = sum(a{t});
    t
end

mer = 0;
mepr = [0 0 0 0];
for t = 1 : 1000
    mer = mer + errorReal{t};
    mepr = mepr + a{t};
end


errorPhase = { sigmawsq*inv( conj(PHI{1}))  sigmawsq*inv( conj(PHI{2}))  sigmawsq*inv( conj(PHI{3}))  sigmawsq*inv( conj(PHI{4}))} 
mepr = mepr/1000
errorTheo = errorPhase{1} + errorPhase{2} + errorPhase{3} + errorPhase{4}
mer = mer/1000



    

