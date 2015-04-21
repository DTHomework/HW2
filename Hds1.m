function [ gi ] = Hds1( wi)
u = 0
a = [ 1 -4.4153 8.6283 -9.4592 6.1051 -1.3542 -3.3622 7.239 -7.9361 5.1221 -1.8401 2.8706e-1 ];
b = [1.3651e-4 8.1905e-4 2.0476e-3 2.7302e-3 2.0476e-3 9.0939e-4 6.7852e-4 1.355e-3 1.8067e-3 1.355e-3 5.3726e-4 6.1818e-5 -7.1294e-5 -9.5058e-5 -7.1294e-5 -2.5505e-5 1.3321e-5 4.5186e-5 6.0248e-5 4.5186e-5 1.8074e-5 3.0124e-6];
[h, ] = impz(b, a);
h_energy = sum(h.^2);
b = b/sqrt(h_energy);     %Normalize the filter to have energy 1
gi = 0;%filter(b, a, wi);
[H,w] = freqz(b,a,10000, 'whole');
% fvtool(b, a);
figure
plot(linspace(0, 0.9999, 10000), 10*log10(H));
axis([0 1 -10 10]);

% wi = [ wi zeros(1, length(H) - length(wi))];
% W = fft(wi);
% Gi = W.'.*H;
% gi = ifft(Gi);


gi = filter(b, a, wi);
length(gi)
% gi = gi(200:210 ); %Transitorio stimato 150, controllare modi del sistema per stima piu accurata
% figure
% plot(0:length(gi)-1,(abs(gi)));

end

