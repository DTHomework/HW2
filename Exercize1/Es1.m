clear all;
close all;

%aggiungiamo cose a caso
filename = 'data.mat';
importfile(filename);

K = length(z); %length of the signal
f = linspace(0, 1, K); %space of the frequencies
f2 = linspace(0, 1, 256);

% %Plot the signal
% figure
% plot(1:1000, z);
% 
% hold on 
% 
% plot(1:1000, imag(z), 'r');
% 
% %Autocorrelation of the signal
% %R_est = zeros([K + 1, 1]);
% R_est= zeros(K/5 + 1, 1);
% 
% for j = 0:1:(K/5)
%     for i = j:1:(K-1)
%     R_est(j+1, 1) = R_est(j+1, 1) + (1/(K - j))*z(i+1)*conj(z(i-j+1));
%     
%     end
% end
% figure
% plot(1:length(R_est), abs(R_est));
% title('Autocorrelation of the signal');
% grid on
%Calculate the DFT
Z = fft(z);
% 
% figure
% subplot(2,1,1);
% plot(f, 20*log10(abs(Z)));  %amplitude
% 
% subplot(2,1,2);
% plot(f, (angle(Z)));    %phase

%CORRELOGRAM
L = K/5; %number of samples to use in the computation of the correlogram
corr = [];
[a, b, R_toep] = armodel(L, K, z);

for i = 1 : L
    corr( i,1 ) = R_toep( 1, L + 1 - i );
end
for i = 1 : L - 1
    corr( i + L,1) = R_toep( i + 1 , 1 );
end
corr = corr .* hamming(length(corr));
COR = zeros([2*L - 1, 1]); 
for l = 1 : 2*L - 1 %DFT fatta a mano
    for n = 1 : 2*L - 1
        COR( l,1 ) = COR( l,1 ) + corr( n,1 )*exp(-2i*pi*(l - 1)*( n - L )/(2*L - 1));
    end
end

% COR = fft( corr ); 

figure
%subplot(2, 1, 1)
%plot( linspace(0, 1, length(COR)) , abs( COR ) );

%subplot(2, 1, 2)
plot( linspace(0, 1, length(COR)) , 10*log10( abs( COR ) ) );
xlabel('Frequency f');
ylabel('PSD(f) dB');
title('PSD Correlogram');
%ARmodel

%fvtool
[a, sigma, matrix] = armodel(3, K, z);
% fvtool(1, a)

[H, ] = freqz(1, a);
% fvtool(1,a);
psd_ar = sigma*(abs(H).^2);


% figure
% subplot(2, 1, 1)
% plot( linspace(0, 1, length(psd_ar)) , abs( psd_ar ) );
% title('PSD AR');
% subplot(2, 1, 2)
% plot( linspace(0, 1, length(psd_ar)) , 10*log( abs( psd_ar ) ) );


%PERIODOGRAM
PPer = (1/K)*(abs(Z)).^2;
PPer_plot = PPer;
figure

%subplot(2,1,1);
%plot(f, PPer);
%ylabel('PPer(f)');


%subplot(2,1,2);
plot(f, 10*log10(PPer));
title('PSD Periodogram');
xlabel('Frequency f');
ylabel('PSD(f) dB');


%ARMODEL
[a, sigma, matrix] = armodel(3, 1000, z);
% fvtool(1, a)

% H = 1/fft(a);
a = [1; a];
[H, ] = freqz(1, a, 'whole');

psd_ar = sigma*(abs(H).^2);

figure
plot( linspace(0,1,length(H)) , 10*log10( abs( psd_ar ) ), 'r' );
title('PSD AR model')
xlabel('Frequency f');
ylabel('PSD(f)dB');

%WELCH PERIODOGRAM
D = K/5; %extracted samples
S = D/2; %overlap samples
Ns = floor( (K-D)/(D-S) +1 ); % number of subsequences
w = hamming( D ); %window
Mw = 0;
for i = 1: D
    Mw = Mw + 1/D*(w(i))*w(i);
end
Pwe = zeros( D, 1);

for s = 0 : Ns-1
    xs = [];
    PPer = [];
    for k = 0 : D-1
        xs( k+1, 1 ) = w( k+1, 1 )*z( k+1 + s*( D - S ));
    end
    Xs = fft( xs );
    PPers = (1/(D*Mw))*(abs(Xs)).^2;
    Pwe = Pwe + 1/Ns*PPers;
end

figure
plot( linspace( 0, 1, D), 10*log10(Pwe));
xlabel('Frequency f');
ylabel('PSD(f) dB');
title('PSD Welch Periodogram')


%Overall plot of PSDs
figure 
%plot( linspace(0, 1, length(COR)) , 10*log10( abs( COR ) ) );
plot(f , 10*log10( abs( resample(COR, K, length(COR))) ), 'r', 'linewidth', 0.25  );
hold on
plot(f, 10*log10(PPer_plot), 'c', 'linewidth', 0.25 );
hold on
plot( linspace(0,1,length(H)) , 10*log10( abs( psd_ar ) ), 'b', 'linewidth', 2 );
hold on
plot( linspace( 0, 1, D), 10*log10(Pwe), 'm', 'linewidth', 2.5);

legend('Correlogram', 'Periodogram', 'AR model', 'Welch', 'Location', 'southeast');


%WHITENING FILTER
white = filter(a, 1, z);
WHITE = fft(white);
figure
% subplot(2, 1, 1)
% plot(f, 10*log10(abs(white)));
% axis([0 1 -10 10]);

title('Whitening');
%subplot(2, 1, 2)
plot(f, 10*log10(abs(WHITE)));
axis([0 1 15 28]);
xlabel('fT_c');
ylabel('Amplitude(dB)');




%Periodogram2
% PPer2 = periodogram(z);
% PPer2 = PPer2(1:K);
% figure
% 
% subplot(2,1,1);
% plot(f, PPer2);
% ylabel('PPer2(f)');
% title('PSD Periodogram2');
% 
% subplot(2,1,2);
% plot(f, 10*log10(PPer2));
% ylabel('PPer2(f)dB');

% %Welch Periodogram
% PWel = pwelch(z,150);
% %PWel = PWel(1:K);
% figure
% 
% subplot(2,1,1);
% plot(f2, PWel);
% ylabel('PWelch(f)');
% title('PSD Welch Periodogram');
% 
% subplot(2,1,2);
% plot(f2, 10*log10(PWel));
% ylabel('PWel(f)dB');
% 
% %CORRELOGRAM
% 
% corr = [];
% [a, b, R_toep] = armodel(1000, 1000, z);
% 
% for i = 1 : 1000
%     corr( i,1 ) = R_toep( 1, 1000 + 1 - i );
% end
% for i = 1 : 999
%     corr( i + 1000,1) = R_toep( i + 1 , 1 );
% end
% 
% COR = zeros([1999, 1]); 
% for l = 1 : 1999 %DFT fatta a mano
%     for n = 1 : 1999
%         COR( l,1 ) = COR( l,1 ) + corr( n,1 )*exp(-2i*pi*(l - 1)*( n - 1000 )/(1999));
%     end
% end
% 
% % COR = fft( corr ); 
% 
% figure;
% title('PSD Correlogram');
% subplot(2, 1, 1)
% plot( linspace(0, 1, length(COR)) , abs( COR ) );
% subplot(2, 1, 2)
% plot( linspace(0, 1, length(COR)) , 10*log10( abs( COR ) ) );
% 
% % %FILTERING
% % b = [0.5, -0.5*exp(1i*pi*0.084)];
% % z1 = filter(b, 1, z);
% % fvtool(b, 1, 256);
% 
% load('C:\Users\Silvia\Desktop\Silvia\Documenti personali\Università\Digital Transmission\coeff1.mat')
% z1 = filter(Num,1,z);
% 
% %Periodogram
% Z1 = fft(z1);
% PPer = (1/K)*(abs(Z1)).^2;
% 
% figure
% 
% subplot(2,1,1);
% plot(f, PPer);
% ylabel('PPer(f)');
% title('PSD Periodogram');
% 
% subplot(2,1,2);
% plot(f, 10*log10(PPer));
% ylabel('PPer(f)dB');






