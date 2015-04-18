clear all
close all
T = 1;      % Sampling time before interpolation
tau_rms = 0.3*T;       
tau = linspace(0, 5, 1000);
PDP_continuous = exp(-tau/tau_rms)./tau_rms;
figure
plot(tau, PDP_continuous);
title('Continuous PDP (T_{sample} = T)');

%I'll now sample the "continuous" PDP
Tc = (0.25*T);      %New sampling time
tau = 0:Tc:5;
samples = 4;
PDP_sampled = exp(-tau/tau_rms)./tau_rms;
PDP_sampled = [PDP_sampled(1:samples) zeros(1, length(PDP_sampled )-samples)] ;
figure
stem(tau, PDP_sampled, 'm');
title('Sampled PDP (T_{sample} = Tc)');



%normalization of the PDP
K = 3; %K in dB
K = 10^(K/10); % K in linear

Md = sum(PDP_sampled);

C = sqrt(K/(K+1));

norm = Md/( 1 - C^2 ); % MdNorm = 1 - c^2

MdNorm = Md/norm; %PDP normalized

sum( MdNorm ) + C^2 
figure
stem(tau, PDP_sampled/norm, 'm');
title('Sampled Normalized PDP (T_{sample} = Tc)');

%filter to create g

%Classical Doppler
%Here we build a proper white noise

fd = 5*10.^(-3);
% lin = linspace(0, 0.999, 1000);
% f = linspace(0, 1, 100001);
% sqrt_D = zeros(1, 50001);
% for i = 1:50001
%     sqrt_D(i) = sqrt((1/(pi*fd*sqrt(1-(f(i)/fd).^2))).*(abs(f(i)) < fd));%ho messo < invece di <= perche per f=fd la funzione diventa infinito
% end
% for i = 1 : 50000 %whole spectrum of Hds in [0,1]
%     sqrt_D( i + 50001 ) = conj( sqrt_D( 50001 - i ) );
% end
% figure
% plot(f, 10*log10(sqrt_D));
% title('DOPPLER SPECTRUM (dB), Classical (Jake) model');
% axis([0 1 8 15]);
% 
% 
% %Chebychev lowpass filter
% n = 10; %order of cheb low pass filter
% Rp = 1;%decibels of peak-to-peak passband ripple
% Wp =  fd;%normalized passband edge frequency
% [b,a] = cheby1(n,Rp,Wp); %returns the coefficent of the TF
% Hds2 = freqz(b,a,100001,'whole');
% figure 
% plot( 0:0.00001:1, 10*log10( abs(Hds2) ));
% title('Frequency response of Chebychev filter (dB)');
% 
% for i = 1:100001
%     Hds(i) = sqrt_D(i)*Hds2(i);
% end
% 
% figure 
% plot( 0:0.00001:1, Hds );%ci sono problemi in dB peche la maggior parte dei valori e' 0
% title('Frequency response of Hds filter (linear)');
% 
% 
% hds = ifft( sqrt_D );
% gi = conv(hds,w_i);
% figure
% plot( 0 : length(gi)-1, gi );
% title('gi');
gtilda = cell(samples, 1);
giInt = cell(samples, 1);

for i = 1:samples
%White noise
w_i = wgn(1, 1000, 1, 'complex');

%NarrowBand Filter
gtilda{i} = hds(w_i);

%Cubic Interpolator
x = 0:length(gtilda{i})-1;
xx = 0: 0.25 : length(gtilda{i})-1;
giInt{i} = spline( x, gtilda{i}, xx);

end
% figure
% plot(x,gi,'o',xx,giInt)

%Sigma_i
sigma_i = zeros(samples, 1);

for i = 1:samples
   sigma_i(i) = sqrt(PDP_sampled(i)/norm); 
end

%g_i
g_i = cell(samples, 1);

for i = 1:samples
    
    g_i{i} = sigma_i(i) * giInt{i};
    
    if i == 1
        g_i{i} = g_i{i} + C;
    end
    figure
    plot(xx, g_i{i});
   
end