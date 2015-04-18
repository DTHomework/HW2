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

%Here we build a proper white noise
w_i = wgn(1, 1000, 1, 'complex');
lin = linspace(0, 0.999, 1000);
fd = 5*10.^(-3);
f = linspace(-fd, fd, 1001);
sqrt_D = zeros(1, 1001);
for i = 1:1001
sqrt_D(i) = sqrt((1/(pi*fd*sqrt(1-(f(i)/fd).^2))).*(abs(f(i)) <= fd));
end
figure
plot(f, 10*log10(sqrt_D));
title('DOPPLER SPECTRUM (dB), Classical (Jake) model');
axis([-fd fd 8 15]);

K = 3; %K in dB
K = 10^(K/10); % K in linear

Md = sum(PDP_sampled) 

C = sqrt(K/(K+1));

norm = Md/( 1 - C^2 ); % MdNorm = 1 - c^2

MdNorm = Md/norm; %PDP normalized

sum( MdNorm )+ C^2 
figure
stem(tau, PDP_sampled/norm, 'm');
title('Sampled PDP (T_{sample} = Tc)');