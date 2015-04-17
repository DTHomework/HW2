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