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
PDP_sampled = exp(-tau/tau_rms)./tau_rms;
figure
stem(tau, PDP_sampled, 'm');
title('Sampled PDP (T_{sample} = Tc)');