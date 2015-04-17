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
subplot(2, 1, 1);
stem(tau, 10*log10(PDP_sampled), 'm');
axis([0 5 -10 6]);
title('Sampled PDP (T_{sample} = Tc)');
subplot(2, 1, 2);
stem(tau, PDP_sampled, 'm');

%Now we build a proper white noise with 0 mean and power = 1
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
