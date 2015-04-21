clear all
close all
T = 4;      % Sampling time before interpolation
tau_rms = 0.3*T;       
tau = linspace(0, 5, 100);
PDP_continuous = exp(-tau/tau_rms)./tau_rms;
figure
plot(tau, PDP_continuous);
title('Continuous PDP (T_{sample} = T)');

%I'll now sample the "continuous" PDP
Tc = (0.25*T);      %New sampling time
tau = 0:Tc:5;
samples = 3;
PDP_sampled = exp(-tau/tau_rms)./tau_rms;
%PDP_sampled = [PDP_sampled(1:samples) zeros(1, length(PDP_sampled )-samples)] ;
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

sum( MdNorm ) + C^2;
PDP_sampled = PDP_sampled/norm;
figure
stem(tau, 10*log10(PDP_sampled), 'm');
title('Sampled Normalized PDP (T_{sample} = Tc)');

%POLYPHASE METHOD

h0 = PDP_sampled(1);
h1 = PDP_sampled(2);
h2 = PDP_sampled(3);

%White noise
inp_length = 1000;
w = wgn(1, inp_length, 0);
%Compute the power of the noise from the one of the signal
Mx = sum(w.^2)
sigmawsq = Mx/40;
%Polyphase
w0 = conv(w, h0);
w1 = conv(w, h1);
w2 = conv(w, h2);
w3 = conv(w, 0);
W = {w0 w1 w2 w3};
x = zeros(1, inp_length*4);

for i = 0:inp_length*4 - 1
    j = mod(i, 4) + 1;
    x(i+1) = W{j}(fix(i/4) + 1);
end

figure
plot(0:0.25:999.75, x);
axis([0 999.75 -1 1]);

%DIRECT METHOD

% w_upsampled = zeros(1, 4000);
% for i = 0:inp_length*4 - 1
%     if (mod(i, 4) == 0)
%         w_upsampled(i + 1) = w(fix(i/4) + 1);
%     end
% end
% x1 = conv(w_upsampled, [h0 h1 h2]);
% x1 = x1(1, 3:length(x1));
% 
% figure
% plot(0:0.25:999.75, x1);

%SYSTEM NOISE
w_sys = wgn(1, inp_length*4, 10*log10(sigmawsq), 'complex');   %This one is complex
d = x.*(w_sys);

figure
plot(0:0.25:999.75, d);
axis([0 999.75 -1 1]);
