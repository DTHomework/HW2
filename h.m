function [ d1 H ] = h( ml )
T = 4;      % Sampling time before interpolation
tau_rms = 0.3*T;       
tau = linspace(0, 5, 100);
PDP_continuous = exp(-tau/tau_rms)./tau_rms;
% figure
% plot(tau, PDP_continuous);
% title('Continuous PDP (T_{sample} = T)');

%I'll now sample the "continuous" PDP
Tc = (0.25*T);      %New sampling time
tau = 0:Tc:5;
samples = 3;
PDP_sampled = exp(-tau/tau_rms)./tau_rms;
%PDP_sampled = [PDP_sampled(1:samples) zeros(1, length(PDP_sampled )-samples)] ;
% figure
% stem(tau, PDP_sampled, 'm');
% title('Sampled PDP (T_{sample} = Tc)');


%normalization of the PDP
K = 3; %K in dB
K = 10^(K/10); % K in linear

Md = sum(PDP_sampled);

C = sqrt(K/(K+1));

norm = Md/( 1 - C^2 ); % MdNorm = 1 - c^2

MdNorm = Md/norm; %PDP normalized

sum( MdNorm ) + C^2;
PDP_sampled = PDP_sampled/norm;
% figure
% stem(tau, 10*log10(PDP_sampled), 'm');
% title('Sampled Normalized PDP (T_{sample} = Tc)');

%POLYPHASE METHOD

h0 = sqrt(PDP_sampled(1)+C);
h1 = sqrt(PDP_sampled(2));
h2 = sqrt(PDP_sampled(3));
H = [h0 h1 h2 0];

%White noise
% inp_length = 1000;
% w = wgn(1, inp_length, 0);
%Compute the power of the noise from the one of the signal
Mx = sum(ml.^2);
sigmawsq = Mx/40;
%Polyphase
w0 = conv(ml, h0);
w1 = conv(ml, h1);
w2 = conv(ml, h2);
w3 = conv(ml, 0);
W = {w0 w1 w2 w3};
x = zeros(1, length(ml)*4);

for i = 0:length(ml)*4 - 1
    j = mod(i, 4) + 1;
    x(i+1) = W{j}(fix(i/4) + 1);
end



%SYSTEM NOISE
w_sys = wgn(1, length(ml)*4, 10*log10(sigmawsq), 'complex');   %This one is complex
d = x+(w_sys);
for i = 1 : 4
    d1( :, i ) = d( i : 4 : end );
end

end