function [ d1 H sigmawsq ] = h( ml )
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

fd = 5*10.^(-3)/T;
inp_length = 10000;

gtilda = cell(samples, 1);
giInt = cell(samples, 1);

for i = 1:samples
%White noise
w_i = wgn(1, inp_length, 0, 'complex');
%NarrowBand Filter
gtilda{i} = Hds1(w_i);

%Cubic Interpolator
x = 0:length(gtilda{i})-1;

xx = 0: (1/80) : length(gtilda{i})-1;

%         figure
%     plot(x, (10*log10(abs(gtilda{i}))));
%     title('Absolute value of gtilda (dB)');
%     xlabel('tau');
%     ylabel('|g_i| (dB)');

giInt{i} = interp1( x, gtilda{i}, xx, 'spline');

end

g_mean = {mean(giInt{1}(:)) mean(giInt{2}(:)) mean(giInt{3}(:))};
g_var = {var(giInt{1}(:)) var(giInt{2}(:)) var(giInt{3}(:))};
% figure
% plot(x,gi,'o',xx,giInt)

%Sigma_i
sigma_i = zeros(samples, 1);

for i = 1:samples
   sigma_i(i) = sqrt(PDP_sampled(i)); 
   
end

%g_i
g_i = cell(samples, 1);

for i = 1:samples
    
    g_i{i} = sigma_i(i) * giInt{i};
    
    if i == 1
        g_i{i} = g_i{i} + C;
    end
end
%POLYPHASE METHOD


%White noise
% inp_length = 1000;
% w = wgn(1, inp_length, 0);
%Compute the power of the noise from the one of the signal
Mx = sum(ml.^2)/length(ml);
sigmawsq = Mx/40;
%Polyphase

for t = 1 : 3
    gi(t,:) = g_i{t}( t : 4 :(length(ml)-1)*4+t);
    w(t, :) = ml.*gi(t,:);
end
w0 = w(1, :); w1 = w(2, :); w2 = w(3, :); w3 = zeros(1, length(w2));
W = {w0 w1 w2 w3};
x = zeros(1, length(ml)*4);
W;
% for i = 0:length(ml)*4 - 1
%     j = mod(i, 4) + 1;
%     x(i+1) = W(j,:)(fix(i/4) + 1);
% end

for i = 1 : 4
    x( i : 4 : length(ml)*4 ) = W{i}( 1 : length( W{i} ));
end
x;


%SYSTEM NOISE
w_sys = wgn(1, length(ml)*4, 10*log10(sigmawsq), 'complex');   %This one is complex
d = x+(w_sys);
for i = 1 : 4
    d1( :, i ) = d( i : 4 : end );
end
H = [ mean(gi(1,:)) mean(gi(2,:)) mean(gi(3,:)) mean(0)];
%[var(gi(1,:)) var(gi(2,:)) var(gi(3,:)) 0];
%H = abs(H);
end