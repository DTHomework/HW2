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

sum( MdNorm ) + C^2 
figure
stem(tau, 10*log10(PDP_sampled/norm), 'm');
title('Sampled Normalized PDP (T_{sample} = Tc)');

%filter to create g

%Classical Doppler
%Here we build a proper white noise

fd = 5*10.^(-3);

gtilda = cell(samples, 1);
giInt = cell(samples, 1);
inp_length = 1000;

for i = 1:samples
%White noise
w_i = wgn(1, inp_length, 0, 'complex');
%NarrowBand Filter
gtilda{i} = hds(w_i);

%Cubic Interpolator
x = 0:length(gtilda{i})-1;
xx = 0: 0.25 : length(gtilda{i})-1;
giInt{i} = spline( x, gtilda{i}, xx);

end

g_mean = {mean(giInt{1}(:)) mean(giInt{2}(:)) mean(giInt{3}(:)) mean(giInt{4}(:))};
g_var = {var(giInt{1}(:)) var(giInt{2}(:)) var(giInt{3}(:)) var(giInt{4}(:))};
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
%     figure
%     plot(xx, 10*log10(abs(g_i{i})));
%     title('Absolute value of g_i (dB)');
%     xlabel('tau');
%     ylabel('|g_i| (dB)');
end

exp_vector = zeros(1, 4);

for i = 1:samples
    exp_vector(i) = mean(abs(g_i{i}).^2); 
end

figure
stem(0:3, 10*log10(exp_vector));

%Here's the behavoiur of |h1(nTc)| for n = 0:1999 (we dropped the transient)
figure
plot(0:0.25:(0.25*1999), 10*log10(abs(g_i{2}(1, 1:2000))));
title('Plot of |g1(nTc)| (dB) for n = 0:1999');

%I'll try to plot the histogram
h1 = g_i{2}(1:inp_length)/sqrt(exp_vector(2));
bin_vector = -2:0.1:2;
n = hist(h1, bin_vector );
energy = sum(n.^2);
n = n/sqrt(energy);
figure
bar(bin_vector, n);
hold on

%Stantard complex normal pdf
z = -3:0.01:3;
pdf = exp(-abs(z).^2)/pi;
plot(z, pdf, 'm', 'LineWidth', 2);





