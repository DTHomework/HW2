clear all
close all
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

%normalization of the PDP
K = 3; %K in dB
K = 10^(K/10); % K in linear

Md = sum(PDP_sampled);

C = sqrt(K/(K+1));

norm = Md/( 1 - C^2 ); % MdNorm = 1 - c^2

MdNorm = Md/norm; %PDP normalized

sum( MdNorm ) + C^2;
PDP_sampled = PDP_sampled/norm;

%2nd method

% T_conf = 1*Tc;
% inp_length = 100000;
% %White noise
% w_i = wgn(1, inp_length, 0, 'complex');
% %NarrowBand Filter
% gtilda = Hds1(w_i);
% 
% %Cubic Interpolator
% x = 0:length(gtilda)-1;
% xx = 0: (1/80)*Tc : length(gtilda)-1;
% giInt = interp1( x, gtilda, xx, 'spline');
% 
% sigma_i = sqrt(PDP_sampled(2));
% 
% %g_i
% 
% g_i = sigma_i * giInt;
% sample_step = 1000;
% vec = 1:sample_step:1000*sample_step;
% h151 = g_i(vec);
% length(h151)

%1st method

iter = 1000
inp_length = 1000;
h151 = zeros(iter, 1);

for i = 1:iter
%White noise
w_i = wgn(1, inp_length, 0, 'complex');
%NarrowBand Filter
gtilda = Hds1(w_i);

%Cubic Interpolator
x = 0:length(gtilda)-1;
xx = 0: (1/80)*Tc : length(gtilda)-1;
giInt = interp1( x, gtilda, xx, 'spline');

sigma_i = sqrt(PDP_sampled(2));

%g_i

g_i = sigma_i * giInt;

h151(i) = g_i(152); 
end


h1 = abs(h151);
radEh1 = sqrt(mean(abs(h151).^2));
f1 = h1./radEh1;
bin_vector = 0:0.1:3;
n = hist(f1, bin_vector );
area = sum(n*0.1);
n = n/area;
figure
bar(bin_vector, n);
hold on

%Stantard complex normal pdf
z = 0:0.01:3;
pdf = 2*z.*exp(-z.^2);
energy = sum(pdf*0.01);
plot(z, pdf, 'm', 'LineWidth', 2);