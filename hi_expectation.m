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
stem(tau, PDP_sampled/norm, 'm');
title('Sampled Normalized PDP (T_{sample} = Tc)');

%filter to create g

%Classical Doppler
%Here we build a proper white noise

fd = 5*10.^(-3);

iter = 500;
trans = 200; % length of the considered transient (see hds)
inp_length = 1000;
out_length = (inp_length - trans)*samples - 3;
gi_matrix = zeros(iter, out_length);
Gi = {gi_matrix gi_matrix gi_matrix gi_matrix};
xx = [];

for j = 1:iter

    gtilda = cell(samples, 1);
    giInt = cell(samples, 1);
    
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
%         figure
%         plot(xx, 10*log10(abs(g_i{i})));
%         title('Absolute value of g_i (dB)');
%         xlabel('tau');
%         ylabel('|g_i| (dB)');
    end
    
%     length(g_i{1})
    for i = 1:samples
    Gi{i}(j, :) = abs(g_i{i}).^2;
    end

end

exp_vector = zeros(1, out_length);
for j = 1:samples
    for i = 1:out_length
        exp_vector(i) = mean(Gi{j}(:, i));
    end
    figure
    plot(xx, 10*log10(exp_vector));
end

