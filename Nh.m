clear all
close all
T = 1;      % Sampling time before interpolation
tau_rms = 0.3*T;       
tau = linspace(0, 5, 1000);
PDP_continuous = exp(-tau/tau_rms)./tau_rms;
n = 1;
iterations = 10000;

%I'll now sample the "continuous" PDP
Tc = (0.25*T);      %New sampling time
tau = 0:Tc:150;
samples = 6;
PDP_sampled = exp(-tau/tau_rms)./tau_rms;
%PDP_sampled = [PDP_sampled(1:samples) zeros(1, length(PDP_sampled )-samples)] ;


%normalization of the PDP
K = 3; %K in dB
K = 10^(K/10); % K in linear

Md = sum(PDP_sampled);

C = sqrt(K/(K+1));
norm = Md/( 1 - C^2 ); % MdNorm = 1 - c^2
PDPreal = PDP_sampled/norm; %PDP normalized


a = [];

for e = 0:6
    Tc = (0.25*T);      %New sampling time
    tau = 0:Tc:e/4;
    PDP_sampled = PDPreal( length(tau)+1: length(PDPreal));
    %PDP_sampled = [PDP_sampled(1:samples) zeros(1, length(PDP_sampled )-samples)] ;
    
    a(e+1) = sum(PDP_sampled);

    
end

a = (a*10).^(-1);
figure 
plot( 1 : 7, 10*log10(a));
    
    
    

