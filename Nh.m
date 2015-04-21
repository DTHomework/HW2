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
PDP_sampled = PDP_sampled/norm; %PDP normalized

for i = 1 : iterations
    realh{i} = PDP_sampled.'/2.*randn(length(PDP_sampled),1) + 1i*PDP_sampled.'/2.*randn(length(PDP_sampled),1) ;
end

a = [];

for e = 0:20
    Tc = (0.25*T);      %New sampling time
    tau = 0:Tc:e/4;
    PDP_sampled = exp(-tau/tau_rms)./tau_rms;
    %PDP_sampled = [PDP_sampled(1:samples) zeros(1, length(PDP_sampled )-samples)] ;
    
    Md = sum(PDP_sampled);

    C = sqrt(K/(K+1));
    norm = Md/( 1 - C^2 ); % MdNorm = 1 - c^2
    PDP_sampled = PDP_sampled/norm; %PDP normalized

    for i = 1 : iterations
        enstimate{i} = [PDP_sampled.'/2.*randn(length(PDP_sampled),1)  + 1i*PDP_sampled.'/2.*randn(length(PDP_sampled),1);  zeros(length(realh{i})-length(PDP_sampled),1)];
    end
    
    for i = 1 : iterations
        delta{i} = abs( realh{i} - enstimate{i} ).^2;
    end
    
    error = 0;
    for i = 1 : iterations
        error = error + sum( delta{i} );
    end
    error = error/iterations;
    a(e+1) = error;
end

a = (a*10).^(-1);
figure 
plot( 1 : 21, a);
    
    
    

