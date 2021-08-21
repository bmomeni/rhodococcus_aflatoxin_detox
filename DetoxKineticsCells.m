function [Ts, Es, tsrng] = DetoxKineticsCells(S,T0,E0,dE,Km,bE,trng)
% detoxification by cells, assuming a fixed enzyme production rate per cell

%% Detox parameters
% bE = 8; % enzyme release rate
% Km = 6; % enzyme's Michaelis Menten coefficient
%% R-pyr
% dE = 0.101; % enzyme degradation rate, 1/hr
% E0 = 0; %1.72; % initial enzyme concentration, converts 1 ug/ml per hr (50.5 uU/ml)
%% R-ery
% dE = 0.134; % enzyme degradation rate, 1/hr
% E0 = 0; %3.02; % initial enzyme concentration, converts 1 ug/ml per hr (50.5 uU/ml)

%% initialize
Nts = 10*length(trng);
tsrng = linspace(min(trng),max(trng),Nts);
dts = trng(2)-trng(1);

Ts = zeros(1,Nts);
Es = zeros(1,Nts);

%% Simulating the dynamics
c = 1;
Ts(c) = T0;
Es(c) = E0;
Ss(c) = S(1);
for t = tsrng(1:Nts-1)
    c = c+1;
    Ss(c) = interp1(trng,S,t);
    Es(c) = Es(c-1) - dts*dE*Es(c-1) + dts*bE*Ss(c);
    Ts(c) = Ts(c-1) - dts*Es(c-1)*Ts(c-1)/(Ts(c-1)+Km);
    if Ts(c)<1e-6
        Ts(c) = 0;
    end
end

return;