function cst = func_EnzymaticDetox_Filtrate2K(xi)
% Modeling toxin enzymatic degradation
% Km = 6;

global texp Texp Km

%% Parameters
dE = 0.1*xi(1); % 0.2, enzyme degradation rate, 1/hr
E0 = xi(2); % initial enzyme concentration, converts 1 ug/ml per hr (50.5 uU/ml)

%% Initial conditions
T0 = mean(Texp(:,1)); % initial toxin concentration, ug/ml

%% time step
dt = 0.02;

%% Simulating the dynamics
trng = min(texp):dt:(max(texp)+dt);
c = 1;
T(c) = T0;
E(c) = E0;
for t = trng(1:length(trng)-1)
    c = c+1;
    E(c) = E(c-1) - dt*dE*E(c-1);
    T(c) = T(c-1) - dt*E(c-1)*T(c-1)/(T(c-1)+Km); % assuming Km = 6
end

Ti = interp1(trng,T,texp);
Res = abs(Ti-Texp);
Nc = 4;
Nt = length(texp);
cst = zeros(1,Nc);
for nc = 1:4
    cst(nc) = sum(Res(nc:Nc:Nt));
end

return;