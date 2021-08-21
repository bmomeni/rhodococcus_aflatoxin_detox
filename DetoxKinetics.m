function [Ts, Es] = DetoxKinetics(E0,dE,Km,T0,trng)

Nt = length(trng);
dts = trng(2)-trng(1);

Ts = zeros(1,Nt);
Es = zeros(1,Nt);

%% Simulating the dynamics
nT = 0;
c = 1;
Ts(c) = T0;
Es(c) = E0;
for t = trng(1:length(trng)-1)
    c = c+1;
    Es(c) = Es(c-1) - dts*dE*Es(c-1);
    Ts(c) = Ts(c-1) - dts*Es(c-1)*Ts(c-1)/(Ts(c-1)+Km);
end
return;