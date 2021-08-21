clear

Km = 10;
Nt = 3*7200+1;
trng = linspace(0,120,Nt);
dts = trng(2)-trng(1);

T = zeros(1,Nt);
E = zeros(1,Nt);

%% Initial conditions
T0rng = 5:5:30; % initial toxin concentration, ug/ml
NT = length(T0rng);

E0rng = 1/50.5*(20:0.005:500);
NE = length(E0rng);

TtD = zeros(NT,NE)*NaN;

%% Simulating the dynamics
nT = 0;
for T0 = T0rng
    nT = nT+1;
    for nE = 1:NE
        c = 1;
        T(c) = T0;
        E(c) = E0rng(nE);
        dE = 1/(17.3-0.0285*50.5*E(1));
        trg = 0;
        for t = trng(1:length(trng)-1)
            c = c+1;
            E(c) = E(c-1) - dts*dE*E(c-1);
            T(c) = T(c-1) - dts*E(c-1)*T(c-1)/(T(c-1)+Km);
            if (T(c)<0.05*T(1))&&(trg==0)
                TtD(nT,nE) = t;
                trg = 1;
            end
        end
    end
    figure(12)
    hold on
    plot(50.5*E0rng,TtD(nT,:))
end
legend('5','10','15','20','25')
xlabel('Initial enz. conc. (\muU/ml)')
ylabel('Time to 95% detox (hrs)')
xlim([0 500])
ylim([0 72])
