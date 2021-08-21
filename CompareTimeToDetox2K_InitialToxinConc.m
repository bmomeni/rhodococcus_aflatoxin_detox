clear

Km = 10;
Nt = 3*7200+1;
trng = linspace(0,3*120,Nt);
dts = trng(2)-trng(1);

T = zeros(1,Nt);
E = zeros(1,Nt);

%% Initial conditions
T0rng = 1:0.005:50; % initial toxin concentration, ug/ml
NT = length(T0rng);

E0rng = 1/50.5*[150:50:400];
NE = length(E0rng);

TtD = zeros(NE,NT)*NaN;

%% Simulating the dynamics
nE = 0;
for E0 = E0rng
    nE = nE+1;
    dE = max(0,1/(17.3-0.0285*50.5*E0));
    for nT = 1:NT
        c = 1;
        T(c) = T0rng(nT);
        E(c) = E0;
        trg = 0;
        for t = trng(1:length(trng)-1)
            c = c+1;
            E(c) = E(c-1) - dts*dE*E(c-1);
            T(c) = T(c-1) - dts*E(c-1)*T(c-1)/(T(c-1)+Km);
            if (T(c)<0.05*T(1))&&(trg==0)
                TtD(nE,nT) = t;
                trg = 1;
            end
        end
    end
    figure(13)
    hold on
    plot(T0rng,TtD(nE,:))
end
legend('150','200','250','300','350','400')
xlabel('Initial AFG_2 conc. (\mug/ml)')
ylabel('Time to 95% detox (hrs)')
xlim([0 30])
ylim([0 72])
