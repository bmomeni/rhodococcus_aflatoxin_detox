clear
T0rng = 5:5:50;
NT0 = length(T0rng);
dts = 0.05;
trng = 0:dts:72;
Nt = length(trng);
E = zeros(1,Nt);
T_Rps = zeros(NT0,Nt);
T_Res = zeros(NT0,Nt);
DegEff_Rp = zeros(1,NT0);
DegEff_Re = zeros(1,NT0);
Km = 0.1;

E0_Rpm = 1.7203;
dE_Rpm = 0.1010;

E0_Rem = 3.0161;
dE_Rem = 0.1342;

cnt = 0;
for T0 = T0rng
    cnt = cnt+1;
    
    
    %% Simulating the dynamics
    % T0 = mean(T_Rp(:,1)); % initial toxin concentration, ug/ml
    c = 1;
    T_Rps(cnt,c) = T0;
    E(c) = E0_Rpm;
    for t = trng(1:length(trng)-1)
        c = c+1;
        E(c) = E(c-1) - dts*dE_Rpm*E(c-1);
        T_Rps(cnt,c) = T_Rps(cnt,c-1) - dts*E(c-1)*T_Rps(cnt,c-1)/(T_Rps(cnt,c-1)+Km);
    end
    DegEff_Rp(cnt) = (1-T_Rps(cnt,Nt)/T0rng(cnt))*100;
    
    % T0 = mean(T_Re(:,1)); % initial toxin concentration, ug/ml
    c = 1;
    T_Res(cnt,c) = T0;
    E(c) = E0_Rem;
    for t = trng(1:length(trng)-1)
        c = c+1;
        E(c) = E(c-1) - dts*dE_Rem*E(c-1);
        T_Res(cnt,c) = T_Res(cnt,c-1) - dts*E(c-1)*T_Res(cnt,c-1)/(T_Res(cnt,c-1)+Km);
    end
    DegEff_Re(cnt) = (1-T_Res(cnt,Nt)/T0rng(cnt))*100;
end

figure
plot(trng,T_Rps,'k:')
text(36,45,'Rpyr')
xlim([0 85])
ylim([0 50])
xlabel('Time (hours)')
ylabel('AFG2 Conc. (\mug/ml)')
for cc = NT0-5:NT0
    text(73,T_Rps(cc,Nt),strcat(num2str(round(DegEff_Rp(cc))),'%'))
end

figure
plot(trng,T_Res,'k:')
text(36,45,'Rery')
xlim([0 85])
ylim([0 50])
xlabel('Time (hours)')
ylabel('AFG2 Conc. (\mug/ml)')
for cc = NT0-5:NT0
    text(73,T_Res(cc,Nt),strcat(num2str(round(DegEff_Re(cc))),'%'))
end

figure
plot(T0rng,DegEff_Rp)
xlim([0 55])
ylim([0 105])
set(gca,'XTick',[0 10 20 30 40 50])
xlabel('Initial AFG2 Conc. (\mug/ml)')
ylabel('Deg. efficiency at 72 hours')

figure
plot(T0rng,DegEff_Re)
xlim([0 55])
ylim([0 105])
set(gca,'XTick',[0 10 20 30 40 50])
xlabel('Initial AFG2 Conc. (\mug/ml)')
ylabel('Deg. efficiency at 72 hours')
