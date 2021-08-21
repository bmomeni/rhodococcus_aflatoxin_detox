clear
%% The issue of multiple acceptable solutions need to be addressed
%% Different pairs of betaT and KE may produce comparable results
% fit2 and func*2 use only two degrees of freedom (dE and kT)
global texp Texp Km

load('Analysis_Background_AFG2_SizeExclusion_Pseudo_AFG2_12112019.mat','T_BG')

infile = 'Filtrate_Cells_TimePoints_AFG2_10252019';
fid = fopen(strcat(infile,'.txt'),'r');
Nr = 865; %number of time points
dt = 5/60; % measurement time-step (hours)
r = 8; % number of rows
c = 12; % number of columns
dg = 3; % data group
FL = ReadDataFromText(infile,Nr,r,c,dg);
FL_AFG2_BG = 262; % estimated from complete degradation cases with live cells

KmRe = 10;
KmRp = 10;

ntrng = 1:Nr; % time points to keep
% %% No cell control
% repl = 2:4; % row # on plate
% ns = 11; % column # on plate
% FL0_R = shiftdim(FL(repl,ns,ntrng),2)';
% T0_R = RFUtoAFG2(FL0_R);

%% Rpyr filtrate AFG2
repl = 2:4; % row # on plate
ns = 3; % column # on plate
FL_Rp = shiftdim(FL(repl,ns,ntrng)-FL_AFG2_BG,2)';
T_Rp = RFUtoAFG2(FL_Rp)./(1/T_BG(1)*T_BG); % normalized to change in AFG2 without filtrate
texp = dt*ntrng;
Km = KmRp;
initial_par = [1 1];
for cnt = 1:3
    Texp = T_Rp(cnt,:);
    opt_par_Rp = lsqnonlin(@func_EnzymaticDetox_Filtrate2K,initial_par);
    dE_Rp(cnt) = 0.1*opt_par_Rp(1); % 0.2, enzyme decay rate, 1/hr
    E0_Rp(cnt) = opt_par_Rp(2); % initial enzyme concentration, converts 1 ug/ml per hr (50.5 U/ml)
end
disp('Rpyr:')
disp([dE_Rp; E0_Rp])
dE_Rpm = mean(dE_Rp);
E0_Rpm = mean(E0_Rp);

%% Rery filtrate AFG2
repl = 2:4; % row # on plate
ns = 6; % column # on plate
FL_Re = shiftdim(FL(repl,ns,ntrng)-FL_AFG2_BG,2)';
T_Re = RFUtoAFG2(FL_Re)./(1/T_BG(3)*T_BG); % normalized to change in AFG2 without filtrate;
T_Re(3,:) = T_Re(3,:) - 1.8; % adjusted for systematic background shift
texp = dt*ntrng;
Km = KmRe;
initial_par = [1 2];
for cnt = 1:3
    Texp = T_Re(cnt,:);
    opt_par_Re = lsqnonlin(@func_EnzymaticDetox_Filtrate2K,initial_par);
    dE_Re(cnt) = 0.1*opt_par_Re(1); % 0.2, enzyme decay rate, 1/hr
    E0_Re(cnt) = opt_par_Re(2); % initial enzyme concentration, converts 1 ug/ml per hr (50.5 U/ml)
end
disp('Rery:')
disp([dE_Re; E0_Re])
dE_Rem = mean(dE_Re);
E0_Rem = mean(E0_Re);

%% plot results
%% Initial conditions
E0 = 0; % initial enzyme concentration, U/ml
%% time step
dts = 0.02;

%% Simulating the dynamics
T0 = mean(T_Rp(:,1)); % initial toxin concentration, ug/ml
trng = min(texp):dts:(max(texp)+dts);
c = 1;
T_Rps(c) = T0;
E(c) = E0_Rpm;
for t = trng(1:length(trng)-1)
    c = c+1;
    E(c) = E(c-1) - dts*dE_Rpm*E(c-1);
    T_Rps(c) = T_Rps(c-1) - dts*E(c-1)*T_Rps(c-1)/(T_Rps(c-1)+KmRp);
end

figure
plot(dt*ntrng,T_Rp)
hold on
plot(trng,T_Rps,'k:')
text(36,20,'Rpyr')
ylim([0 30])
xlabel('Time (hours)')
ylabel('AFG2 Conc. (\mug/ml)')

T0 = mean(T_Re(:,1)); % initial toxin concentration, ug/ml
trng = min(texp):dts:(max(texp)+dts);
c = 1;
T_Res(c) = T0;
E(c) = E0_Rem;
for t = trng(1:length(trng)-1)
    c = c+1;
    E(c) = E(c-1) - dts*dE_Rem*E(c-1);
    T_Res(c) = T_Res(c-1) - dts*E(c-1)*T_Res(c-1)/(T_Res(c-1)+KmRe);
end

figure
plot(dt*ntrng,T_Re)
hold on
plot(trng,T_Res,'k:')
text(36,20,'Rery')
ylim([0 30])
xlabel('Time (hours)')
ylabel('AFG2 Conc. (\mug/ml)')

figure
bar([1 1.5],1./[mean(dE_Rp) mean(dE_Re)],'FaceColor',[0.4 0.4 0.4])
hold on
errorbar([1 1.5],1./[mean(dE_Rp) mean(dE_Re)],1./[mean(dE_Rp) mean(dE_Re)].^2.*[std(dE_Rp) std(dE_Re)],'ko')
set(gca,'XTick', [1 1.5], 'XTickLabel',{'Rp','Re'})
xlim([0.5 2])
ylim([0 12])
ylabel('Enzyme lifetime (hrs)')
disp('T test (dE):')
[h,p] = ttest2(1./dE_Rp,1./dE_Re)
% pnp = ranksum(1./dE_Rp,1./dE_Re)

figure
bar([1 1.5],50.5*[mean(E0_Rp) mean(E0_Re)],'FaceColor',[0.4 0.4 0.4])
hold on
errorbar([1 1.5],50.5*[mean(E0_Rp) mean(E0_Re)],50.5*[std(E0_Rp) std(E0_Re)],'ko')
set(gca,'XTick', [1 1.5], 'XTickLabel',{'Rp','Re'})
xlim([0.5 2])
ylim([0 300])
ylabel('Initial enz. conc. (\muU/ml)')
disp('T test (E0):')
[h,p] = ttest2(E0_Rp,E0_Re)
% pnp = ranksum(E0_Rp,E0_Re)

% figure
% errorbar([1 1.5],[mean(KE_Rp) mean(KE_Re)],[std(KE_Rp) std(KE_Re)],'o')
% set(gca,'XTick', [1 1.5], 'XTickLabel',{'Rp','Re'})
% xlim([0.5 2])
% ylim([0 0.18])
% ylabel('Enz. sat. conc. (\mug/ml)')

DetoxEff_Rp = 1 - T_Rp(:,Nr)./T_Rp(:,1);
DetoxEff_Re = 1 - T_Re(:,Nr)./T_Re(:,1);
disp('Detox efficiency in 72 hours:')
disp([mean(DetoxEff_Rp) std(DetoxEff_Rp)])
disp([mean(DetoxEff_Re) std(DetoxEff_Re)])

save(strcat('Fit3_Filtrate_AFG2_Km_',infile,'.mat'))
