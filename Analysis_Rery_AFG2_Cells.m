clear

load('Analysis_Background_AFG2_SizeExclusion_Pseudo_AFG2_12112019.mat','T_BG')

infile = 'Clonal_ERY_WT_G2_B1_03032020';
Nr = 865; %number of time points
dt = 5/60; % measurement time-step (hours)
r = 8; % number of rows
c = 12; % number of columns
dg = 1; % data group
OD = ReadDataFromText(infile,Nr,r,c,dg);
dg = 2; % data group
FL = ReadDataFromText(infile,Nr,r,c,dg);

ntrng = 1:Nr; % time points to keep
%% No cell control
repl = 2:4; % row # on plate
ns = 11; % column # on plate
OD0_R = shiftdim(OD(repl,ns,ntrng),2)';
FL0_R = shiftdim(FL(repl,ns,ntrng),2)';
T0_R = RFUtoAFG2(FL0_R);

%% Rpyr filtrate AFG2
repl = 2:4; % row # on plate
ns = 8; % column # on plate
OD_R = shiftdim(OD(repl,ns,ntrng),2)'-mean(mean(OD0_R(:,10:20)))+0.01;
FL_R = shiftdim(FL(repl,ns,ntrng)-mean(FL(repl,ns,Nr-50:Nr),3),2)';
T_R = RFUtoAFG2(FL_R)./(1/T_BG(3)*T_BG); % normalized to change in AFG2 without filtrate;

%% plot results
figure
plot(dt*ntrng,OD_R)
xlabel('Time (hours)')
ylabel('OD')
xlim([0 48])

figure
plot(dt*ntrng,T_R)
% hold on
% plot(dt*ntrng,T0_R,'k')
ylim([0 30])
xlabel('Time (hours)')
ylabel('AFG_2 conc. (\mug/ml)')
ylim([0 30])
xlim([0 48])

indx = (mean(T_R)<mean(T0_R));
figure
plot(OD_R(:,indx)',T_R(:,indx)','.')
xlabel('OD')
ylabel('AFG2 Conc. (\mug/ml)')
ylim([0 30])

rng = [450 450 450];
% rngT = [250 200 300];
for nfit = 1:3
    nODrng = 1:rng(nfit);
    % nTrng = rngT(nfit):Nr;
    ft = fittype('L + (U-L)/(1 + exp(-4*log(3)*(x-tmid)/tscale80))','indep','x');
    mdl = fit(dt*ntrng(nODrng)',OD_R(nfit,nODrng)',ft,'start',[0.01, 1, 30,3]);
    % tox = fit(dt*ntrng(nTrng)',T_R(nfit,nTrng)',ft,'start',[30, -30, 30,3]);
    t = dt*ntrng;
    fit_OD(nfit,:) = mdl.L + (mdl.U-mdl.L)./(1 + exp(-4*log(3)*(t-mdl.tmid)/mdl.tscale80));
    % fit_T(nfit,:) = tox.L + (tox.U-tox.L)./(1 + exp(-4*log(3)*(t-tox.tmid)/tox.tscale80));
end

DetoxRate = -1/dt*diff(T_R')'./fit_OD(:,1:Nr-1);
% DetoxRateFit = -1/dt*diff(fit_T')'./fit_OD(1:Nr-1);
figure
plot(T_R(:,indx(1:Nr-1)),DetoxRate(:,indx(1:Nr-1)),'.')
xlabel('AFG2 Conc. (\mug/ml)')
ylabel('Detox rate (per cell)')
ylim([-2 21])
xlim([0 10])

figure
plot(dt*ntrng,fit_OD)
xlabel('Time (hours)')
ylabel('OD')
xlim([0 48])

trng = dt*ntrng;
S = mean(fit_OD);
T0 = mean(T_R(:,1));
E0 = 0;
dE = 0; %0.1342;
Km = 6;
bE = 0.3/50.5;
[Ts1, Es1, ~] = DetoxKineticsCells(S,T0,E0,dE,Km,bE,trng);
bE = 1/50.5;
[Ts2, Es2, ~] = DetoxKineticsCells(S,T0,E0,dE,Km,bE,trng);
bE = 3/50.5;
[Ts3, Es3, ~] = DetoxKineticsCells(S,T0,E0,dE,Km,bE,trng);
bE = 10/50.5;
[Ts4, Es4, tsrng] = DetoxKineticsCells(S,T0,E0,dE,Km,bE,trng);
figure
plot(trng,T_R)
hold on
plot(tsrng,Ts1,':','color',[0.3 0.3 0.3])
plot(tsrng,Ts2,':','color',[0.4 0.4 0.4])
plot(tsrng,Ts3,':','color',[0.5 0.5 0.5])
plot(tsrng,Ts4,':','color',[0.6 0.6 0.6])
ylim([0 30])
xlabel('Time (hours)')
ylabel('AFG_2 conc. (\mug/ml)')
ylim([0 30])
xlim([0 48])

save(strcat('Analysis_Rery_AFG2_Cells_',infile,'.mat'))

