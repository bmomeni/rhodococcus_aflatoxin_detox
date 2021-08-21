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

% AFG2 mol weight = 330 g/mol
% 1 U catalyzes 1 umol per min or 19.8 mg/hr
% Amount of enzyme that catalyzes 1 ug/ml/hr is 50.5 uU/ml
ECF = 50.5; % Conversion factor, converts 1 ug/ml per hr (50.5 uU/ml)

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

% rng = [450 450 450];
% % rngT = [250 200 300];
% for nfit = 1:3
%     nODrng = 1:rng(nfit);
%     % nTrng = rngT(nfit):Nr;
%     ft = fittype('L + (U-L)/(1 + exp(-4*log(3)*(x-tmid)/tscale80))','indep','x');
%     mdl = fit(dt*ntrng(nODrng)',OD_R(nfit,nODrng)',ft,'start',[0.01, 1, 30,3]);
%     % tox = fit(dt*ntrng(nTrng)',T_R(nfit,nTrng)',ft,'start',[30, -30, 30,3]);
%     t = dt*ntrng;
%     fitOD(nfit,:) = mdl.L + (mdl.U-mdl.L)./(1 + exp(-4*log(3)*(t-mdl.tmid)/mdl.tscale80));
%     % fit_T(nfit,:) = tox.L + (tox.U-tox.L)./(1 + exp(-4*log(3)*(t-tox.tmid)/tox.tscale80));
% end

DetoxRate = -1/dt*diff(T_R')'./OD_R(:,1:Nr-1); % ug/ml/hr/OD
% DetoxRateFit = -1/dt*diff(fit_T')'./fit_OD(1:Nr-1);
% figure
% plot(T_R(:,indx(1:Nr-1)),DetoxRate(:,indx(1:Nr-1)),'.')
% xlabel('AFG2 Conc. (\mug/ml)')
% ylabel('Detox rate (per cell)')
% ylim([-2 21])
% xlim([0 10])

load('Fit2K_Filtrate_AFG2_Filtrate_Cells_TimePoints_AFG2_10252019.mat','E0_Re','E0_Rp','dE_Re','dE_Rp')

KT = 10;
cs = 10;
for ct = 1:3
    OD_RS(ct,:) = smooth(OD_R(ct,:),10);
end  
for cnt = 1:3
    T_RS(cnt,:) = smooth(T_R(cnt,:),cs);
end
Et = -(1+KT./T_RS(:,1:Nr-1)).*diff(T_RS')';
for cnt = 1:3
    EtS(cnt,:) = smooth(Et(cnt,:),cs);
end
bE = 1./OD_RS(:,2:Nr-1).*(diff(EtS')' + diag(dE_Re)*EtS(:,2:Nr-1));
bEo = 1./OD_R(:,2:Nr-1).*(diff(EtS')' + diag(dE_Re)*EtS(:,2:Nr-1));
for cnt = 1:3
    bES(cnt,:) = smooth(bE(cnt,:),cs);
    bEoS(cnt,:) = smooth(bEo(cnt,:),cs);
end
tmrng = 1:36;
for tm = tmrng
    tmi = (dt*ntrng<tm)&(dt*ntrng>=tm-1);
    bESm(tm) = mean(mean(bES(:,tmi(2:Nr-1))));
    bESsd(tm) = std(reshape(bES(:,tmi(2:Nr-1)),1,3*sum(tmi(2:Nr-1))));
    bEoSm(tm) = mean(mean(bEoS(:,tmi(2:Nr-1))));
    bEoSsd(tm) = std(reshape(bEoS(:,tmi(2:Nr-1)),1,3*sum(tmi(2:Nr-1))));
end

%% plot results
figure
plot(dt*ntrng,OD_RS)
set(gca,'YScale','log')
xlabel('Time (hours)')
ylabel('OD')
xlim([0 36])
ylim([0.01 1])

nts = 288; % sampling time point
figure
plot(dt*ntrng(1:Nr-1),ECF*EtS,'.')
% hold on
% errorbar(tmrng(ts)-0.5,mean(ECF*EtS(:,nts)),std(ECF*EtS(:,nts)),'ko-')
xlabel('Time (hours)')
ylabel('Enzyme conc. (U/ml)')
xlim([0 48])
ylim([0 40])

% figure
% plot(dt*ntrng(2:Nr-1),ECF*bES,'.')
% xlabel('Time (hours)')
% ylabel('Enzyme release rate (U/OD/hr)')
% xlim([0 36])
% ylim([0 30])

figure
% errorbar(tmrng-0.5,bEoSm,bEoSsd,'ko-')
% hold on
plot(dt*ntrng(2:Nr-1),ECF*bEoS,'.')
xlabel('Time (hours)')
ylabel('Enzyme release rate (U/OD/hr)')
xlim([0 48])
ylim([0 30])

Nod = 20;
bEoSsdOD = zeros(1,Nod-1);
odrng = logspace(log10(0.05),0,Nod);
for ctod = 1:Nod-1
    binbEoSOD = [];
    for ct = 1:3
        tmi = (OD_R(ct,:)<odrng(ctod+1))&(OD_R(ct,:)>=odrng(ctod));
        bEoSmOD(ct,ctod) = mean(mean(bEoS(ct,tmi(2:Nr-1))));
        binbEoSOD = [binbEoSOD, bEoS(ct,tmi(2:Nr-1))];
    end
    bEoSsdOD(ctod) = std(binbEoSOD);
end

figure
errorbar(odrng(1:Nod-1),ECF*mean(bEoSmOD),ECF*bEoSsdOD,'o-')
xlabel('Cell density (OD)')
ylabel('Enzyme release rate (\muU/OD/hr)')
set(gca,'XScale','log')
xlim([0.04 1.2])
ylim([0 30])

Nr36 = 36/dt;
tcum = dt*ntrng(1:Nr36-1); % time
ECcum = dt*cumsum(ECF*EtS(:,1:Nr36-1),2); % in uU/ml
figure
plot(tcum,ECcum)
xlabel('Time (hours)')
ylabel('Total enzyme released (\muU/ml)')
xlim([0 40])
ylim([0 400])

ECcum_decay(:,1) = dt*ECF*EtS(:,1);
d0 = 0.008; % intrinsic enzyme decay rate
for ct = 2:Nr36-1
    ECcum_decay(:,ct) = ECcum_decay(:,ct-1)*exp(-d0*dt) + dt*ECF*EtS(:,ct);
end
figure
plot(tcum,ECcum_decay)
xlabel('Time (hours)')
ylabel('Total enzyme released (\muU/ml)')
xlim([0 36])
ylim([0 400])

Nr36 = 36/dt;
disp([mean(ECcum(:,Nr36-1)), std(ECcum(:,Nr36-1))])
disp([mean(ECcum_decay(:,Nr36-1)), std(ECcum_decay(:,Nr36-1))])

save(strcat('EstimateEnzReleaseRateK_Rery_AFG2_Cells_',infile,'.mat'))

