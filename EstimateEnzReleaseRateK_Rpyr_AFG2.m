clear

load('Analysis_Background_AFG2_SizeExclusion_Pseudo_AFG2_12112019.mat','T_BG')

infile = 'Clonal_PYR_WT_G2_B1_01302020';
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
repl = 2:3; % row # on plate
ns = 7; % column # on plate
OD_R(1:2,:) = shiftdim(OD(repl,ns,ntrng),2)'-mean(mean(OD0_R(:,10:20)))+0.015;
FL_R(1:2,:) = shiftdim(FL(repl,ns,ntrng)-mean(FL(repl,ns,Nr-50:Nr),3),2)';
repl = 2; % row # on plate
ns = 8; % column # on plate
OD_R(3,:) = shiftdim(OD(repl,ns,ntrng),2)'-mean(mean(OD0_R(:,10:20)))+0.015;
FL_R(3,:) = shiftdim(FL(repl,ns,ntrng)-mean(FL(repl,ns,Nr-50:Nr),3),2)';

T_R = RFUtoAFG2(FL_R)./(1/T_BG(3)*T_BG); % normalized to change in AFG2 without filtrate;
% T_R(T_R<0) = 0;

for ct = 1:3
    OD_RS(ct,:) = smooth(OD_R(ct,:),10);
end  
%% plot results
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
%     disp([mdl.L mdl.U mdl.tmid mdl.tscale80])
%     % fit_T(nfit,:) = tox.L + (tox.U-tox.L)./(1 + exp(-4*log(3)*(t-tox.tmid)/tox.tscale80));
% end

DetoxRate = -1/dt*diff(T_R')'./OD_R(:,1:Nr-1);
% DetoxRateFit = -1/dt*diff(fit_T')'./fit_OD(1:Nr-1);
figure
plot(T_R(:,indx(1:Nr-1)),DetoxRate(:,indx(1:Nr-1)),'.')
xlabel('AFG2 Conc. (\mug/ml)')
ylabel('Detox rate (per cell)')
ylim([-2 21])
xlim([0 10])

load('Fit2K_Filtrate_AFG2_Filtrate_Cells_TimePoints_AFG2_10252019.mat','E0_Re','E0_Rp','dE_Re','dE_Rp')

KT = 10;
cs = 10;
for cnt = 1:3
    T_RS(cnt,:) = smooth(T_R(cnt,:),cs);
end
Et = -(1+KT./T_RS(:,1:Nr-1)).*diff(T_RS')'; % enz. concentration
for cnt = 1:3
    EtS(cnt,:) = smooth(Et(cnt,:),cs);
end
bE = 1./OD_RS(:,2:Nr-1).*(diff(EtS')' + diag(dE_Rp)*EtS(:,2:Nr-1)); % enz. release rate, smooth OD
bEo = 1./OD_R(:,2:Nr-1).*(diff(EtS')' + diag(dE_Rp)*EtS(:,2:Nr-1)); % enz. release rate, original
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
figure
plot(dt*ntrng,OD_RS)
xlabel('Time (hours)')
ylabel('Cell density (OD)')
xlim([0 48])
ylim([0 1.6])

figure
plot(dt*ntrng,OD_RS)
hold on
errorbar(tmrng-0.5,bESm,bESsd,'ko-')
set(gca,'YScale','log')
xlabel('Time (hours)')
ylabel('Cell density (OD)')
xlim([0 36])
ylim([0 1.5])

figure
plot(dt*ntrng,T_R)
hold on
plot(dt*ntrng,T_RS,'k')
ylim([0 30])
xlabel('Time (hours)')
ylabel('AFG_2 conc. (\mug/ml)')
ylim([0 30])
xlim([0 48])

figure
plot(dt*ntrng(1:Nr-1),ECF*EtS,'.')
% % hold on
% plot(dt*ntrng(1:Nr-1),ECF*Et,'.')
xlabel('Time (hours)')
ylabel('Enzyme conc. (uU/ml)')
xlim([0 48])
ylim([0 30])

figure
plot(dt*ntrng(2:Nr-1),ECF*bES,'.')
xlabel('Time (hours)')
ylabel('Enzyme release rate (uU/OD/hr)')
xlim([0 48])
ylim([0 40])

figure
% errorbar(tmrng-0.5,bEoSm,bEoSsd,'ko-')
% hold on
plot(dt*ntrng(2:Nr-1),ECF*bEoS,'.')
xlabel('Time (hours)')
ylabel('Enzyme release rate (uU/OD/hr)')
xlim([0 48])
ylim([0 40])

figure
plot(OD_R(2:Nr-1),ECF*bEo,'o')
xlabel('Cell density (OD)')
ylabel('Enzyme release rate (uU/OD/hr)')
set(gca,'XScale','log')
xlim([0.04 1.5])
ylim([0 35])

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
ylabel('Enzyme release rate (uU/ml/OD/hr)')
set(gca,'XScale','log')
xlim([0.04 1.2])
ylim([0 20])

Nr24 = 24/dt;
tcum = dt*ntrng(1:Nr24-1); % time
ECcum = dt*cumsum(ECF*EtS(:,1:Nr24-1),2); % in U/ml
figure
plot(tcum,ECcum)
xlabel('Time (hours)')
ylabel('Total enzyme released (uU/ml)')
xlim([0 24])
ylim([0 400])

ECcum_decay(:,1) = dt*ECF*EtS(:,1);
d0 = 0.008; % intrinsic enzyme decay rate
for ct = 2:Nr24-1
    ECcum_decay(:,ct) = ECcum_decay(:,ct-1)*exp(-d0*dt) + dt*ECF*EtS(:,ct);
end
figure
plot(tcum,ECcum_decay)
xlabel('Time (hours)')
ylabel('Total enzyme released (\muU/ml)')
xlim([0 24])
ylim([0 400])

disp([mean(ECcum(:,Nr24-1)), std(ECcum(:,Nr24-1))])
disp([mean(ECcum_decay(:,Nr24-1)), std(ECcum_decay(:,Nr24-1))])

save(strcat('EstimateEnzReleaseRate2K_K10_Rpyr_AFG2_Cells_',infile,'.mat'))

