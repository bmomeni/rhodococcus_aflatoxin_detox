clear
%% The issue of multiple acceptable solutions need to be addressed
%% Different pairs of betaT and KE may produce comparable results
% fit2 and func*2 use only two degrees of freedom (dE and kT)
global texp Texp Km

load('Analysis_Background_AFG2_SizeExclusion_Pseudo_AFG2_12112019.mat','T_BG')

Nr = 577; %number of time points
dt = 5/60; % measurement time-step (hours)
r = 8; % number of rows
c = 12; % number of columns
FL_AFG2_BG = 262; % estimated from complete degradation cases with live cells

KmRp = 10; % ug/ml, estimated from DetoxByFiltrate_AFG2_glucose_InitialToxConc.m
KmRe = 10; % ug/ml, estimated from DetoxByFiltrate_AFG2_glucose_InitialToxConc.m

ntrng = 1:Nr; % time points to keep
T_BG = T_BG(1:Nr);
% %% No cell control
% repl = 2:4; % row # on plate
% ns = 11; % column # on plate
% FL0_R = shiftdim(FL(repl,ns,ntrng),2)';
% T0_R = RFUtoAFG2(FL0_R);

%% Rpyr filtrate AFG2, 0 hr
repl = 5:7; % row # on plate
ns = 7; % column # on plate
infile = 'StoredFiltrate_AFG2_09202020';
fid = fopen(strcat(infile,'.txt'),'r');
dg = 2; % data group
FL = ReadDataFromText(infile,Nr,r,c,dg);
fclose(fid);
FL_Rp0 = shiftdim(FL(repl,ns,ntrng)-FL_AFG2_BG,2)';
T_Rp0 = RFUtoAFG2(FL_Rp0)./(1/T_BG(1)*T_BG); % normalized to change in AFG2 without filtrate
texp = dt*ntrng;
Km = KmRp;
initial_par = [1 1];
for cnt = 1:3
    Texp = T_Rp0(cnt,:);
    opt_par_Rp = lsqnonlin(@func_EnzymaticDetox_Filtrate3,initial_par);
    dE_Rp0(cnt) = 0.1*opt_par_Rp(1); % 0.2, enzyme decay rate, 1/hr
    E0_Rp0(cnt) = opt_par_Rp(2); % initial enzyme concentration, converts 1 ug/ml per hr (50.5 U/ml)
end
disp('Rpyr0:')
disp([dE_Rp0; E0_Rp0])
dE_Rpm0 = mean(dE_Rp0);
dE_Rps0 = std(dE_Rp0);
E0_Rpm0 = mean(E0_Rp0);
E0_Rps0 = std(E0_Rp0);

%% Rpyr filtrate AFG2, 48 hr
repl = 5:7; % row # on plate
ns = 8; % column # on plate
infile = 'StoredFiltrate_AFG2_09202020';
fid = fopen(strcat(infile,'.txt'),'r');
dg = 4; % data group
FL = ReadDataFromText(infile,Nr,r,c,dg);
fclose(fid);
FL_Rp48 = shiftdim(FL(repl,ns,ntrng)-FL_AFG2_BG,2)';
T_Rp48 = RFUtoAFG2(FL_Rp48)./(1/T_BG(1)*T_BG); % normalized to change in AFG2 without filtrate
texp = dt*ntrng;
Km = KmRp;
initial_par = [1 1];
for cnt = 1:3
    Texp = T_Rp48(cnt,:);
    opt_par_Rp = lsqnonlin(@func_EnzymaticDetox_Filtrate3,initial_par);
    dE_Rp48(cnt) = 0.1*opt_par_Rp(1); % 0.2, enzyme decay rate, 1/hr
    E0_Rp48(cnt) = opt_par_Rp(2); % initial enzyme concentration, converts 1 ug/ml per hr (50.5 U/ml)
end
disp('Rpyr48:')
disp([dE_Rp48; E0_Rp48])
dE_Rpm48 = mean(dE_Rp48);
dE_Rps48 = std(dE_Rp48);
E0_Rpm48 = mean(E0_Rp48);
E0_Rps48 = std(E0_Rp48);

%% Rpyr filtrate AFG2, 96 hr
repl = 5:7; % row # on plate
ns = 9; % column # on plate
infile = 'StoredFiltrate_AFG2_09202020';
fid = fopen(strcat(infile,'.txt'),'r');
dg = 6; % data group
FL = ReadDataFromText(infile,Nr,r,c,dg);
fclose(fid);
FL_Rp96 = shiftdim(FL(repl,ns,ntrng)-FL_AFG2_BG,2)';
T_Rp96 = RFUtoAFG2(FL_Rp96)./(1/T_BG(1)*T_BG); % normalized to change in AFG2 without filtrate
texp = dt*ntrng;
Km = KmRp;
initial_par = [1 1];
for cnt = 1:3
    Texp = T_Rp96(cnt,:);
    opt_par_Rp = lsqnonlin(@func_EnzymaticDetox_Filtrate3,initial_par);
    dE_Rp96(cnt) = 0.1*opt_par_Rp(1); % 0.2, enzyme decay rate, 1/hr
    E0_Rp96(cnt) = opt_par_Rp(2); % initial enzyme concentration, converts 1 ug/ml per hr (50.5 U/ml)
end
disp('Rpyr96:')
disp([dE_Rp96; E0_Rp96])
dE_Rpm96 = mean(dE_Rp96);
dE_Rps96 = std(dE_Rp96);
E0_Rpm96 = mean(E0_Rp96);
E0_Rps96 = std(E0_Rp96);

%% Rery filtrate AFG2, 0 hr
repl = 2:4; % row # on plate
ns = 7; % column # on plate
infile = 'StoredFiltrate_AFG2_09202020';
fid = fopen(strcat(infile,'.txt'),'r');
dg = 2; % data group
FL = ReadDataFromText(infile,Nr,r,c,dg);
fclose(fid);
FL_Re0 = shiftdim(FL(repl,ns,ntrng)-FL_AFG2_BG,2)';
T_Re0 = RFUtoAFG2(FL_Re0)./(1/T_BG(1)*T_BG); % normalized to change in AFG2 without filtrate
texp = dt*ntrng;
Km = KmRe;
initial_par = [1 1];
for cnt = 1:3
    Texp = T_Re0(cnt,:);
    opt_par_Re = lsqnonlin(@func_EnzymaticDetox_Filtrate3,initial_par);
    dE_Re0(cnt) = 0.1*opt_par_Re(1); % 0.2, enzyme decay rate, 1/hr
    E0_Re0(cnt) = opt_par_Re(2); % initial enzyme concentration, converts 1 ug/ml per hr (50.5 U/ml)
end
disp('Rery0:')
disp([dE_Re0; E0_Re0])
dE_Rem0 = mean(dE_Re0);
dE_Res0 = std(dE_Re0);
E0_Rem0 = mean(E0_Re0);
E0_Res0 = std(E0_Re0);

%% Reyr filtrate AFG2, 48 hr
repl = 2:4; % row # on plate
ns = 8; % column # on plate
infile = 'StoredFiltrate_AFG2_09202020';
fid = fopen(strcat(infile,'.txt'),'r');
dg = 4; % data group
FL = ReadDataFromText(infile,Nr,r,c,dg);
fclose(fid);
FL_Re48 = shiftdim(FL(repl,ns,ntrng)-FL_AFG2_BG,2)';
T_Re48 = RFUtoAFG2(FL_Re48)./(1/T_BG(1)*T_BG); % normalized to change in AFG2 without filtrate
texp = dt*ntrng;
Km = KmRe;
initial_par = [1 1];
for cnt = 1:3
    Texp = T_Re48(cnt,:);
    opt_par_Re = lsqnonlin(@func_EnzymaticDetox_Filtrate3,initial_par);
    dE_Re48(cnt) = 0.1*opt_par_Re(1); % 0.2, enzyme decay rate, 1/hr
    E0_Re48(cnt) = opt_par_Re(2); % initial enzyme concentration, converts 1 ug/ml per hr (50.5 U/ml)
end
disp('Rery48:')
disp([dE_Re48; E0_Re48])
dE_Rem48 = mean(dE_Re48);
dE_Res48 = std(dE_Re48);
E0_Rem48 = mean(E0_Re48);
E0_Res48 = std(E0_Re48);

%% Reyr filtrate AFG2, 96 hr
repl = 2:4; % row # on plate
ns = 9; % column # on plate
infile = 'StoredFiltrate_AFG2_09202020';
fid = fopen(strcat(infile,'.txt'),'r');
dg = 6; % data group
FL = ReadDataFromText(infile,Nr,r,c,dg);
fclose(fid);
FL_Re96 = shiftdim(FL(repl,ns,ntrng)-FL_AFG2_BG,2)';
T_Re96 = RFUtoAFG2(FL_Re96)./(1/T_BG(1)*T_BG); % normalized to change in AFG2 without filtrate
texp = dt*ntrng;
Km = KmRe;
initial_par = [1 1];
for cnt = 1:3
    Texp = T_Re96(cnt,:);
    opt_par_Re = lsqnonlin(@func_EnzymaticDetox_Filtrate3,initial_par);
    dE_Re96(cnt) = 0.1*opt_par_Re(1); % 0.2, enzyme decay rate, 1/hr
    E0_Re96(cnt) = opt_par_Re(2); % initial enzyme concentration, converts 1 ug/ml per hr (50.5 U/ml)
end
disp('Rery96:')
disp([dE_Re96; E0_Re96])
dE_Rem96 = mean(dE_Re96);
dE_Res96 = std(dE_Re96);
E0_Rem96 = mean(E0_Re96);
E0_Res96 = std(E0_Re96);

%% plot results
figure
plot(dt*ntrng,T_Rp0,'r')
hold on
plot(dt*ntrng,T_Rp48,'color',[0.8 0.4 0.4])
plot(dt*ntrng,T_Rp96,'color',[0.5 0 0])
text(36,15,'Rpyr')
ylim([0 20])
xlabel('Time (hours)')
ylabel('AFG2 Conc. (\mug/ml)')

figure
plot(dt*ntrng,T_Re0,'b')
hold on
plot(dt*ntrng,T_Re48,'color',[0.4 0.4 0.8])
plot(dt*ntrng,T_Re96,'color',[0 0 0.5])
text(36,15,'Rery')
ylim([0 20])
xlabel('Time (hours)')
ylabel('AFG2 Conc. (\mug/ml)')

figure
bar(1,1./dE_Rpm0,0.3,'FaceColor',[1 0 0])
hold on
bar(1.5,1./dE_Rpm48,0.3,'FaceColor',[0.8 0.4 0.4])
bar(2,1./dE_Rpm96,0.3,'FaceColor',[0.5 0 0])
errorbar([1 1.5 2],1./[dE_Rpm0 dE_Rpm48 dE_Rpm96],1./[dE_Rpm0 dE_Rpm48 dE_Rpm96].^2.*[dE_Rps0 dE_Rps48 dE_Rps96],'ko')
set(gca,'XTick', [1 1.5 2], 'XTickLabel',{'Rp0','Rp48','Rp96'})
xlim([0.5 2.5])
ylim([0 15])
ylabel('Enz. lifetime (hr)')

figure
bar(1,50.5*E0_Rpm0,0.3,'FaceColor',[1 0 0])
hold on
bar(1.5,50.5*E0_Rpm48,0.3,'FaceColor',[0.8 0.4 0.4])
bar(2,50.5*E0_Rpm96,0.3,'FaceColor',[0.5 0 0])
errorbar([1 1.5 2],50.5*[E0_Rpm0 E0_Rpm48 E0_Rpm96],50.5*[E0_Rps0 E0_Rps48 E0_Rps96],'ko')
set(gca,'XTick', [1 1.5 2], 'XTickLabel',{'Rp0','Rp48','Rp96'})
xlim([0.5 2.5])
ylim([0 180])
ylabel('Initial enz. conc. (U/ml)')

figure
bar(1,1./dE_Rem0,0.3,'FaceColor',[0 0 1])
hold on
bar(1.5,1./dE_Rem48,0.3,'FaceColor',[0.4 0.4 0.8])
bar(2,1./dE_Rem96,0.3,'FaceColor',[0 0 0.5])
errorbar([1 1.5 2],1./[dE_Rem0 dE_Rem48 dE_Rem96],1./[dE_Rem0 dE_Rem48 dE_Rem96].^2.*[dE_Res0 dE_Res48 dE_Res96],'ko')
set(gca,'XTick', [1 1.5 2], 'XTickLabel',{'Re0','Re48','Re96'})
xlim([0.5 2.5])
ylim([0 15])
ylabel('Enz. lifetime (hr)')

figure
bar(1,50.5*E0_Rem0,0.3,'FaceColor',[0 0 1])
hold on
bar(1.5,50.5*E0_Rem48,0.3,'FaceColor',[0.4 0.4 0.8])
bar(2,50.5*E0_Rem96,0.3,'FaceColor',[0 0 0.5])
errorbar([1 1.5 2],50.5*[E0_Rem0 E0_Rem48 E0_Rem96],50.5*[E0_Res0 E0_Res48 E0_Res96],'ko')
set(gca,'XTick', [1 1.5 2], 'XTickLabel',{'Re0','Re48','Re96'})
xlim([0.5 2.5])
ylim([0 180])
ylabel('Initial enz. conc. (U/ml)')

save(strcat('StoredFiltrate_AFG2_Km10_',infile,'.mat'))
