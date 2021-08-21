clear
%% The issue of multiple acceptable solutions need to be addressed
%% Different pairs of betaT and KE may produce comparable results
% fit2 and func*2 use only two degrees of freedom (dE and kT)
global texp Texp Km

load('Analysis_Background_AFG2_SizeExclusion_Pseudo_AFG2_12112019.mat','T_BG')

%% Set 1
infile = 'Filtrate_TP_AFG2_11082019';
fid = fopen(strcat(infile,'.txt'),'r');
Nr = 865; %number of time points
dt = 5/60; % measurement time-step (hours)
r = 8; % number of rows
c = 12; % number of columns
dg = 2; % data group
FL = ReadDataFromText(infile,Nr,r,c,dg);
FL_AFG2_BG = 262; % estimated from complete degradation cases with live cells

ntrng = 1:Nr; % time points to keep
KmRp = 10;

%% Rpyr starch filtrate AFG2
repl = 2:4; % row # on plate
ns = 8; % column # on plate
FL_Rp = shiftdim(FL(repl,ns,ntrng)-FL_AFG2_BG,2)';
T_Rp = RFUtoAFG2(FL_Rp)./(1/T_BG(3)*T_BG(1:length(repl),1:Nr)); % normalized to change in AFG2 without filtrate;
texp = dt*ntrng;
Km = KmRp;

initial_par = [1 4];
for cnt = 1:3
    Texp = T_Rp(cnt,:);
    opt_par_Rp = lsqnonlin(@func_EnzymaticDetox_Filtrate2K,initial_par);
    dE_Rp(cnt) = 0.1*opt_par_Rp(1); % 0.2, enzyme decay rate, 1/hr
    E0_Rp(cnt) = opt_par_Rp(2); % initial enzyme concentration, converts 1 ug/ml per hr (50.5 U/ml)
end
disp('Rpyr1:')
disp(T_Rp(:,1)')
disp([dE_Rp; E0_Rp])

DetoxEff_Rp = 1 - T_Rp(:,Nr)./T_Rp(:,1);
disp('Detox efficiency in 72 hours (Exp):')
disp(DetoxEff_Rp)

ind = 1:3;
T(ind) = T_Rp(1:3,1)';
D(:,ind) = [dE_Rp(1:3); E0_Rp(1:3)];
DE(ind) = DetoxEff_Rp(1:3)';
I(ind) = 1;

%% Set 2
infile = 'Filtrate_TP_ProteaseArrest_G2_12132019';
fid = fopen(strcat(infile,'.txt'),'r');
Nr = 865; %number of time points
dt = 5/60; % measurement time-step (hours)
r = 8; % number of rows
c = 12; % number of columns
dg = 2; % data group
FL = ReadDataFromText(infile,Nr,r,c,dg);
FL_AFG2_BG = 262; % estimated from complete degradation cases with live cells

ntrng = 1:Nr; % time points to keep

%% Rpyr starch filtrate AFG2
repl = 3; % row # on plate
ns = 10:11; % column # on plate
FL_Rp = shiftdim(FL(repl,ns,ntrng)-FL_AFG2_BG,1);
T_Rp = RFUtoAFG2(FL_Rp)./(1/T_BG(3)*T_BG(1:length(ns),1:Nr)); % normalized to change in AFG2 without filtrate;
texp = dt*ntrng;

initial_par = [1 4];
for cnt = 1:length(ns)
    Texp = T_Rp(cnt,:);
    opt_par_Rp = lsqnonlin(@func_EnzymaticDetox_Filtrate2K,initial_par);
    dE_Rp(cnt) = 0.1*opt_par_Rp(1); % 0.2, enzyme decay rate, 1/hr
    E0_Rp(cnt) = opt_par_Rp(2); % initial enzyme concentration, converts 1 ug/ml per hr (50.5 U/ml)
end
disp('Rpyr:')
disp(T_Rp(:,1)')
disp([dE_Rp; E0_Rp])

DetoxEff_Rp = 1 - T_Rp(:,Nr)./T_Rp(:,1);
disp('Detox efficiency in 72 hours (Exp):')
disp(DetoxEff_Rp)

ind = 4:5;
T(ind) = T_Rp(1:2,1)';
D(:,ind) = [dE_Rp(1:2); E0_Rp(1:2)];
DE(ind) = DetoxEff_Rp(1:2)';
I(ind) = 2;

infile = 'Filtrate_Cells_TimePoints_AFG2_10252019';
fid = fopen(strcat(infile,'.txt'),'r');
Nr = 865; %number of time points
dt = 5/60; % measurement time-step (hours)
r = 8; % number of rows
c = 12; % number of columns
dg = 3; % data group
FL = ReadDataFromText(infile,Nr,r,c,dg);
FL_AFG2_BG = 262; % estimated from complete degradation cases with live cells

ntrng = 1:Nr; % time points to keep

%% Rpyr starch filtrate AFG2
repl = 3; % row # on plate
ns = 8:9; % column # on plate
FL_Rp = shiftdim(FL(repl,ns,ntrng)-FL_AFG2_BG,1);
T_Rp = RFUtoAFG2(FL_Rp)./(1/T_BG(3)*T_BG(1:length(ns),ntrng)); % normalized to change in AFG2 without filtrate;
texp = dt*ntrng;

initial_par = [1 4];
for cnt = 1:length(ns)
    Texp = T_Rp(cnt,:);
    opt_par_Rp = lsqnonlin(@func_EnzymaticDetox_Filtrate2K,initial_par);
    dE_Rp(cnt) = 0.1*opt_par_Rp(1); % 0.2, enzyme decay rate, 1/hr
    E0_Rp(cnt) = opt_par_Rp(2); % initial enzyme concentration, converts 1 ug/ml per hr (50.5 U/ml)
end
disp('Rpyr:')
disp(T_Rp(:,1)')
disp([dE_Rp; E0_Rp])

DetoxEff_Rp = 1 - T_Rp(:,Nr)./T_Rp(:,1);
disp('Detox efficiency in 72 hours (Exp):')
disp(DetoxEff_Rp)

% infile = 'Filtrate_Cells_TimePoints_AFG2_10252019';
% repl = 3; % row # on plate
% ns = 8:9; % column # on plate
ind = 6:7;
T(ind) = [T_Rp(1:2,1)'];
D(:,ind) = [dE_Rp(1:2); E0_Rp(1:2)];
DE(ind) = DetoxEff_Rp(1:2)';
I(ind) = 3;

figure
plot(T,D(1,:),'sr')
xlabel('AFG2 Conc. (\mug/ml)')
ylabel('Enz. decay rate')

figure
plot(T,DE,'sr')
xlabel('AFG2 Conc. (\mug/ml)')
ylabel('Detox efficiency')

figure
plot(50.5*D(2,:),D(1,:),'sr')
xlabel('Initial enz. conc. (\muU/ml)')
ylabel('Enz. decay rate')

%% plot results
% figure
% plot(texp,T_Rp)
% hold on
% plot(trng,T_Rps,':')
% text(36,20,'Rpyr')
% ylim([0 30])
% xlabel('Time (hours)')
% ylabel('AFG2 Conc. (\mug/ml)')

% %% verify fit
% E0_Rpst = 5;
% dE_Rpst = 0.22;
% initial_par = [1 4];
% T0 = T_Rp(1,1); % initial toxin concentration, ug/ml
% c = 1;
% T_Rpst(c) = T0;
% E(c) = E0_Rpst;
% for t = trng(1:length(trng)-1)
%     c = c+1;
%     E(c) = E(c-1) - dts*dE_Rpst*E(c-1);
%     T_Rpst(c) = T_Rpst(c-1) - dts*E(c-1)*T_Rpst(c-1)/(T_Rpst(c-1)+Km);
% end
% 
% figure
% plot(trng,T_Rpst,'r')
% text(36,20,'Rpyr')
% text(74,T_Rpst(c),num2str(round(100*(1-T_Rpst(c)/T_Rpst(1)))))
% ylim([0 30])
% xlabel('Time (hours)')
% ylabel('AFG2 Conc. (\mug/ml)')

save(strcat('Fit2K_Filtrate_AFG2_RpyrAll_starch.mat'))
