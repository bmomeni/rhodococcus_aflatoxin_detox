clear
%% The issue of multiple acceptable solutions need to be addressed
%% Different pairs of betaT and KE may produce comparable results
% fit2 and func*2 use only two degrees of freedom (dE and kT)
global texp Texp Km

load('Analysis_Background_AFG2_SizeExclusion_Pseudo_AFG2_12112019.mat','T_BG')

%% Rery glucose filtrate AFG2

%% Set 1
infile = 'SizeExclusion_AFG2_11122019';
fid = fopen(strcat(infile,'.txt'),'r');
Nr = 490; %number of time points
dt = 5/60; % measurement time-step (hours)
r = 8; % number of rows
c = 12; % number of columns
dg = 2; % data group
FL = ReadDataFromText(infile,Nr,r,c,dg);
FL_AFG2_BG = 262; % estimated from complete degradation cases with live cells

ntrng = 1:Nr; % time points to keep
KmRe = 10;
Km = KmRe;

repl = 2:4; % row # on plate
ns = 2; % column # on plate
FL_Re = shiftdim(FL(repl,ns,ntrng)-FL_AFG2_BG,2)';
T_Re = RFUtoAFG2(FL_Re)./(1/T_BG(3)*T_BG(1:length(repl),1:Nr)); % normalized to change in AFG2 without filtrate;
texp = dt*ntrng;

initial_par = [1 2];
for cnt = 1:3
    Texp = T_Re(cnt,:);
    opt_par_Re = lsqnonlin(@func_EnzymaticDetox_Filtrate2K,initial_par);
    dE_Re(cnt) = 0.1*opt_par_Re(1); % 0.2, enzyme decay rate, 1/hr
    E0_Re(cnt) = opt_par_Re(2); % initial enzyme concentration, converts 1 ug/ml per hr (50.5 U/ml)
end
disp('Rery1:')
disp(T_Re(:,1)')
disp([dE_Re; E0_Re])

DetoxEff_Re = 1 - T_Re(:,Nr)./T_Re(:,1);
disp('Detox efficiency in 72 hours (Exp):')
disp(DetoxEff_Re)

% infile = 'Filtrate_TP_AFG2_11082019';
% repl = 2:4; % row # on plate
% ns = 8; % column # on plate
ind = 1:3;
T(ind) = T_Re(1:3,1)';
D(:,ind) = [dE_Re(1:3); E0_Re(1:3)];
DE(ind) = DetoxEff_Re(1:3)';
I(ind) = 1;

%% Set 2
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

%% Rpyr starch filtrate AFG2
repl = 5:7; % row # on plate
ns = 3; % column # on plate
FL_Re = shiftdim(FL(repl,ns,ntrng)-FL_AFG2_BG,2)';
T_Re = RFUtoAFG2(FL_Re)./(1/T_BG(3)*T_BG(1:length(repl),1:Nr)); % normalized to change in AFG2 without filtrate;
texp = dt*ntrng;

initial_par = [1 2];
for cnt = 1:3
    Texp = T_Re(cnt,:);
    opt_par_Re = lsqnonlin(@func_EnzymaticDetox_Filtrate2K,initial_par);
    dE_Re(cnt) = 0.1*opt_par_Re(1); % 0.2, enzyme decay rate, 1/hr
    E0_Re(cnt) = opt_par_Re(2); % initial enzyme concentration, converts 1 ug/ml per hr (50.5 U/ml)
end
disp('Rpyr:')
disp(T_Re(:,1)')
disp([dE_Re; E0_Re])

DetoxEff_Re = 1 - T_Re(:,Nr)./T_Re(:,1);
disp('Detox efficiency in 72 hours (Exp):')
disp(DetoxEff_Re)

% infile = 'Filtrate_TP_ProteaseArrest_G2_12132019';
% repl = 3; % row # on plate
% ns = 10:11; % column # on plate
ind = 4:6;
T(ind) = T_Re(1:3,1)';
D(:,ind) = [dE_Re(1:3); E0_Re(1:3)];
DE(ind) = DetoxEff_Re(1:3)';
I(ind) = 2;

%% Set 3
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

%% Rpyr glucose filtrate AFG2
repl = 2:4; % row # on plate
ns = 6; % column # on plate
FL_Re = shiftdim(FL(repl,ns,ntrng)-FL_AFG2_BG,2)';
T_Re = RFUtoAFG2(FL_Re)./(1/T_BG(3)*T_BG(1:length(repl),1:Nr)); % normalized to change in AFG2 without filtrate;
texp = dt*ntrng;

initial_par = [1 2];
for cnt = 1:3
    Texp = T_Re(cnt,:);
    opt_par_Re = lsqnonlin(@func_EnzymaticDetox_Filtrate2K,initial_par);
    dE_Re(cnt) = 0.1*opt_par_Re(1); % 0.2, enzyme decay rate, 1/hr
    E0_Re(cnt) = opt_par_Re(2); % initial enzyme concentration, converts 1 ug/ml per hr (50.5 U/ml)
end
disp('Rpyr2:')
disp(T_Re(:,1)')
disp([dE_Re; E0_Re])

DetoxEff_Re = 1 - T_Re(:,Nr)./T_Re(:,1);
disp('Detox efficiency in 72 hours (Exp):')
disp(DetoxEff_Re)

% infile = 'Filtrate_TP_ProteaseArrest_G2_12132019';
% repl = 3; % row # on plate
% ns = 10:11; % column # on plate
ind = 7:9;
T(ind) = T_Re(1:3,1)';
D(:,ind) = [dE_Re(1:3); E0_Re(1:3)];
DE(ind) = DetoxEff_Re(1:3)';
I(ind) = 3;

% infile = 'SizeExclusion_AFG2_11122019';
% fid = fopen(strcat(infile,'.txt'),'r');
% Nr = 865; %number of time points
% dt = 5/60; % measurement time-step (hours)
% r = 8; % number of rows
% c = 12; % number of columns
% dg = 2; % data group
% FL = ReadDataFromText(infile,Nr,r,c,dg);
% FL_AFG2_BG = 262; % estimated from complete degradation cases with live cells
% 
% ntrng = 1:Nr; % time points to keep
% 
% %% Rpyr starch filtrate AFG2
% repl = 5:7; % row # on plate
% ns = 2; % column # on plate
% FL_Re = shiftdim(FL(repl,ns,ntrng)-FL_AFG2_BG,2)';
% T_Re = RFUtoAFG2(FL_Re)./(1/T_BG(3)*T_BG(1:length(repl),ntrng)); % normalized to change in AFG2 without filtrate;
% texp = dt*ntrng;
% 
% initial_par = [1.5 3];
% for cnt = 1:length(ns)
%     Texp = T_Re(cnt,:);
%     opt_par_Re = lsqnonlin(@func_EnzymaticDetox_Filtrate2,initial_par);
%     dE_Re(cnt) = 0.1*opt_par_Re(1); % 0.2, enzyme decay rate, 1/hr
%     E0_Re(cnt) = opt_par_Re(2); % initial enzyme concentration, converts 1 ug/ml per hr (50.5 U/ml)
% end
% disp('Rpyr:')
% disp(T_Re(:,1)')
% disp([dE_Re; E0_Re])
% 
% DetoxEff_Re = 1 - T_Re(:,Nr)./T_Re(:,1);
% disp('Detox efficiency in 72 hours (Exp):')
% disp(DetoxEff_Re)
% 
% % infile = 'Filtrate_Cells_TimePoints_AFG2_10252019';
% % repl = 3; % row # on plate
% % ns = 8:9; % column # on plate
% ind = 6:8;
% T(ind) = T_Re(1:3,1)';
% D(:,ind) = [dE_Re(1:3); E0_Re(1:3)];
% DE(ind) = DetoxEff_Re(1:3)';
% I(ind) = 3;

figure
plot(T,D(1,:),'bo')
xlabel('AFG2 Conc. (\mug/ml)')
ylabel('Enz. decay rate')

figure
plot(T,DE,'bo')
xlabel('AFG2 Conc. (\mug/ml)')
ylabel('Detox efficiency')

figure
plot(50.5*D(2,:),D(1,:),'bo')
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

save(strcat('Fit2K_Filtrate_AFG2_ReryAll_glucose_',infile,'ns',num2str(ns),'.mat'))
