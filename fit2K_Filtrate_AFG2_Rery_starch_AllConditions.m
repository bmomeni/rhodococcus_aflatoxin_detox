clear
%% The issue of multiple acceptable solutions need to be addressed
%% Different pairs of betaT and KE may produce comparable results
% fit2 and func*2 use only two degrees of freedom (dE and kT)
global texp Texp Km

load('Analysis_Background_AFG2_SizeExclusion_Pseudo_AFG2_12112019.mat','T_BG')

%% set 1
infile = 'StarchCultures_Supernatant_AFG2_06212019';
fid = fopen(strcat(infile,'.txt'),'r');
Nr = 865; %number of time points
dt = 5/60; % measurement time-step (hours)
r = 8; % number of rows
c = 12; % number of columns
dg = 2; % data group
FL = ReadDataFromText(infile,Nr,r,c,dg);
FL_AFG2_BG = 262; % estimated from complete degradation cases with live cells

ntrng = 1:Nr; % time points to keep
KmRe = 10;
Km = KmRe;

%% Rpyr starch filtrate AFG2
repl = 5:7; % row # on plate
ns = 6; % column # on plate
FL_Re = shiftdim(FL(repl,ns,ntrng)-FL_AFG2_BG,2)';
T_Re = RFUtoAFG2(FL_Re)./(1/T_BG(3)*T_BG(1:length(repl),1:Nr)); % normalized to change in AFG2 without filtrate;
texp = dt*ntrng;

initial_par = [1.5 3];
for cnt = 1:3
    Texp = T_Re(cnt,:);
    opt_par_Re = lsqnonlin(@func_EnzymaticDetox_Filtrate2K,initial_par);
    dE_Re(cnt) = 0.1*opt_par_Re(1); % 0.2, enzyme decay rate, 1/hr
    E0_Re(cnt) = opt_par_Re(2); % initial enzyme concentration, converts 1 ug/ml per hr (50.5 U/ml)
end
DetoxEff_Re = 1 - T_Re(:,Nr)./T_Re(:,1);
disp('Detox efficiency in 72 hours (Exp):')
disp(DetoxEff_Re)

% infile = 'Filtrate_TP_AFG2_11082019';
% repl = 2:4; % row # on plate
% ns = 8; % column # on plate
ind = 1:2;
T(ind) = T_Re(2:3,1)';
D(:,ind) = [dE_Re(2:3); E0_Re(2:3)];
DE(ind) = DetoxEff_Re(2:3)';
I(ind) = 1;

%% set 2
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
ns = 8; % column # on plate
FL_Re = shiftdim(FL(repl,ns,ntrng)-FL_AFG2_BG,2)';
T_Re = RFUtoAFG2(FL_Re)./(1/T_BG(3)*T_BG(1:length(repl),1:Nr)); % normalized to change in AFG2 without filtrate;
texp = dt*ntrng;

initial_par = [1.5 3];
for cnt = 1:3
    Texp = T_Re(cnt,:);
    opt_par_Re = lsqnonlin(@func_EnzymaticDetox_Filtrate2K,initial_par);
    dE_Re(cnt) = 0.1*opt_par_Re(1); % 0.2, enzyme decay rate, 1/hr
    E0_Re(cnt) = opt_par_Re(2); % initial enzyme concentration, converts 1 ug/ml per hr (50.5 U/ml)
end
DetoxEff_Re = 1 - T_Re(:,Nr)./T_Re(:,1);
disp('Detox efficiency in 72 hours (Exp):')
disp(DetoxEff_Re)

ind = 3:5;
T(ind) = T_Re(1:3,1)';
D(:,ind) = [dE_Re(1:3); E0_Re(1:3)];
DE(ind) = DetoxEff_Re(1:3)';
I(ind) = 2;

%% set 3
infile = 'SizeExclusion_Pseudo_AFG2_12112019';
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
ns = 2; % column # on plate
FL_Re = shiftdim(FL(repl,ns,ntrng)-FL_AFG2_BG,2)';
T_Re = RFUtoAFG2(FL_Re)./(1/T_BG(3)*T_BG(1:length(repl),1:Nr)); % normalized to change in AFG2 without filtrate;
texp = dt*ntrng;

initial_par = [1.5 3];
for cnt = 1:3
    Texp = T_Re(cnt,:);
    opt_par_Re = lsqnonlin(@func_EnzymaticDetox_Filtrate2K,initial_par);
    dE_Re(cnt) = 0.1*opt_par_Re(1); % 0.2, enzyme decay rate, 1/hr
    E0_Re(cnt) = opt_par_Re(2); % initial enzyme concentration, converts 1 ug/ml per hr (50.5 U/ml)
end
DetoxEff_Re = 1 - T_Re(:,Nr)./T_Re(:,1);
disp('Detox efficiency in 72 hours (Exp):')
disp(DetoxEff_Re)

ind = 6:8;
T(ind) = T_Re(1:3,1)';
D(:,ind) = [dE_Re(1:3); E0_Re(1:3)];
DE(ind) = DetoxEff_Re(1:3)';
I(ind) = 3;

%% set 4
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

%% Rpyr starch filtrate AFG2
repl = 5:7; % row # on plate
ns = 2; % column # on plate
FL_Re = shiftdim(FL(repl,ns,ntrng)-FL_AFG2_BG,2)';
T_Re = RFUtoAFG2(FL_Re)./(1/T_BG(3)*T_BG(1:length(repl),1:Nr)); % normalized to change in AFG2 without filtrate;
texp = dt*ntrng;

initial_par = [1.5 3];
for cnt = 1:3
    Texp = T_Re(cnt,:);
    opt_par_Re = lsqnonlin(@func_EnzymaticDetox_Filtrate2K,initial_par);
    dE_Re(cnt) = 0.1*opt_par_Re(1); % 0.2, enzyme decay rate, 1/hr
    E0_Re(cnt) = opt_par_Re(2); % initial enzyme concentration, converts 1 ug/ml per hr (50.5 U/ml)
end
DetoxEff_Re = 1 - T_Re(:,Nr)./T_Re(:,1);
disp('Detox efficiency in 72 hours (Exp):')
disp(DetoxEff_Re)

ind = 9:10;
T(ind) = T_Re(1:2:3,1)';
D(:,ind) = [dE_Re(1:2:3); E0_Re(1:2:3)];
DE(ind) = DetoxEff_Re(1:2:3)';
I(ind) = 4;

figure
plot(T,D(1,:),'sb')
xlabel('AFG2 Conc. (\mug/ml)')
ylabel('Enz. decay rate')

figure
plot(T,DE,'sb')
xlabel('AFG2 Conc. (\mug/ml)')
ylabel('Detox efficiency')

figure
plot(50.5*D(2,:),D(1,:),'sb')
xlabel('Initial enz. conc. (\muU/ml)')
ylabel('Enz. decay rate')
title('Rery, starch')

save(strcat('Fit2K_Filtrate_AFG2_ReryAll_starch_',infile,'ns',num2str(ns),'.mat'))
