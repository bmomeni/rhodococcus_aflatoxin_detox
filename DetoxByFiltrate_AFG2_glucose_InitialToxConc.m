clear
%% The issue of multiple acceptable solutions need to be addressed
%% Different pairs of betaT and KE may produce comparable results
% fit2 and func*2 use only two degrees of freedom (dE and kT)

load('Analysis_Background_AFG2_SizeExclusion_Pseudo_AFG2_12112019.mat','T_BG')

E0_Rp0 = [2.1 2.1 2.1 2.1 2.1]; % independently estimated between the two experiments (filtrates not always the same)
E0_Re0 = [3.96 2 2 2 2]; % independently estimated between the two experiments (filtrates not always the same)

%% Rpyr glucose filtrate AFG2
infile = 'DiffInitialToxConc_07132021';
fid = fopen(strcat(infile,'.txt'),'r');
Nr = 865; %number of time points
dt = 5/60; % measurement time-step (hours)
r = 8; % number of rows
c = 12; % number of columns
dg = 2; % data group
FL = ReadDataFromText(infile,Nr,r,c,dg);
FL_AFG2_BG = 262; % estimated from complete degradation cases with live cells

ntrng = 1:Nr; % time points to keep
texp = dt*ntrng;
SlopeRange = 7:18; % range to estimate detox slope, hours 0.5-1.5
%% Rpyr
repl = 5:7; % row # on plate
for ns = 2:5 % column # on plate
    for cnt = 1:length(repl)
        initial_par = [1.5 2];
        FL_Rp = shiftdim(FL(repl(cnt),ns,ntrng)-FL_AFG2_BG,1);
        T_Rp = RFUtoAFG2(FL_Rp)./(1/T_BG(1)*T_BG(1,1:Nr)); % normalized to change in AFG2 without filtrate;
        DetoxEff_Rp(ns,cnt) = 1 - T_Rp(Nr)./T_Rp(1);
        InitT_Rp(ns,cnt) = mean(T_Rp(1:3));
        p = polyfit(texp(SlopeRange),T_Rp(SlopeRange),1);
        DetoxSlope_Rp(ns,cnt) = p(1);
        DetoxMod_Rp(ns,cnt) = -1/E0_Rp0(ns)*DetoxSlope_Rp(ns,cnt);
        figure(11)
        hold on
        plot(texp,T_Rp,'color',(1-(ns-1)/10)*[1 0 0])
        xlabel('Time (hrs)')
        ylabel('AFG2 conc. (\mug/ml)')
    end
end

%% Rery
repl = 2:4; % row # on plate
for ns = 2:5 % column # on plate
    for cnt = 1:length(repl)
        initial_par = [1.5 2];
        FL_Re = shiftdim(FL(repl(cnt),ns,ntrng)-FL_AFG2_BG,1);
        T_Re = RFUtoAFG2(FL_Re)./(1/T_BG(1)*T_BG(1,1:Nr)); % normalized to change in AFG2 without filtrate;
        DetoxEff_Re(ns,cnt) = 1 - T_Re(Nr)./T_Re(1);
        InitT_Re(ns,cnt) = mean(T_Re(1:3));
        p = polyfit(texp(SlopeRange),T_Re(SlopeRange),1);
        DetoxSlope_Re(ns,cnt) = p(1);
        DetoxMod_Re(ns,cnt) = -1/E0_Re0(ns)*DetoxSlope_Re(ns,cnt);
        figure(12)
        hold on
        plot(texp,T_Re,'color',(1-(ns-1)/10)*[0 0 1])
        xlabel('Time (hrs)')
        ylabel('AFG2 conc. (\mug/ml)')
    end
end
fclose(fid);

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

%% Rpyr filtrate AFG2, T ~ 30
repl = 2:4; % row # on plate
ns = 3; % column # on plate
% FL_Rp = shiftdim(FL(repl,ns,ntrng)-FL_AFG2_BG,2)';
% T_Rp = RFUtoAFG2(FL_Rp)./(1/T_BG(1)*T_BG); % normalized to change in AFG2 without filtrate
for cnt = 1:length(repl)
    initial_par = [1.5 2];
    FL_Rp = shiftdim(FL(repl(cnt),ns,ntrng)-FL_AFG2_BG,1);
    T_Rp = RFUtoAFG2(FL_Rp)./(1/T_BG(1)*T_BG(1,1:Nr)); % normalized to change in AFG2 without filtrate;
    DetoxEff_Rp(1,cnt) = 1 - T_Rp(Nr)./T_Rp(1);
    InitT_Rp(1,cnt) = mean(T_Rp(1:3));
    p = polyfit(texp(SlopeRange),T_Rp(SlopeRange),1);
    DetoxSlope_Rp(1,cnt) = p(1);
    DetoxMod_Rp(1,cnt) = -1/E0_Rp0(1)*DetoxSlope_Rp(1,cnt);
    figure(11)
    hold on
    plot(texp,T_Rp,'color',(1-(ns-1)/10)*[1 0 0])
    xlabel('Time (hrs)')
    ylabel('AFG2 conc. (\mug/ml)')
end

%% Rery filtrate AFG2, T ~ 30
repl = 2:4; % row # on plate
ns = 6; % column # on plate
% FL_Re = shiftdim(FL(repl,ns,ntrng)-FL_AFG2_BG,2)';
% T_Re = RFUtoAFG2(FL_Re)./(1/T_BG(3)*T_BG); % normalized to change in AFG2 without filtrate;
for cnt = 1:length(repl)
    initial_par = [1.5 2];
    FL_Re = shiftdim(FL(repl(cnt),ns,ntrng)-FL_AFG2_BG,1);
    T_Re = RFUtoAFG2(FL_Re)./(1/T_BG(1)*T_BG(1,1:Nr)); % normalized to change in AFG2 without filtrate;
    T_Re = T_Re - 1.8*(cnt==3); % adjusted for systematic background shift
    DetoxEff_Re(1,cnt) = 1 - T_Re(Nr)./T_Re(1);
    InitT_Re(1,cnt) = mean(T_Re(1:3));
    p = polyfit(texp(SlopeRange),T_Re(SlopeRange),1);
    DetoxSlope_Re(1,cnt) = p(1);
    DetoxMod_Re(1,cnt) = -1/E0_Re0(1)*DetoxSlope_Re(1,cnt);
    figure(12)
    hold on
    plot(texp,T_Re,'color',(1-(ns-1)/10)*[0 0 1])
    xlabel('Time (hrs)')
    ylabel('AFG2 conc. (\mug/ml)')
end
fclose(fid);

TK = linspace(0,30,200);
KRp = 10;
yRp = TK./(TK+KRp);
figure
plot(InitT_Rp,DetoxMod_Rp,'ro')
hold on
plot(TK,yRp,'k:')
xlabel('Initial AFG_2 conc. (\mug/ml)')
ylabel('T_0/(K_m+T_0)')

% figure
% plot(InitT_Rp,DetoxSlope_Rp,'ro')
% xlabel('Initial AFG_2 conc. (\mug/ml)')
% ylabel('Detox slope, Rp')

TK = linspace(0,30,200);
KRe = 12;
yRe = TK./(TK+KRe);
figure
plot(InitT_Re,DetoxMod_Re,'bo')
hold on
plot(TK,yRe,'k:')
xlabel('Initial AFG_2 conc. (\mug/ml)')
ylabel('T_0/(K_m+T_0)')

% figure
% plot(InitT_Re,DetoxSlope_Re,'bo')
% xlabel('Initial AFG_2 conc. (\mug/ml)')
% ylabel('Detox slope (Re)')

save(strcat('DetoxByFiltrate_AFG2_glucose_InitialToxConc_',infile,'.mat'))
