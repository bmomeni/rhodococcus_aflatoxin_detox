% Plot degradation comparisons

load('Fit2K_Filtrate_AFG2_RpyrAll_glucose.mat')
TRpG = T;
DRpG = D;
DERpG = DE;

load('Fit2K_Filtrate_AFG2_ReryAll_glucose.mat')
TReG = T;
DReG = D;
DEReG = DE;

load('Fit2K_Filtrate_AFG2_RpyrAll_starch.mat')
TRpS = T;
DRpS = D;
DERpS = DE;

load('Fit2K_Filtrate_AFG2_ReryAll_starch.mat')
TReS = T;
DReS = D;
DEReS = DE;

figure
bar(1,mean(1./DRpG(1,:)),'FaceColor','r')
hold on
errorbar(1,mean(1./DRpG(1,:)),std(1./DRpG(1,:)),'k')
bar(2,mean(1./DReG(1,:)),'FaceColor','b')
errorbar(2,mean(1./DReG(1,:)),std(1./DReG(1,:)),'k')
ylabel('Enz. lifetime (hrs)')

disp('Enz. lifetime comparison on glucose (Rp vs. Re):')
HGR = ranksum(1./DRpG(1,:),1./DReG(1,:)) % lifetime statistical difference? Rp vs Re on glucose
[HGT,pGT] = ttest2(1./DRpG(1,:),1./DReG(1,:)) % lifetime statistical difference? Rp vs Re on glucose

figure
bar(1,50.5*mean(DRpG(2,:)),'FaceColor','r')
hold on
errorbar(1,50.5*mean(DRpG(2,:)),50.5*std(DRpG(2,:)),'k')
bar(2,50.5*mean(DReG(2,:)),'FaceColor','b')
errorbar(2,50.5*mean(DReG(2,:)),50.5*std(DReG(2,:)),'k')
ylabel('Initial enz. conc. (\muU/ml)')

disp('Initial enz. conc. comparison on glucose (Rp vs. Re):')
HGR = ranksum(DRpG(2,:),DReG(2,:)) % initial enz. conc. statistical difference? Rp vs Re on glucose
[HGT,pGT] = ttest2(DRpG(2,:),DReG(2,:)) % initial enz. conc. statistical difference? Rp vs Re on glucose

disp('Enz. lifetime comparison, Rp (glu. vs. starch):')
pTRpR = ranksum(1./DRpG(1,:),1./DRpS(1,:)) % lifetime statistical difference? Rp vs Re on glucose
% HGT = ttest2(1./DRpG(1,:),1./DRpS(1,:)) % lifetime statistical difference? Rp vs Re on glucose

disp('Initia enz. conc. comparison, Rp (glu. vs. starch):')
pCRpR = ranksum(DRpG(2,:),DRpS(2,:)) % lifetime statistical difference? Rp vs Re on glucose
[HCRpT,pCRpT] = ttest2(DRpG(2,:),DRpS(2,:)) % lifetime statistical difference? Rp vs Re on glucose

disp('Enz. lifetime comparison, Re (glu. vs. starch):')
pTReR = ranksum(1./DReG(1,:),1./DReS(1,:)) % lifetime statistical difference? Rp vs Re on glucose
% [HTReT,pTReT] = ttest2(1./DReG(1,:),1./DReS(1,:)) % lifetime statistical difference? Rp vs Re on glucose

disp('Initia enz. conc. comparison, Re (glu. vs. starch):')
pCReR = ranksum(DReG(2,:),DReS(2,:)) % lifetime statistical difference? Rp vs Re on glucose
[HCReT,pCReT] = ttest2(DReG(2,:),DReS(2,:)) % lifetime statistical difference? Rp vs Re on glucose

figure
bar(1,mean(1./DRpG(1,:)),0.5,'FaceColor','r')
hold on
errorbar(1,mean(1./DRpG(1,:)),std(1./DRpG(1,:)),'k')
bar(1.5,mean(1./DRpS(1,:)),0.5,'FaceColor','r')
errorbar(1.5,mean(1./DRpS(1,:)),std(1./DRpS(1,:)),'k')
bar(2.5,mean(1./DReG(1,:)),0.5,'FaceColor','b')
errorbar(2.5,mean(1./DReG(1,:)),std(1./DReG(1,:)),'k')
bar(3,mean(1./DReS(1,:)),0.5,'FaceColor','b')
errorbar(3,mean(1./DReS(1,:)),std(1./DReS(1,:)),'k')
ylabel('Enz. lifetime (hrs)')
ylim([0 16])

figure
bar(1,mean(50.5*DRpG(2,:)),0.5,'FaceColor','r')
hold on
errorbar(1,mean(50.5*DRpG(2,:)),std(50.5*DRpG(2,:)),'k')
bar(1.5,mean(50.5*DRpS(2,:)),0.5,'FaceColor','r')
errorbar(1.5,mean(50.5*DRpS(2,:)),std(50.5*DRpS(2,:)),'k')
bar(2.5,mean(50.5*DReG(2,:)),0.5,'FaceColor','b')
errorbar(2.5,mean(50.5*DReG(2,:)),std(50.5*DReG(2,:)),'k')
bar(3,mean(50.5*DReS(2,:)),0.5,'FaceColor','b')
errorbar(3,mean(50.5*DReS(2,:)),std(50.5*DReS(2,:)),'k')
ylabel('Filtrate enz. conc. (\muU/ml)')
set(gca,'XTick',[1 1.5 2.5 3])
% ylim([0 11])

figure
plot(1./DRpG(1,:),50.5*DRpG(2,:),'ro')
hold on
plot(1./DRpS(1,:),50.5*DRpS(2,:),'rs')
plot(1./DReG(1,:),50.5*DReG(2,:),'bo')
plot(1./DReS(1,:),50.5*DReS(2,:),'bs')
xlabel('Enz. lifetime (hrs)')
ylabel('Filtrate enz. conc. (\muU/ml)')
legend('Rp, glucose','Rp, starch','Re, glucose','Re, starch')
xlim([5 18])
ylim([0 400])

% trade-off
DT = [DRpG, DReG, DReS, DRpS];
figure
plot(1./DT(1,:),50.5*DT(2,:),'o')
xlabel('Enz. lifetime (hrs)')
ylabel('Filtrate enz. conc. (\muU/ml)')

figure
plot(50.5*DT(2,:),1./DT(1,:),'o')
ylabel('Enz. lifetime (hrs)')
xlabel('Filtrate enz. conc. (\muU/ml)')