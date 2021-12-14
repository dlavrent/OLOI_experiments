%% plots ORNs
clear all
close all

glomsToPlot{1}='DC2';
glomsToPlot{2}='DL5';
glomsToPlot{3}='DM1';
glomsToPlot{4}='DM2';
glomsToPlot{5}='DM3';

load('ornexample_volumes.mat')
load('ornexample_volumes_manualLabels.mat')

todelete=zeros(1,length(clusterLabels));
glomorder=zeros(1,length(glomsToPlot));
for i=1:length(clusterLabels)
    glomok=0;
   for j=1:length(glomsToPlot)
       if strcmp(clusterLabels{i},glomsToPlot{j})
           glomok=1;
           glomorder(j)=i;
          break 
       end
   end
   
   if ~glomok
       todelete(i)=1;
   end
end
deleter=find(todelete);
clusterLabels(deleter)=[];
clusterVolU(deleter)=[];
clusterInfoU(deleter)=[];
grnResponse(:,deleter,:)=[];

[nada inds]=sort(glomorder);

% reorder gloms so they are in correct order
clusterLabels=clusterLabels(inds);
clusterVolU=clusterVolU(inds);
clusterInfoU=clusterInfoU(inds);
grnResponse=grnResponse(:,inds,:);


% FIG 1e
% figure 1
showClusters(clusterVolU,clusterInfoU,clusterLabels)
axis off

grnResponseO=grnResponse;


% PNs
load('pnexample_volumes.mat')
load('pnexample_volumes_manualLabels.mat')


todelete=zeros(1,length(clusterLabels));
glomorder=zeros(1,length(glomsToPlot));
for i=1:length(clusterLabels)
    glomok=0;
   for j=1:length(glomsToPlot)
       if strcmp(clusterLabels{i},glomsToPlot{j})
           glomok=1;
           glomorder(j)=i;
          break 
       end
   end
   
   if ~glomok
       todelete(i)=1;
   end
end
deleter=find(todelete);
clusterLabels(deleter)=[];
clusterVolU(deleter)=[];
clusterInfoU(deleter)=[];
grnResponse(:,deleter,:)=[];

[nada inds]=sort(glomorder);

% reorder gloms so they are in correct order
clusterLabels=clusterLabels(inds);
clusterVolU=clusterVolU(inds);
clusterInfoU=clusterInfoU(inds);
grnResponse=grnResponse(:,inds,:);

% FIG 1e
% figure 2
showClusters(clusterVolU,clusterInfoU,clusterLabels)
axis off

for i=1:13
    for j=1:length(glomsToPlot)
   grnResponseO(i,j,:)=grnResponseO(i,j,:)-mean(grnResponseO(i,j,1:3)); 
   grnResponse(i,j,:)=grnResponse(i,j,:)-mean(grnResponse(i,j,1:3)); 
    end
end

tofill=[2 4 11 9 7];
for i=1:length(tofill)
grnResponseO(tofill(i),:,19:20)=grnResponseO(tofill(i),:,17:18);
grnResponse(tofill(i),:,19:20)=grnResponse(tofill(i),:,17:18);
end

pncolor=[0 0.8 0];
orncolor=[0.9 0 0.7];

figure %3
plot(1.2*(1:20),squeeze(nanmean(grnResponseO(:,1,:),1)),'Color',orncolor,'LineWidth',3)
hold on
plot(1.2*(1:20),squeeze(nanmean(grnResponse(:,1,:),1)),'Color',pncolor,'LineWidth',3)
xlabel('time (s)')
ylabel('dF/F')
legend('ORNs','PNs')
legend boxoff
box off
set(gca,'FontSize',15)


figure %4
plot(1.2*(1:20),squeeze(grnResponse(2,1,:)),'Color',[0.4 0.3 0.8],'LineWidth',3)
hold on
plot(1.2*(1:20),squeeze(grnResponse(11,1,:)),'--','Color',[0.8 0.3 0.2],'LineWidth',3)
xlabel('time (s)')
ylabel('dF/F')
legend('3-octanol','4-methylcyclohexanol')
legend boxoff
box off
set(gca,'FontSize',15)


figure %5
plot(1.2*(1:20),squeeze(grnResponse(6,1,:)),'Color',[0.8 0.4 0.2],'LineWidth',3)
hold on
plot(1.2*(1:20),squeeze(grnResponse(6,4,:)),'--','Color',[0.3 0.8 0.4],'LineWidth',3)
xlabel('time (s)')
ylabel('dF/F')
legend('DC2','DM2')
legend boxoff
box off
set(gca,'FontSize',15)

