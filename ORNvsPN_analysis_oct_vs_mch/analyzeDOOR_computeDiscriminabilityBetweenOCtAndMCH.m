%% analyze DoOR data. Compute Discriminability between OCT and MCH for all glomerulus pairs
clear all
close all

publishedOdorPath='/Users/mattchurgin/Desktop/dblab/mattFunctions/odorpanelResponsesDoOR/odorPanel_12/odorPanel_12_DoORData.mat';
load(publishedOdorPath);
publishedORNresponse=publishedOR.gh146response;
gh146glomnames=publishedOR.gh146glomerulusNames;

oct=publishedORNresponse(:,1);
mch=publishedORNresponse(:,10);

octmchdiff=(oct-mch);

discrim=zeros(length(octmchdiff),length(octmchdiff));
for i = 1:length(octmchdiff)
    for j = 1:length(octmchdiff)
        discrim(i,j)=abs(octmchdiff(i)-octmchdiff(j));
    end
end

figure;
imagesc(discrim)
set(gca,'xtick',1:length(octmchdiff),'xticklabel',[string(gh146glomnames)],'FontSize',5)
xtickangle(90)
set(gca,'ytick',1:length(octmchdiff),'yticklabel',[string(gh146glomnames)])
colorbar
title('OCT-MCH discriminability')
xlabel('glomerulus','FontSize',15)
ylabel('glomerulus','FontSize',15)

figure;
[y1 x1]=hist(discrim(:),20);
plot(x1,y1,'k','LineWidth',3)
hold on
% draw line where DC2-DM2 lies
plot(discrim(6,13)*ones(1,5),linspace(0,max(y1),5),'k--','LineWidth',2);
text(discrim(6,13),max(y1),['DM2-DC2'],'FontSize',10)
xlabel('discriminability')
ylabel('counts')
box off
set(gca,'FontSize',15)