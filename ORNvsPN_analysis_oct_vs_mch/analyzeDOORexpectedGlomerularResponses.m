%% analyze DoOR data to determine likely glomeruli to extract
clear all
close all

publishedOdorPath='/Users/mattchurgin/Desktop/dblab/mattFunctions/odorpanelResponsesDoOR/odorPanel_12/odorPanel_12_DoORData.mat';
load(publishedOdorPath);
publishedORNresponse=publishedOR.gh146response;
 
% sort gh146 by mean response to odor panel
meanPublishedResponse=nanmean(publishedORNresponse,2);
[ix iy]=sort(meanPublishedResponse,'descend');
iy(isnan(ix))=[];
ix(isnan(ix))=[];

sortedgh146glomnames=publishedOR.gh146glomerulusNames(iy);

figure
plot(1:length(iy),ix,'ko','LineWidth',3)
ylabel('mean activation to odor panel (DoOR)')
set(gca,'xtick',1:length(iy),'xticklabel',[string(sortedgh146glomnames)])
xtickangle(30)
box off

figure
% sort gh146 by mean response to octanol
octPublishedResponse=publishedORNresponse(:,1);
[ix iy]=sort(octPublishedResponse,'descend');
iy(isnan(ix))=[];
ix(isnan(ix))=[];

octsortedgh146glomnames=publishedOR.gh146glomerulusNames(iy);

subplot(1,2,1)
plot(1:length(iy),ix,'ko','LineWidth',3)
ylabel('OCT activation (DoOR)')
set(gca,'FontSize',10)
set(gca,'xtick',1:length(iy),'xticklabel',[string(octsortedgh146glomnames)])
xtickangle(30)
box off

% sort gh146 by mean response to MCH
mchPublishedResponse=publishedORNresponse(:,10);
[ix iy]=sort(mchPublishedResponse,'descend');
iy(isnan(ix))=[];
ix(isnan(ix))=[];

mchsortedgh146glomnames=publishedOR.gh146glomerulusNames(iy);

subplot(1,2,2)
plot(1:length(iy),ix,'ko','LineWidth',3)
ylabel('MCH activation (DoOR)')
set(gca,'xtick',1:length(iy),'xticklabel',[string(mchsortedgh146glomnames)])
xtickangle(30)
set(gca,'FontSize',10)
box off


%% now for orco
publishedORNresponse=publishedOR.orcoresponse;
% sort orco by mean response to odor panel
meanPublishedResponse=nanmean(publishedORNresponse,2);
[ix iy]=sort(meanPublishedResponse,'descend');
iy(isnan(ix))=[];
ix(isnan(ix))=[];

sortedorcoglomnames=publishedOR.orcoGlomerulusNames(iy);

figure
plot(1:length(iy),ix,'ko','LineWidth',3)
ylabel('mean activation to odor panel (DoOR)')
set(gca,'xtick',1:length(iy),'xticklabel',[string(sortedorcoglomnames)])
xtickangle(30)
box off

figure
% sort orco by mean response to octanol
octPublishedResponse=publishedORNresponse(:,1);
[ix iy]=sort(octPublishedResponse,'descend');
iy(isnan(ix))=[];
ix(isnan(ix))=[];

octsortedorcoglomnames=publishedOR.orcoGlomerulusNames(iy);

subplot(1,2,1)
plot(1:length(iy),ix,'ko','LineWidth',3)
ylabel('OCT activation (DoOR)')
set(gca,'FontSize',10)
set(gca,'xtick',1:length(iy),'xticklabel',[string(octsortedorcoglomnames)])
xtickangle(30)
box off

% sort orco by mean response to MCH
mchPublishedResponse=publishedORNresponse(:,10);
[ix iy]=sort(mchPublishedResponse,'descend');
iy(isnan(ix))=[];
ix(isnan(ix))=[];

mchsortedorcoglomnames=publishedOR.orcoGlomerulusNames(iy);

subplot(1,2,2)
plot(1:length(iy),ix,'ko','LineWidth',3)
ylabel('MCH activation (DoOR)')
set(gca,'xtick',1:length(iy),'xticklabel',[string(mchsortedorcoglomnames)])
xtickangle(30)
set(gca,'FontSize',10)
box off
