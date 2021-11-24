% Compare PN and ORN data
% use PN "test set" flies acquired during similar time period as ORNs

% load PNs
clear all
close all

allColorMaps

pncolor=[0 0.8 0];
orncolor=[0.9 0 0.7];
rng('default')

manualLabelHome='/Users/mattchurgin/Dropbox (Harvard University)/flyimaging/analysis/ORNvsPN_analysis_oct_vs_air/PN_alldata';

publishedOdorPath='/Users/mattchurgin/Desktop/dblab/mattFunctions/odorpanelResponsesDoOR/odorPanel_12/odorPanel_12_DoORData.mat';
load(publishedOdorPath);

manualLabelledFolders=dir(manualLabelHome);
manualLabelledFolders=manualLabelledFolders(3:end);

labels=cell(1,length(manualLabelledFolders));
glomsannotated=zeros(1,length(manualLabelledFolders));
totalPutativeGloms=zeros(1,length(manualLabelledFolders));

nodors=13;
odortimes=[6:9]; % hard-coded specific time interval for summarizing odor response
odortimesfortimecourse=[1:18];
responsesTensorFullyOrthogonal=NaN*zeros(22,2,39,nodors,length(odortimesfortimecourse)); % flynum HARD CODED
responsesGlomByOdor=[];
responsesTimeCourseGlomByOdor=[];
flyodors=[];
flygloms=[];
glombyodorflies=[];
behaviorOcc=[];
behaviorpreOcc=[];
leftRightSplitVal=1000;
firstSecondPanelSplitVal=500;
flyNum=0;
glomfound=zeros(1,39);

oldFlyNum=9000;
oldDate=[];
oldLobe=[];

for i=1:length(manualLabelledFolders)
    currname=manualLabelledFolders(i).name;
    if ~strcmp(currname(1),'.')
        if strcmp(currname(end),'1')
            currVol='Volumes';
        else
            currVol='Volumes2';
        end
        if strcmp(currname(end-2),'r')
            currLobe='rightLobe';
            cL=2;
        else
            currLobe='leftLobe';
            cL=1;
        end
        if strcmp(currname(1),'H')
            currDate=currname(8:13);
        else
            currDate=currname(1:6);
        end
        underscores=strfind(currname,'_');
        currFlyNum=currname((underscores(1)+4):(underscores(2)-1));
        
        if ~strcmp(currDate,oldDate) || ~strcmp(currFlyNum,oldFlyNum)
            flyNum=flyNum+1;
        end
        
        currFiles=dir([manualLabelHome '/' manualLabelledFolders(i).name]);
        currFiles=currFiles(3:end);
        for j=1:length(currFiles)
            cfcname=currFiles(j).name;
            if ~strcmp(cfcname(1),'.')
                % load manual label file, response data, and behavior data
                load([manualLabelHome '/' manualLabelledFolders(i).name '/' currFiles(j).name])
            end
        end
        
        %behaviorOcc=[behaviorOcc occ_zeroed];
        %behaviorpreOcc=[behaviorpreOcc preocc_zeroed];
        
        behaviorOcc=[behaviorOcc occ-preocc];
        behaviorpreOcc=[behaviorpreOcc preocc];
        
        
        manualClusterLabels=clusterLabels;
        totalPutativeGloms(i)=length(manualClusterLabels);
        
        
        gs=median(grnResponse(:,:,odortimes),3); % use median
         gs=zeros(13,size(grnResponse,2));
        for j=1:size(grnResponse,2)
            temp = squeeze(grnResponse(:,j,odortimes)-nanmedian(grnResponse(:,j,1:5),3));
            gs(:,j)= transpose(max(temp'));
        end
        %gs=prctile(grnResponse(:,:,odortimes),95,3); % use percentile
        responseTemp=NaN*zeros(length(publishedOR.gh146glomerulusNames),nodors);
        responseTempTimeCourse=NaN*zeros(length(publishedOR.gh146glomerulusNames)*length(odortimesfortimecourse),size(grnResponse,1));
        responsesTensorTemp=NaN*zeros(length(publishedOR.gh146glomerulusNames),size(grnResponse,1),length(odortimesfortimecourse));
        labels{i}.manual=cell(1,length(manualClusterLabels));
        % record manual labels
        for j=1:length(manualClusterLabels)
            labels{i}.manual{j}=manualClusterLabels{j};
        end
        
        % fill in response matrix
        for j=1:length(manualClusterLabels)
            for k=1:length(publishedOR.gh146glomerulusNames)
                if strcmp(labels{i}.manual{j},publishedOR.gh146glomerulusNames{k})
                    responseTemp(k,:)=gs(:,j);
                    glomsannotated(i)=glomsannotated(i)+1;
                    glomfound(k)=glomfound(k)+1;
                    for oo=1:nodors
                        responseTempTimeCourse((((k-1)*length(odortimesfortimecourse)+1):k*length(odortimesfortimecourse)),oo)=grnResponse(oo,j,odortimesfortimecourse);
                        responsesTensorTemp(k,oo,:)=grnResponse(oo,j,odortimesfortimecourse);
                    end
                    break
                end
            end
        end
        
        responseTempT=responseTemp';
        responseTempTimeCourseT=responseTempTimeCourse';
        
        responsesTensorFullyOrthogonal(flyNum,cL,:,:,:)=responsesTensorTemp;
        
        responsesGlomByOdor=[responsesGlomByOdor responseTempT(:)];
        responsesTimeCourseGlomByOdor=[responsesTimeCourseGlomByOdor responseTempTimeCourseT(:)];
        glombyodorflies=[glombyodorflies (flyNum+(cL-1)*leftRightSplitVal)];
        
        oldFlyNum=currFlyNum;
        oldDate=currDate;
        oldLobe=currLobe;
    end
end

% % remove antenna lobes with fewer than X glomeruli annotated
% glomthreshold=5;
% toremove=find(glomsannotated<=glomthreshold);
% responsesGlomByOdor(:,toremove)=[];
% responsesTimeCourseGlomByOdor(:,toremove)=[];
% glombyodorflies(toremove)=[];

% data points that correspond to each fly
flyindices=cell(1,flyNum);
flyindicesL=cell(1,flyNum);
flyindicesR=cell(1,flyNum);
for i=1:flyNum
    [temp temp2]=find(glombyodorflies==i);
    [tem3 temp4]=find(glombyodorflies==(i+leftRightSplitVal));
    flyindices{i}=[temp2 temp4];
    flyindicesL{i}=temp2;
    flyindicesR{i}=temp4;
    
    
end

mycmap=distinguishable_colors(flyNum);

% figure;
% for i=1:length(glomsannotated)
%     plot([1 2],[glomsannotated(i) totalPutativeGloms(i)],'Marker','o','LineStyle','--','LineWidth',2)
%     hold on
% end
% set(gca,'XTick',[1:2])
% set(gca,'XTickLabel',[{'Annotated'},{'Available'}])
% axis([0.5 2.5 0 40])
% ylabel('Glomeruli')
% box off
% set(gca,'FontSize',15)

% perform pca on PN responses

clear responsesNoResponseRemoved

fracIn=0.25;% best results when fracIn is high, ~0.5, only using high confidence glomeruli

responsesNoResponseRemoved=responsesGlomByOdor;
responsesNoResponseRemovedT=responsesTimeCourseGlomByOdor;
responsesTensor=responsesTensorFullyOrthogonal;

gNames=publishedOR.gh146glomerulusNames;
glomsFound=glomfound;
numFinite=sum(isfinite(responsesNoResponseRemoved),2);
toRemove=find(numFinite/size(responsesNoResponseRemoved,2)<=fracIn);
responsesNoResponseRemoved(toRemove,:)=[];
numFiniteT=sum(isfinite(responsesNoResponseRemovedT),2);
toRemoveT=find(numFiniteT/size(responsesNoResponseRemovedT,2)<=fracIn);
responsesNoResponseRemovedT(toRemoveT,:)=[];

temp=nodors;
fp=toRemove(find(mod(toRemove,temp)==1));
glomsremoved=((fp-1)/temp)+1;
gNames(glomsremoved)=[];
glomsFound(glomsremoved)=[];
responsesTensor(:,:,glomsremoved,:,:)=[];


% manually remove D
gNames(1)=[];
responsesNoResponseRemoved(1:nodors,:)=[];
responsesNoResponseRemovedT(1:((nodors)*length(odortimesfortimecourse)),:)=[];
responsesTensor(:,:,1,:,:)=[];

gh146rawdata=responsesNoResponseRemoved;

% % fill nans with mean
for i=1:size(responsesNoResponseRemoved,1)
    for j=1:size(responsesNoResponseRemoved,2)
        if isnan(responsesNoResponseRemoved(i,j))
            responsesNoResponseRemoved(i,j)=nanmean(responsesNoResponseRemoved(i,:));
        end
    end
end

% fill missing values with linear interpolation
% data=responsesNoResponseRemoved';
% dz=data;
% dataFilled=fillWithRegressedValues(dz);
% responsesNoResponseRemoved=dataFilled';

yesZscore=0;
if yesZscore
    responsesNoResponseRemoved=responsesNoResponseRemoved';
    responsesNoResponseRemoved=(responsesNoResponseRemoved-mean(responsesNoResponseRemoved))./std(responsesNoResponseRemoved);
    responsesNoResponseRemoved=responsesNoResponseRemoved';
end
% load ORNs

manualLabelHome='/Users/mattchurgin/Dropbox (Harvard University)/flyimaging/analysis/ORNvsPN_analysis_oct_vs_air/ORN_alldata';

publishedOdorPath='/Users/mattchurgin/Desktop/dblab/mattFunctions/odorpanelResponsesDoOR/odorPanel_12/odorPanel_12_DoORData.mat';
load(publishedOdorPath);
manualLabelledFolders=dir(manualLabelHome);
manualLabelledFolders=manualLabelledFolders(3:end);

labels=cell(1,length(manualLabelledFolders));
glomsannotatedO=zeros(1,length(manualLabelledFolders));
totalPutativeGloms=zeros(1,length(manualLabelledFolders));

nodors=13;
odortimes=[6:9]; % hard-coded specific time interval for summarizing odor response
odortimesfortimecourse=[1:18];
responsesGlomByOdorO=[];
responsesTimeCourseGlomByOdorO=[];
responsesTensorFullyOrthogonalO=NaN*zeros(30,2,39,nodors,length(odortimesfortimecourse)); % flynum HARD CODED
flyodors=[];
flygloms=[];
glombyodorfliesO=[];
behaviorOccO=[];
behaviorpreOccO=[];
leftRightSplitVal=1000;
firstSecondPanelSplitVal=500;
flyNumO=0;
glomfoundO=zeros(1,39);

oldFlyNum=9000;
oldDate=[];
oldLobe=[];

for i=1:length(manualLabelledFolders)
    currname=manualLabelledFolders(i).name;
    if ~strcmp(currname(1),'.')
        if strcmp(currname(end),'1')
            currVol='Volumes';
        else
            currVol='Volumes2';
        end
        if strcmp(currname(end-2),'r')
            currLobe='rightLobe';
            cL=2;
        else
            currLobe='leftLobe';
            cL=1;
        end
        if strcmp(currname(1),'H')
            currDate=currname(8:13);
        else
            currDate=currname(1:6);
        end
        underscores=strfind(currname,'_');
        currFlyNum=currname((underscores(1)+4):(underscores(2)-1));
        
        if ~strcmp(currDate,oldDate) || ~strcmp(currFlyNum,oldFlyNum)
            flyNumO=flyNumO+1;
        end
        
        currFiles=dir([manualLabelHome '/' manualLabelledFolders(i).name]);
        currFiles=currFiles(3:end);
        for j=1:length(currFiles)
            cfcname=currFiles(j).name;
            if ~strcmp(cfcname(1),'.')
                % load manual label file, response data, and behavior data
                load([manualLabelHome '/' manualLabelledFolders(i).name '/' currFiles(j).name])
            end
        end
        
        %behaviorOcc=[behaviorOcc occ_zeroed];
        %behaviorpreOcc=[behaviorpreOcc preocc_zeroed];
        behaviorOccO=[behaviorOccO occ-preocc];
        %behaviorOcc=[behaviorOcc occ];
        behaviorpreOccO=[behaviorpreOccO preocc];
        
        manualClusterLabels=clusterLabels;
        totalPutativeGloms(i)=length(manualClusterLabels);
        
        
        gs=median(grnResponse(:,:,odortimes),3); % use median
         gs=zeros(13,size(grnResponse,2));
        for j=1:size(grnResponse,2)
            temp = squeeze(grnResponse(:,j,odortimes)-nanmedian(grnResponse(:,j,1:5),3));
            gs(:,j)= transpose(max(temp'));
        end
        %gs=prctile(grnResponse(:,:,odortimes),95,3); % use percentile
        responseTemp=NaN*zeros(length(publishedOR.gh146glomerulusNames),nodors);
        responseTempTimeCourse=NaN*zeros(length(publishedOR.gh146glomerulusNames)*length(odortimesfortimecourse),size(grnResponse,1));
        responsesTensorTemp=NaN*zeros(length(publishedOR.gh146glomerulusNames),size(grnResponse,1),length(odortimesfortimecourse));
        labels{i}.manual=cell(1,length(manualClusterLabels));
        % record manual labels
        for j=1:length(manualClusterLabels)
            labels{i}.manual{j}=manualClusterLabels{j};
        end
        
        % fill in response matrix
        for j=1:length(manualClusterLabels)
            for k=1:length(publishedOR.gh146glomerulusNames)
                if strcmp(labels{i}.manual{j},publishedOR.gh146glomerulusNames{k})
                    responseTemp(k,:)=gs(:,j);
                    glomsannotatedO(i)=glomsannotatedO(i)+1;
                    glomfoundO(k)=glomfoundO(k)+1;
                    for oo=1:nodors
                        responseTempTimeCourse((((k-1)*length(odortimesfortimecourse)+1):k*length(odortimesfortimecourse)),oo)=grnResponse(oo,j,odortimesfortimecourse);
                        responsesTensorTemp(k,oo,:)=grnResponse(oo,j,odortimesfortimecourse);
                    end
                    break
                end
            end
        end
        
        responsesTensorFullyOrthogonalO(flyNumO,cL,:,:,:)=responsesTensorTemp;
        
        responseTempT=responseTemp';
        responseTempTimeCourseT=responseTempTimeCourse';
        
        responsesGlomByOdorO=[responsesGlomByOdorO responseTempT(:)];
        responsesTimeCourseGlomByOdorO=[responsesTimeCourseGlomByOdorO responseTempTimeCourseT(:)];
        glombyodorfliesO=[glombyodorfliesO (flyNumO+(cL-1)*leftRightSplitVal)];
        
        oldFlyNum=currFlyNum;
        oldDate=currDate;
        oldLobe=currLobe;
    end
end

% % remove antenna lobes with fewer than X glomeruli annotated
% glomthreshold=5;
% toremove=find(glomsannotated<=glomthreshold);
% responsesGlomByOdor(:,toremove)=[];
% responsesTimeCourseGlomByOdor(:,toremove)=[];
% glombyodorflies(toremove)=[];

% data points that correspond to each fly
flyindicesO=cell(1,flyNumO);
flyindicesOL=cell(1,flyNumO);
flyindicesOR=cell(1,flyNumO);
for i=1:flyNumO
    [temp temp2]=find(glombyodorfliesO==i);
    [tem3 temp4]=find(glombyodorfliesO==(i+leftRightSplitVal));
    flyindicesO{i}=[temp2 temp4];
    flyindicesOL{i}=temp2;
    flyindicesOR{i}=temp4;
end


mycmapO=distinguishable_colors(flyNumO);

% perform pca on responses on ORNs

clear responsesNoResponseRemovedO

fracIn=0.3; % best results when fracIn is high, ~0.5, only using high confidence glomeruli

responsesNoResponseRemovedO=responsesGlomByOdorO;
responsesNoResponseRemovedTO=responsesTimeCourseGlomByOdorO;
responsesTensorO=responsesTensorFullyOrthogonalO;

% remove fully empty rows
%totallyempty=~any(isfinite(responsesNoResponseRemoved),2);
%responsesNoResponseRemoved(totallyempty,:)=[];

gNamesO=publishedOR.gh146glomerulusNames;
glomsFoundO=glomfoundO;
numFinite=sum(isfinite(responsesNoResponseRemovedO),2);
toRemove=find(numFinite/size(responsesNoResponseRemovedO,2)<=fracIn);
responsesNoResponseRemovedO(toRemove,:)=[];
numFiniteT=sum(isfinite(responsesNoResponseRemovedTO),2);
toRemoveT=find(numFiniteT/size(responsesNoResponseRemovedTO,2)<=fracIn);
responsesNoResponseRemovedTO(toRemoveT,:)=[];

temp=nodors;
fp=toRemove(find(mod(toRemove,temp)==1));
glomsremoved=((fp-1)/temp)+1;
gNamesO(glomsremoved)=[];
glomsFoundO(glomsremoved)=[];
responsesTensorO(:,:,glomsremoved,:,:)=[];

orcoRawData=responsesNoResponseRemovedO;

% % fill nans with mean
for i=1:size(responsesNoResponseRemovedO,1)
    for j=1:size(responsesNoResponseRemovedO,2)
        if isnan(responsesNoResponseRemovedO(i,j))
            responsesNoResponseRemovedO(i,j)=nanmean(responsesNoResponseRemovedO(i,:));
        end
    end
end

% remove air and ethanol
%responsesNoResponseRemoved(1:nodors:end,:)=0;
%responsesNoResponseRemoved(8:nodors:end,:)=0;

% fill missing values with linear interpolation
% data=responsesNoResponseRemoved';
% dz=data;
% dataFilled=fillWithRegressedValues(dz);
% responsesNoResponseRemoved=dataFilled';

yesZscore=0;
if yesZscore
    responsesNoResponseRemovedO=responsesNoResponseRemovedO';
    responsesNoResponseRemovedO=(responsesNoResponseRemovedO-mean(responsesNoResponseRemovedO))./std(responsesNoResponseRemovedO);
    responsesNoResponseRemovedO=responsesNoResponseRemovedO';
end


% plot all fly x trials
figure;
subplot(1,2,1)
imagesc(gh146rawdata,[-0.5 3])
xlabel('fly x trial')
ylabel('glom x odor')
set(gca,'ytick',6:nodors:size(gh146rawdata,1))
set(gca,'yticklabel',[string(gNames)])
ytickangle(30)
title('PNs')
axis tight
set(gca,'FontSize',15)

subplot(1,2,2)
imagesc(orcoRawData,[-0.5 3])
xlabel('fly x trial')
title('ORNs')
set(gca,'ytick',6:nodors:size(gh146rawdata,1))
set(gca,'yticklabel',[''])
set(gca,'FontSize',15)
axis tight
hcb=colorbar;
title(hcb,'dF/F')


% plot average response for each fly over all trials
gh146flyaverage=zeros(size(gh146rawdata,1),length(flyindices));
orcoflyaverage=zeros(size(orcoRawData,1),length(flyindicesO));
for i=1:length(flyindices)
    gh146flyaverage(:,i)=nanmean(gh146rawdata(:,flyindices{i}),2);
end
for i=1:length(flyindicesO)
    orcoflyaverage(:,i)=nanmean(orcoRawData(:,flyindicesO{i}),2);
end
for i=1:size(gh146flyaverage,1)
    temp=find(isnan(gh146flyaverage(i,:)));
    gh146flyaverage(i,temp)=nanmean(gh146flyaverage(i,:));
end
for i=1:size(orcoflyaverage,1)
    temp=find(isnan(orcoflyaverage(i,:)));
    orcoflyaverage(i,temp)=nanmean(orcoflyaverage(i,:));
end

ctrst=[0 1.75];

figure;
subplot(1,2,1)
imagesc(gh146flyaverage,ctrst)
colormap(hot)
%colormap(cm.egoalley)
xlabel('fly #')
ylabel('glom x odor')
set(gca,'ytick',6:nodors:size(gh146flyaverage,1))
set(gca,'yticklabel',[string(gNames)])
ytickangle(30)
title('PNs')
axis tight
set(gca,'FontSize',15)

subplot(1,2,2)
imagesc(orcoflyaverage,ctrst)
%colormap(cm.egoalley)
colormap(hot)
xlabel('fly #')
title('ORNs')
set(gca,'ytick',6:nodors:size(orcoflyaverage,1))
set(gca,'yticklabel',[''])
set(gca,'FontSize',15)
axis tight
hcb=colorbar;
title(hcb,'dF/F')

% plot GH146 vs Orco activation
gh=nanmean(gh146rawdata,2);
orc=nanmean(orcoRawData,2);
figure
plot(orc,gh,'k.','LineWidth',3,'MarkerSize',15)
xlabel('ORN activation (dF/F)')
ylabel('PN activation (dF/F)')
set(gca,'FontSize',15)
box off

% show heat map of average odor response
ghmatrix=transpose(reshape(gh,nodors,length(gNames)));
ormatrix=transpose(reshape(orc,nodors,length(gNamesO)));

odorNames{1}='air';
odorNames{2}='3-octanol';
odorNames{3}='1-hexanol';
odorNames{4}='ethyl lactate';
odorNames{5}='citronella';
odorNames{6}='2-heptanone';
odorNames{7}='1-pentanol';
odorNames{8}='ethanol';
odorNames{9}='geranyl acetate';
odorNames{10}='hexyl acetate';
odorNames{11}='4-methylcyclohexanol';
odorNames{12}='pentyl acetate';
odorNames{13}='1-butanol';


ctrst=[0 1.65];

figure;
imagesc(ghmatrix,ctrst)
%colormap(cm.egoalley)
colormap(hot)
ylabel('glomerulus','FontSize',15)
set(gca,'ytick',1:length(gNames))
set(gca,'yticklabel',[string(gNames)])
set(gca,'xtick',1:nodors)
set(gca,'xticklabel',[string(odorNames)])
ytickangle(30)
xtickangle(30)
title('PNs')
set(gca,'FontSize',12)
hcb=colorbar;
title(hcb,'peak dF/F')

figure
imagesc(ormatrix,ctrst)
%colormap(cm.egoalley)
colormap(hot)
ylabel('glomerulus','FontSize',15)
set(gca,'ytick',1:length(gNamesO))
set(gca,'yticklabel',[string(gNamesO)])
set(gca,'xtick',1:nodors)
set(gca,'xticklabel',[string(odorNames)])
ytickangle(30)
xtickangle(30)
title('ORNs')
set(gca,'FontSize',12)
hcb=colorbar;
title(hcb,'peak dF/F')

% rescale heatmaps for vectorizing heat map
g2=ghmatrix;
g2=g2-min(g2(:));
g2=255*g2/max(g2(:));
g2=g2+1;
vectorPixels(g2,hot(256),[0 0 0])
o2=ormatrix;
o2=o2-min(o2(:));
o2=255*o2/max(o2(:));
o2=o2+1;
vectorPixels(o2,hot(256),[0 0 0])

% for making colorbar
figure;
imagesc(1:256)
colormap(hot)

figure
for j=1:length(gNames)
    subplot(1,5,j)
plot(orc(j*(1:nodors)),gh(j*(1:nodors)),'k.','LineWidth',3,'MarkerSize',15)
xlabel('ORN activation (dF/F)')
ylabel('PN activation (dF/F)')
set(gca,'FontSize',15)
title(gNames{j})
box off
end

% plot response correlation between glomerulus measured in ORNs and PNs
figure
imagesc(corr(ormatrix',ghmatrix'))
set(gca,'ytick',1:length(gNamesO))
set(gca,'yticklabel',[string(gNamesO)])
set(gca,'xtick',1:length(gNames))
set(gca,'xticklabel',[string(gNames)])
h=colorbar;
xlabel('PN')
ylabel('ORN')
title(h,'r')
set(gca,'FontSize',15)
% 
% figure;
% [oy ox]=hist(orc,[-.1:0.1:1.5]);
% [gy gx]=hist(gh,[-.1:0.1:1.5]);
% plot(ox,oy/sum(oy),'Color',orncolor,'LineWidth',3)
% hold on
% plot(gx,gy/sum(gy),'Color',pncolor,'LineWidth',3)
% xlabel('activation (dF/F)')
% ylabel('fraction of responses')
% legend('ORNs','PNs')
% legend boxoff
% set(gca,'FontSize',15)
% box off
% 
% 
% figure;
% [oy ox]=hist(orcoRawData(:),[-.5:0.05:4.5]);
% [gy gx]=hist(gh146rawdata(:),[-.5:0.05:4.5]);
% plot(ox,oy/sum(oy),'Color',orncolor,'LineWidth',3)
% hold on
% plot(gx,gy/sum(gy),'Color',pncolor,'LineWidth',3)
% xlabel('activation (dF/F)')
% ylabel('fraction of responses')
% legend('ORNs','PNs')
% legend boxoff
% set(gca,'FontSize',15)
% box off

% compare odor rank correlations between gh146 and orco with published ORN
% data
publishedORNresponsegh146=publishedOR.gh146response; % keep gloms expressed in gh146
publishedORNresponseorco=publishedOR.orcoresponse; % keep gloms expressed in orco

publishedORNidxgh146=zeros(1,length(gNames));
publishedORNidxorco=zeros(1,length(gNames));
for i=1:length(gNames)
    for j=1:length(publishedOR.gh146glomerulusNames)
        if strcmp(publishedOR.gh146glomerulusNames{j},gNames{i})
            publishedORNidxgh146(i)=j;
        end
    end
       for j=1:length(publishedOR.orcoGlomerulusNames)
        if strcmp(publishedOR.orcoGlomerulusNames{j},gNames{i})
            publishedORNidxorco(i)=j;
        end
    end
end

meanPublishedResponsegh146=NaN*zeros(1,length(gNames)*nodors);
meanPublishedResponseorco=NaN*zeros(1,length(gNames)*nodors);
% create mean published response for available glomeruli
for i=1:length(gNames)
    meanPublishedResponsegh146(((i-1)*nodors+2):(i)*nodors)=publishedORNresponsegh146(publishedORNidxgh146(i),:);
    meanPublishedResponseorco(((i-1)*nodors+2):(i)*nodors)=publishedORNresponseorco(publishedORNidxorco(i),:);
end
meanPublishedResponsegh146(1:nodors:end)=NaN; % remove air
meanPublishedResponseorco(1:nodors:end)=NaN;  % remove air
orcorr=corrcoef(meanPublishedResponseorco(isfinite(meanPublishedResponseorco)),orc(isfinite(meanPublishedResponseorco)));
ghcorr=corrcoef(meanPublishedResponsegh146(isfinite(meanPublishedResponsegh146)),gh(isfinite(meanPublishedResponsegh146)));

% generate shuffled correlations
nshuffles=1000;
for ii=1:nshuffles
    publishedORNidxSgh146=zeros(1,length(gNames));
    publishedORNidxSorco=zeros(1,length(gNames));
    shuffledgnamesgh146=randperm(length(publishedOR.gh146glomerulusNames));
    shuffledgnamesorco=randperm(length(publishedOR.orcoGlomerulusNames));
    for i=1:length(gNames)
        for j=1:length(publishedOR.gh146glomerulusNames)
            if strcmp(publishedOR.gh146glomerulusNames{shuffledgnamesgh146(j)},gNames{i})
                publishedORNidxSgh146(i)=j;
            end
        end
        for j=1:length(publishedOR.orcoGlomerulusNames)
            if strcmp(publishedOR.orcoGlomerulusNames{shuffledgnamesorco(j)},gNames{i})
                publishedORNidxSorco(i)=j;
            end
        end
    end
    
    meanPublishedResponseSgh146=NaN*zeros(1,length(gNames)*nodors);
    meanPublishedResponseSorco=NaN*zeros(1,length(gNames)*nodors);
    % create mean published response for available glomeruli
    for i=1:length(gNames)
        meanPublishedResponseSgh146(((i-1)*nodors+2):(i)*nodors)=publishedORNresponsegh146(publishedORNidxSgh146(i),:);
        meanPublishedResponseSorco(((i-1)*nodors+2):(i)*nodors)=publishedORNresponseorco(publishedORNidxSorco(i),:);
    end
    orcorrS=corrcoef(meanPublishedResponseSorco(isfinite(meanPublishedResponseSorco)),orc(isfinite(meanPublishedResponseSorco)));
    ghcorrS=corrcoef(meanPublishedResponseSgh146(isfinite(meanPublishedResponseSgh146)),gh(isfinite(meanPublishedResponseSgh146)));
    
    orcorrShuffled(ii)=orcorrS(1,2);
    ghcorrShuffled(ii)=ghcorrS(1,2);
    if mod(ii,10)==0
        disp(['iter ' num2str(ii)])
    end
end

figure;
plot(meanPublishedResponseorco,orc,'o','Color',orncolor,'LineWidth',3)
hold on
plot(meanPublishedResponsegh146,gh,'*','Color',pncolor,'LineWidth',3)
xlabel('DoOR ORN data')
ylabel('calcium response (dF/F)')
legend('ORNs','PNs')
legend boxoff
box off
text(0,1,['ORN r = ' num2str(orcorr(1,2))],'FontSize',15)
text(0,0,['PN r = ' num2str(ghcorr(1,2))],'FontSize',15)
set(gca,'FontSize',15)

odornames=[{'air'},{'3-octanol'},{'1-hexanol'},{'ethyl lactate'},{'citronella'},{'2-heptanone'},{'1-pentanol'},{'ethanol'},{'geranyl acetate'},{'hexyl acetate'},{'4-methylcyclohexanol'},{'pentyl acetate'},{'1-butanol'}];
% generate odor response distance maps
ghOdorResponseMapDistance=zeros(nodors);
orOdorResponseMapDistance=zeros(nodors);
for i=1:nodors
   for j=1:nodors
       if j>=i
           %ghOdorResponseMapDistance(i,j)=mean(((ghmatrix(:,i)-ghmatrix(:,j)).^2)./((ghmatrix(:,i)+ghmatrix(:,j)).^2));
           %orOdorResponseMapDistance(i,j)=mean(((ormatrix(:,i)-ormatrix(:,j)).^2)./((ormatrix(:,i)+ormatrix(:,j)).^2));
           
           %ghOdorResponseMapDistance(i,j)=corr(ghmatrix(:,i),ghmatrix(:,j));
           %orOdorResponseMapDistance(i,j)=corr(ormatrix(:,i),ormatrix(:,j));
           
           for k=1:size(gh146flyaverage,2)
               ghOdorResponseMapDistance(i,j)=ghOdorResponseMapDistance(i,j)+corr(gh146flyaverage(i:nodors:end,k),gh146flyaverage(j:nodors:end,k));
           end
           ghOdorResponseMapDistance(i,j)=ghOdorResponseMapDistance(i,j)/size(gh146flyaverage,2);
           for k=1:size(orcoflyaverage,2)
               orOdorResponseMapDistance(i,j)=orOdorResponseMapDistance(i,j)+corr(orcoflyaverage(i:nodors:end,k),orcoflyaverage(j:nodors:end,k));
           end
           orOdorResponseMapDistance(i,j)=orOdorResponseMapDistance(i,j)/size(orcoflyaverage,2);
       end
   end
end
figure;
imagesc(ghOdorResponseMapDistance,[-1 1])
set(gca,'XTick',[1:nodors])
set(gca,'XTickLabel',string(odornames))
xtickangle(30)
set(gca,'YTick',[1:nodors])
set(gca,'YTickLabel',string(odornames))
ytickangle(30)
title('PN')
set(gca,'FontSize',15)

figure;
imagesc(orOdorResponseMapDistance,[-1 1])
set(gca,'XTick',[1:nodors])
set(gca,'XTickLabel',string(odornames))
xtickangle(30)
set(gca,'YTick',[1:nodors])
set(gca,'YTickLabel',string(odornames))
ytickangle(30)
title('ORN')
set(gca,'FontSize',15)


% get glomerulus-glomerulus distances
ghGlomResponseMapDistance=zeros(length(gNames));
orGlomResponseMapDistance=zeros(length(gNamesO));

for i=1:length(gNames)
   for j=1:length(gNamesO)
       if j>=i
           
           for k=1:size(gh146flyaverage,2)
               ghGlomResponseMapDistance(i,j)=ghGlomResponseMapDistance(i,j)+corr(gh146flyaverage(((i-1)*nodors+1):((i)*nodors),k),gh146flyaverage(((j-1)*nodors+1):((j)*nodors),k));
           end
           ghGlomResponseMapDistance(i,j)=ghGlomResponseMapDistance(i,j)/size(gh146flyaverage,2);
           for k=1:size(orcoflyaverage,2)
               orGlomResponseMapDistance(i,j)=orGlomResponseMapDistance(i,j)+corr(orcoflyaverage(((i-1)*nodors+1):((i)*nodors),k),orcoflyaverage(((j-1)*nodors+1):((j)*nodors),k));
           end
           orGlomResponseMapDistance(i,j)=orGlomResponseMapDistance(i,j)/size(orcoflyaverage,2);

       end
   end
end
figure;
imagesc(ghGlomResponseMapDistance,[-1 1])
set(gca,'YTick',[1:length(gNames)])
set(gca,'YTickLabel',string(gNames))
xtickangle(30)
set(gca,'XTick',[1:length(gNames)])
set(gca,'XTickLabel',string(gNames))
ytickangle(30)
title('PN')
set(gca,'FontSize',15)

figure;
imagesc(orGlomResponseMapDistance,[-1 1])
set(gca,'XTick',[1:length(gNamesO)])
set(gca,'XTickLabel',string(gNames))
xtickangle(30)
set(gca,'YTick',[1:length(gNamesO)])
set(gca,'YTickLabel',string(gNamesO))
ytickangle(30)
title('ORN')
set(gca,'FontSize',15)


% plot distribution of glomerulus activations across individuals
glomactivation=zeros(length(gNames),nodors*flyNum);
glomactivationO=zeros(length(gNamesO),nodors*flyNumO);
for i=1:length(gNames)
    temp=gh146flyaverage(((i-1)*nodors+1):((i)*nodors),:);
   glomactivation(i,:)=temp(:);
end
for i=1:length(gNamesO)
    temp=orcoflyaverage(((i-1)*nodors+1):((i)*nodors),:);
   glomactivationO(i,:)=temp(:);
end

figure
subplot(1,2,1)
distributionPlot(glomactivation','histOpt',1,'colormap',1-gray(64),'showMM',0)
xlabel('PNs')
ylabel('df/f')
axis([0 length(gNames)+1 min(glomactivation(:)) max(glomactivation(:))])
box off
set(gca,'XTick',1:length(gNames))
set(gca,'XTickLabel',gNames)
xtickangle(30)
set(gca,'FontSize',15)
subplot(1,2,2)
distributionPlot(glomactivationO','histOpt',1,'colormap',1-gray(64),'showMM',0)
xlabel('ORNs')
ylabel('df/f')
axis([0 length(gNames)+1 min(glomactivationO(:)) max(glomactivationO(:))])
box off
set(gca,'XTick',1:length(gNamesO))
set(gca,'XTickLabel',gNamesO)
xtickangle(30)
set(gca,'FontSize',15)

% plot std vs mean of each glom df/f
figure
plot(mean(glomactivation,2),std(glomactivation'),'k*','LineWidth',3,'MarkerSize',10)
hold on
text(mean(glomactivation,2)+0.01,std(glomactivation'),gNames,'FontSize',15)
x=[0.1:0.1:0.6];
xlabel('pn glomerulus df/f \mu')
ylabel('pn glomerulus df/f \sigma')
plot(x,x,'--','Color',[0.65 0.65 0.65],'LineWidth',2)
box off
set(gca,'FontSize',15)

figure
plot(mean(glomactivationO,2),std(glomactivationO'),'k*','LineWidth',3,'MarkerSize',10)
hold on
text(mean(glomactivationO,2)+0.01,std(glomactivationO'),gNames,'FontSize',15)
x=[0.0:0.1:0.7];
xlabel('orn glomerulus df/f \mu')
ylabel('orn glomerulus df/f \sigma')
plot(x,x,'--','Color',[0.65 0.65 0.65],'LineWidth',2)
box off
set(gca,'FontSize',15)

%% plot odors as observations in pc space
gh146reshaped = zeros(13,size(responsesNoResponseRemoved,1)*size(responsesNoResponseRemoved,2)/13);
orcoreshaped = zeros(13,size(responsesNoResponseRemovedO,1)*size(responsesNoResponseRemovedO,2)/13);
for i = 1:13
    temp=responsesNoResponseRemoved(i:13:end,:);
    gh146reshaped(i,:)=temp(:);
    
    temp=responsesNoResponseRemovedO(i:13:end,:);
    orcoreshaped(i,:)=temp(:);
end

opt = statset('pca');
opt.Display='iter';
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(gh146reshaped,'Options',opt);
[COEFFO, SCOREO, LATENTO, TSQUAREDO, EXPLAINEDO] = pca(orcoreshaped,'Options',opt);

mycolors=hsv(13);
figure
scatter(SCORE(:,1),SCORE(:,2),70,mycolors,'filled')
text(SCORE(:,1),SCORE(:,2),odornames)
title('PNs')
xlabel('PC 1')
ylabel('PC 2')
box off
set(gca,'FontSize',15)

figure
scatter(SCOREO(:,1),SCOREO(:,2),70,mycolors,'filled')
text(SCOREO(:,1),SCOREO(:,2),odornames)
title('ORNs')
xlabel('PC 1')
ylabel('PC 2')
box off
set(gca,'FontSize',15)


%% plot responses of oct and mch

% ASSUMING 5 glomeruli = DC2, DL5, DM1, DM2, and DM3
 
% pns
dc2anddm2=[gh146flyaverage(41,:)' gh146flyaverage(50,:)' gh146flyaverage(40,:)' gh146flyaverage(2,:)' gh146flyaverage(11,:)' gh146flyaverage(1,:)' ];

figure
distributionPlot(dc2anddm2,'histOpt',1,'colormap',1-gray(64),'showMM',0)
set(gca,'XTick',1:6)
set(gca,'XTickLabel',[{'DM2 - oct'}, {'DM2 - mch'}, {'DM2 - air'}, {'DC2 - oct'}, {'DC2 - mch'}, {'DC2 - air'}])
xtickangle(30)
ylabel('df/f')
set(gca,'FontSize',15)

dc2anddm2o=[orcoflyaverage(41,:)' orcoflyaverage(50,:)' orcoflyaverage(40,:)' orcoflyaverage(2,:)' orcoflyaverage(11,:)' orcoflyaverage(1,:)' ];


% calculate difference between oct and mch for each glomerulus
octminusmch=zeros(length(gNames),length(gNames));
for i=1:length(gNames)
    for j=1:length(gNames)
        if i<j
            o1=mean(gh146flyaverage((i-1)*nodors+2,:));
            o2=mean(gh146flyaverage((j-1)*nodors+2,:));
            m1=mean(gh146flyaverage((i-1)*nodors+11,:));
            m2=mean(gh146flyaverage((j-1)*nodors+11,:));
            
            octminusmch(i,j)=((o1-o2)-(m1-m2));
        end
    end
end

figure
imagesc(abs(octminusmch))


% % orns
% figure
% distributionPlot(dc2anddm2o,'histOpt',1,'colormap',1-gray(64),'showMM',0)
% set(gca,'XTick',1:6)
% set(gca,'XTickLabel',[{'DM2 - oct'}, {'DM2 - mch'}, {'DM2 - air'}, {'DC2 - oct'}, {'DC2 - mch'}, {'DC2 - air'}])
% xtickangle(30)
% ylabel('df/f')
% set(gca,'FontSize',15)
%% if desired, keep only trials with no missing data

deleteTrialsWithMissingData=0; % 1: delete trials with missing data

if deleteTrialsWithMissingData
    ghToDelete=isnan(mean(gh146rawdata));
    orcoToDelete=isnan(mean(orcoRawData));
    responsesNoResponseRemoved(:,ghToDelete)=[];
    responsesNoResponseRemovedO(:,orcoToDelete)=[];
    glombyodorflies(ghToDelete)=[];
    glombyodorfliesO(orcoToDelete)=[];
end
%% 

% PNs

opt = statset('pca');
opt.Display='iter';
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(responsesNoResponseRemoved','Options',opt);


% ORNs
opt = statset('pca');
opt.Display='iter';
[COEFFO, SCOREO, LATENTO, TSQUAREDO, EXPLAINEDO] = pca(responsesNoResponseRemovedO','Options',opt);


figure;
plot(EXPLAINED,'-','Color',pncolor,'LineWidth',3)
hold on
plot(EXPLAINEDO,'o-','Color',orncolor,'LineWidth',3)
ylabel('variance explained (%)')
xlabel('PC #')
legend('PNs','ORNs')
legend boxoff
box off
set(gca,'FontSize',20)

figure;
semilogy(EXPLAINED,'-','Color',pncolor,'LineWidth',3)
hold on
semilogy(EXPLAINEDO,'o-','Color',orncolor,'LineWidth',3)
ylabel('variance explained (%)')
xlabel('PC #')
legend('PNs','ORNs')
legend boxoff
box off
set(gca,'FontSize',20)

% plot pc loadings
figure;
for i=1:8
    subplot(2,4,i)
    plot(COEFF(:,i),'*','Color',pncolor,'LineWidth',2,'MarkerSize',8)
    hold on
    plot(zeros(1,length(COEFF(:,i))),'k--','LineWidth',3)
    j=1;
    for ii=1:nodors:length(COEFF)
        plot((ii-0.5)*ones(1,5), linspace(min(COEFF(:,i)),max(COEFF(:,i)),5),'k--','LineWidth',2)
        %text(i+floor(nodors/3),min(PCContribution),num2str(glomsFound(j)),'FontSize',15)
        j=j+1;
    end
    
    set(gca,'xtick',(1:nodors:length(COEFF(:,i)))+floor(nodors/2),'xticklabel',string(gNames),'FontSize',10)
    xtickangle(30)
    ylabel(['PC ' num2str(i) ' loadings'])
    box off
    set(gca,'FontSize',15)
end

figure;
for i=1:8
    subplot(2,4,i)
    plot(COEFFO(:,i),'*','Color',orncolor,'LineWidth',2,'MarkerSize',8)
    hold on
    plot(zeros(1,length(COEFFO(:,i))),'k--','LineWidth',3)
    j=1;
    for ii=1:nodors:length(COEFFO)
        plot((ii-0.5)*ones(1,5), linspace(min(COEFFO(:,i)),max(COEFFO(:,i)),5),'k--','LineWidth',2)
        %text(i+floor(nodors/3),min(PCContribution),num2str(glomsFound(j)),'FontSize',15)
        j=j+1;
    end
    
    set(gca,'xtick',(1:nodors:length(COEFFO(:,i)))+floor(nodors/2),'xticklabel',string(gNamesO),'FontSize',10)
    xtickangle(30)
    ylabel(['PC ' num2str(i) ' loadings'])
    box off
    set(gca,'FontSize',15)
end

%% plot PN data

varianceToKeep=100; % percent of variance to keep


co = SCORE;
totalVarianceExplained=cumsum(EXPLAINED);
pcsWithinVariance=find(totalVarianceExplained<varianceToKeep);
pcstouse=pcsWithinVariance(end);
disp(['Using ' num2str(pcstouse) ' PCs'])


withinLLobe=NaN*zeros(1,flyNum);
withinRLobe=NaN*zeros(1,flyNum);
withinDifferentLobe=NaN*zeros(1,flyNum);
acrossLLobe=NaN*zeros(1,flyNum);
acrossRLobe=NaN*zeros(1,flyNum);
acrossAllLobe=NaN*zeros(1,flyNum);
% plot an odor in odor space for each fly and each lobe
allFlies=1:(flyNum);


for j=1:(flyNum)
    
    ltemps=find(glombyodorflies==(j));
    
    lcurr=ltemps;
    
    
    rtemps=find(glombyodorflies==(j+leftRightSplitVal));
    
    rcurr=rtemps;
    
    % calculate within-fly distances
    if length(lcurr)>1
        withinLLobe(j)=sqrt(sum((co(lcurr(1),1:pcstouse)-co(lcurr(2),1:pcstouse)).^2));
    end
    if length(rcurr)>1
        withinRLobe(j)=sqrt(sum((co(rcurr(1),1:pcstouse)-co(rcurr(2),1:pcstouse)).^2));
    end
    if length(lcurr)>0 && length(rcurr)>0
        withinDifferentLobe(j)=0;
        for kk=1:length(lcurr)
            for jj=1:length(rcurr)
                withinDifferentLobe(j)=withinDifferentLobe(j)+sqrt(sum((co(lcurr(kk),1:pcstouse)-co(rcurr(jj),1:pcstouse)).^2));
            end
        end
        withinDifferentLobe(j)=withinDifferentLobe(j)/(length(lcurr)*length(rcurr));
    end
    
    % calculate across-fly distances
    ltempsAcross=find(glombyodorflies~=(j));
    ltempsAcross=ltempsAcross(find(glombyodorflies(ltempsAcross)~=(j+leftRightSplitVal)));
    ltempsAcross=ltempsAcross(find(glombyodorflies(ltempsAcross)<leftRightSplitVal));
    lcurrAcross=ltempsAcross;
    
    rtempsAcross=find(glombyodorflies~=(j+leftRightSplitVal));
    rtempsAcross=rtempsAcross(find(glombyodorflies(rtempsAcross)~=(j)));
    rtempsAcross=rtempsAcross(find(glombyodorflies(rtempsAcross)>leftRightSplitVal));
    rcurrAcross=rtempsAcross;
    
    if length(lcurr)>0
        acrossLLobe(j)=0;
        for jj=1:length(lcurrAcross)
            for kk=1:length(lcurr)
                acrossLLobe(j)=acrossLLobe(j)+sqrt(sum((co(lcurr(kk),1:pcstouse)-(co(lcurrAcross(jj),1:pcstouse))).^2));
            end
        end
        acrossLLobe(j)=(acrossLLobe(j))/(length(lcurrAcross)*length(lcurr));
    end
    
    if length(rcurr)>0
        acrossRLobe(j)=0;
        for jj=1:length(rcurrAcross)
            for kk=1:length(rcurr)
                acrossRLobe(j)=acrossRLobe(j)+sqrt(sum((co(rcurr(kk),1:pcstouse)-(co(rcurrAcross(jj),1:pcstouse))).^2));
            end
        end
        acrossRLobe(j)=(acrossRLobe(j))/(length(rcurrAcross)*length(rcurr));
    end
    
    allAcross=[lcurrAcross rcurrAcross];
    allCurr=[lcurr rcurr];
    
    if length(lcurr)>0 && length(rcurr)>0
        acrossAllLobe(j)=0;
        for jj=1:length(allAcross)
            for kk=1:length(allCurr)
                acrossAllLobe(j)=acrossAllLobe(j)+sqrt(sum((co(allCurr(kk),1:pcstouse)-(co(allAcross(jj),1:pcstouse))).^2));
            end
        end
        acrossAllLobe(j)=acrossAllLobe(j)/(length(allAcross)*length(allCurr));
    end
end

msize=10;
lsize=3;
%
% % plot average across all odors for flies with both lobe data
figure
hold on

% plot(0,0,'Marker',m{1},'Color','k')
% plot(0,0,'Marker',m{2},'Color','k')
%legend('Left Lobe','Right Lobe')
for j=1:(flyNum)
    
    ltemps=find(glombyodorflies==(j));
    lcurr=ltemps(1:end);
    rtemps=find(glombyodorflies==(j+1000));
    rcurr=rtemps(1:end);
 
    if length(lcurr)>0  
        plot(co(lcurr(1:size(lcurr,2)),1),co(lcurr(1:size(lcurr,2)),2),'Color',mycmap(j,:),'Marker','*','LineWidth',lsize,'MarkerSize',msize);
        
    end
    if length(rcurr)>0
        plot(co(rcurr(1:size(rcurr,2)),1),co(rcurr(1:size(rcurr,2)),2),'Color',mycmap(j,:),'Marker','o','LineWidth',lsize,'MarkerSize',msize);
    end
    
    if length(lcurr)>0  && length(rcurr)>0  
        plot([co(lcurr(1:size(lcurr,2)),1)' co(rcurr(1:size(rcurr,2)),1)' co(lcurr(1),1)],[co(lcurr(1:size(lcurr,2)),2)' co(rcurr(1:size(rcurr,2)),2)' co(lcurr(1),2)],'Color',mycmap(j,:),'LineStyle','--','LineWidth',lsize,'MarkerSize',msize);
    end
    
end
xlabel('PC 1 Score')
ylabel('PC 2 Score')
set(gca,'FontSize',15)

figure
hold on
for j=1:(flyNum)
    
    ltemps=find(glombyodorflies==(j));
    lcurr=ltemps(1:end);
    rtemps=find(glombyodorflies==(j+1000));
    rcurr=rtemps(1:end);
    
    
    if length(lcurr)>0
        h3=plot(mean(co(lcurr(1:size(lcurr,2)),1)),mean(co(lcurr(1:size(lcurr,2)),2)),'Color',mycmap(j,:),'Marker','*','LineWidth',lsize,'MarkerSize',msize);
        %text(mean(co(lcurr(1:size(lcurr,2)),1)),mean(co(lcurr(1:size(lcurr,2)),2)),[num2str(j)],'FontSize',15)
    end
    if  length(rcurr)>0
        h4=plot(mean(co(rcurr(1:size(rcurr,2)),1)),mean(co(rcurr(1:size(rcurr,2)),2)),'Color',mycmap(j,:),'Marker','o','LineWidth',lsize,'MarkerSize',msize);
    end
    if  length(lcurr)>0 && length(rcurr)>0
        h5=plot([mean(co(lcurr(1:size(lcurr,2)),1))' mean(co(rcurr(1:size(rcurr,2)),1))'],[mean(co(lcurr(1:size(lcurr,2)),2))' mean(co(rcurr(1:size(rcurr,2)),2))'],'Color',mycmap(j,:),'LineStyle','--','LineWidth',lsize,'MarkerSize',1);
    end
    
end
h1=plot(-100,-100,'Marker','o','MarkerSize',msize,'LineStyle','none','LineWidth',lsize,'Color','k');
h2=plot(-100,-100,'Marker','*','MarkerSize',msize,'LineStyle','none','LineWidth',lsize,'Color','k');
legend([h1 h2],{'Left Lobe','Right Lobe'})
legend boxoff
box off
axis([min(co(:,1)) max(co(:,1)) min(co(:,2)) max(co(:,2))])
xlabel('PC 1 Score')
ylabel('PC 2 Score')
set(gca,'FontSize',15)



%Boxplots of distances
% summarize across flies

withinleft=withinLLobe;
withinright=withinRLobe;
withinacross=withinDifferentLobe;
acrossleft=acrossLLobe;
acrossright=acrossRLobe;
acrossall=acrossAllLobe;

figure
boxplot([withinleft(:) withinright(:) withinacross(:) acrossleft(:) acrossright(:) acrossall(:)])
ylabel('Distance in Coding Space')
xlabels{1}='Within Fly (Left Lobe)';
xlabels{2}='Within Fly (Right Lobe)';
xlabels{3}='Within Fly (Opposite Lobes)';

xlabels{4}='Across Fly (Left Lobe)';
xlabels{5}='Across Fly (Right Lobe)';
xlabels{6}='Across Fly (Both Lobes)';
set(gca,'xtick',1:6,'xticklabel',xlabels,'FontSize',10)
xtickangle(30)
box off
set(gca,'FontSize',15)

%% plot ORN data

co = SCOREO;
totalVarianceExplained=cumsum(EXPLAINEDO);
pcsWithinVariance=find(totalVarianceExplained<varianceToKeep);
pcstouse=pcsWithinVariance(end);
disp(['Using ' num2str(pcstouse) ' PCs'])


withinLLobe=NaN*zeros(1,flyNumO);
withinRLobe=NaN*zeros(1,flyNumO);
withinDifferentLobe=NaN*zeros(1,flyNumO);
acrossLLobe=NaN*zeros(1,flyNumO);
acrossRLobe=NaN*zeros(1,flyNumO);
acrossAllLobe=NaN*zeros(1,flyNumO);
% plot an odor in odor space for each fly and each lobe
allFlies=1:(flyNumO);


for j=1:(flyNumO)
    
    ltemps=find(glombyodorfliesO==(j));
    
    lcurr=ltemps;
    
    
    rtemps=find(glombyodorfliesO==(j+leftRightSplitVal));
    
    rcurr=rtemps;
    
    % calculate within-fly distances
    if length(lcurr)>1
        withinLLobe(j)=sqrt(sum((co(lcurr(1),1:pcstouse)-co(lcurr(2),1:pcstouse)).^2));
    end
    if length(rcurr)>1
        withinRLobe(j)=sqrt(sum((co(rcurr(1),1:pcstouse)-co(rcurr(2),1:pcstouse)).^2));
    end
    if length(lcurr)>0 && length(rcurr)>0
        withinDifferentLobe(j)=0;
        for kk=1:length(lcurr)
            for jj=1:length(rcurr)
                withinDifferentLobe(j)=withinDifferentLobe(j)+sqrt(sum((co(lcurr(kk),1:pcstouse)-co(rcurr(jj),1:pcstouse)).^2));
            end
        end
        withinDifferentLobe(j)=withinDifferentLobe(j)/(length(lcurr)*length(rcurr));
    end
    
    % calculate across-fly distances
    ltempsAcross=find(glombyodorfliesO~=(j));
    ltempsAcross=ltempsAcross(find(glombyodorfliesO(ltempsAcross)~=(j+leftRightSplitVal)));
    ltempsAcross=ltempsAcross(find(glombyodorfliesO(ltempsAcross)<leftRightSplitVal));
    lcurrAcross=ltempsAcross;
    
    rtempsAcross=find(glombyodorfliesO~=(j+leftRightSplitVal));
    rtempsAcross=rtempsAcross(find(glombyodorfliesO(rtempsAcross)~=(j)));
    rtempsAcross=rtempsAcross(find(glombyodorfliesO(rtempsAcross)>leftRightSplitVal));
    rcurrAcross=rtempsAcross;
    
    if length(lcurr)>0
        acrossLLobe(j)=0;
        for jj=1:length(lcurrAcross)
            for kk=1:length(lcurr)
                acrossLLobe(j)=acrossLLobe(j)+sqrt(sum((co(lcurr(kk),1:pcstouse)-(co(lcurrAcross(jj),1:pcstouse))).^2));
            end
        end
        acrossLLobe(j)=(acrossLLobe(j))/(length(lcurrAcross)*length(lcurr));
    end
    
    if length(rcurr)>0
        acrossRLobe(j)=0;
        for jj=1:length(rcurrAcross)
            for kk=1:length(rcurr)
                acrossRLobe(j)=acrossRLobe(j)+sqrt(sum((co(rcurr(kk),1:pcstouse)-(co(rcurrAcross(jj),1:pcstouse))).^2));
            end
        end
        acrossRLobe(j)=(acrossRLobe(j))/(length(rcurrAcross)*length(rcurr));
    end
    
    allAcross=[lcurrAcross rcurrAcross];
    allCurr=[lcurr rcurr];
    
    if length(lcurr)>0 && length(rcurr)>0
        acrossAllLobe(j)=0;
        for jj=1:length(allAcross)
            for kk=1:length(allCurr)
                acrossAllLobe(j)=acrossAllLobe(j)+sqrt(sum((co(allCurr(kk),1:pcstouse)-(co(allAcross(jj),1:pcstouse))).^2));
            end
        end
        acrossAllLobe(j)=acrossAllLobe(j)/(length(allAcross)*length(allCurr));
    end
end

msize=10;
lsize=3;
%
% % plot average across all odors for flies with both lobe data
figure
hold on

% plot(0,0,'Marker',m{1},'Color','k')
% plot(0,0,'Marker',m{2},'Color','k')
%legend('Left Lobe','Right Lobe')
for j=1:(flyNumO)
    
    ltemps=find(glombyodorfliesO==(j));
    lcurr=ltemps(1:end);
    rtemps=find(glombyodorfliesO==(j+1000));
    rcurr=rtemps(1:end);
 
    if length(lcurr)>0  
        plot(co(lcurr(1:size(lcurr,2)),1),co(lcurr(1:size(lcurr,2)),2),'Color',mycmapO(j,:),'Marker','*','LineWidth',lsize,'MarkerSize',msize);
        
    end
    if length(rcurr)>0
        plot(co(rcurr(1:size(rcurr,2)),1),co(rcurr(1:size(rcurr,2)),2),'Color',mycmapO(j,:),'Marker','o','LineWidth',lsize,'MarkerSize',msize);
    end
    
    if length(lcurr)>0  && length(rcurr)>0  
        plot([co(lcurr(1:size(lcurr,2)),1)' co(rcurr(1:size(rcurr,2)),1)' co(lcurr(1),1)],[co(lcurr(1:size(lcurr,2)),2)' co(rcurr(1:size(rcurr,2)),2)' co(lcurr(1),2)],'Color',mycmapO(j,:),'LineStyle','--','LineWidth',lsize,'MarkerSize',msize);
    end
    
end
xlabel('PC 1 Score')
ylabel('PC 2 Score')
set(gca,'FontSize',15)

figure
hold on
for j=1:(flyNumO)
    
    ltemps=find(glombyodorfliesO==(j));
    lcurr=ltemps(1:end);
    rtemps=find(glombyodorfliesO==(j+1000));
    rcurr=rtemps(1:end);
    
    
    if length(lcurr)>0
        h3=plot(mean(co(lcurr(1:size(lcurr,2)),1)),mean(co(lcurr(1:size(lcurr,2)),2)),'Color',mycmapO(j,:),'Marker','*','LineWidth',lsize,'MarkerSize',msize);
        %text(mean(co(lcurr(1:size(lcurr,2)),1)),mean(co(lcurr(1:size(lcurr,2)),2)),[num2str(j)],'FontSize',15)
    end
    if  length(rcurr)>0
        h4=plot(mean(co(rcurr(1:size(rcurr,2)),1)),mean(co(rcurr(1:size(rcurr,2)),2)),'Color',mycmapO(j,:),'Marker','o','LineWidth',lsize,'MarkerSize',msize);
    end
    if  length(lcurr)>0 && length(rcurr)>0
        h5=plot([mean(co(lcurr(1:size(lcurr,2)),1))' mean(co(rcurr(1:size(rcurr,2)),1))'],[mean(co(lcurr(1:size(lcurr,2)),2))' mean(co(rcurr(1:size(rcurr,2)),2))'],'Color',mycmapO(j,:),'LineStyle','--','LineWidth',lsize,'MarkerSize',1);
    end
    
end
h1=plot(-100,-100,'Marker','o','MarkerSize',msize,'LineStyle','none','LineWidth',lsize,'Color','k');
h2=plot(-100,-100,'Marker','*','MarkerSize',msize,'LineStyle','none','LineWidth',lsize,'Color','k');
legend([h1 h2],{'Left Lobe','Right Lobe'})
legend boxoff
box off
axis([min(co(:,1)) max(co(:,1)) min(co(:,2)) max(co(:,2))])
xlabel('PC 1 Score')
ylabel('PC 2 Score')
set(gca,'FontSize',15)


%Boxplots of distances
% summarize across flies

withinleftO=withinLLobe;
withinrightO=withinRLobe;
withinacrossO=withinDifferentLobe;
acrossleftO=acrossLLobe;
acrossrightO=acrossRLobe;
acrossallO=acrossAllLobe;

figure
boxplot([withinleftO(:) withinrightO(:) withinacrossO(:) acrossleftO(:) acrossrightO(:) acrossallO(:)])
ylabel('Distance in Coding Space')
xlabels{1}='Within Fly (Left Lobe)';
xlabels{2}='Within Fly (Right Lobe)';
xlabels{3}='Within Fly (Opposite Lobes)';

xlabels{4}='Across Fly (Left Lobe)';
xlabels{5}='Across Fly (Right Lobe)';
xlabels{6}='Across Fly (Both Lobes)';
set(gca,'xtick',1:6,'xticklabel',xlabels,'FontSize',10)
xtickangle(30)
box off
set(gca,'FontSize',15)

%%



PNtoboxplot=[withinleft(:) withinright(:) withinacross(:) acrossleft(:) acrossright(:) acrossall(:)];
ORNtoboxplot=[ withinleftO(:) withinrightO(:) withinacrossO(:) acrossleftO(:) acrossrightO(:) acrossallO(:)];

combined=[PNtoboxplot(:); ORNtoboxplot(:)];
groupings=[];
colors=[];
for i=1:6
    groupings=[groupings i*ones(1,flyNum)];
    colors(i,:)=[0 0.8 0];
end
for i=7:12
    groupings=[groupings i*ones(1,flyNumO)];
    colors(i,:)=[0.9 0 0.7];
end


figure
boxplot(combined,groupings,'PlotStyle','compact','Colors',colors)
ylabel('coding space distance')
xlabels{1}='within (left)';
xlabels{2}='within (right)';
xlabels{3}='within (opposite)';

xlabels{4}='across (left)';
xlabels{5}='across (right)';
xlabels{6}='across (both)';

xlabels{7}='within (left)';
xlabels{8}='within (right)';
xlabels{9}='within (opposite)';

xlabels{10}='across (left)';
xlabels{11}='across (right)';
xlabels{12}='across (both)';

set(gca,'xtick',1:12,'xticklabel',xlabels,'FontSize',10)
xtickangle(30)
box off
set(gca,'FontSize',15)
% 
% 
% cc=cell(1,12);
% for i=1:6
%     cc{i}=[0 0 0];
% end
% for i=7:12
%     cc{i}=[0.65 0.65 0.65];
% end
% figure
% distributionPlot(combined,'groups',groupings,'histOpt',1,'color',cc,'showMM',0);





PNtoboxplot=[withinleft(:) withinright(:) withinacross(:)];
ORNtoboxplot=[ withinleftO(:) withinrightO(:) withinacrossO(:)];

combined=[PNtoboxplot(:); ORNtoboxplot(:)];
groupings=[];
colors=[];
for i=1:3
    groupings=[groupings i*ones(1,flyNum)];
    colors(i,:)=[0 0.8 0];
end
for i=4:6
    groupings=[groupings i*ones(1,flyNumO)];
    colors(i,:)=[0.9 0.0 0.7];
end

figure
boxplot(combined,groupings,'PlotStyle','compact','Colors',colors)
ylabel('Distance in Coding Space')
xlabels{1}='Within Fly (Left Lobe)';
xlabels{2}='Within Fly (Right Lobe)';
xlabels{3}='Within Fly (Opposite Lobes)';

xlabels{4}='Within Fly (Left Lobe)';
xlabels{5}='Within Fly (Right Lobe)';
xlabels{6}='Within Fly (Opposite Lobes)';

set(gca,'xtick',1:6,'xticklabel',xlabels,'FontSize',10)
xtickangle(30)
box off
set(gca,'FontSize',15)


% 
% cc=cell(1,6);
% for i=1:3
%     cc{i}=[0 0 0];
% end
% for i=4:6
%     cc{i}=[0.65 0.65 0.65];
% end
% figure
% distributionPlot(combined,'groups',groupings,'histOpt',1,'color',cc,'showMM',0);



PNtoboxplot=[acrossleft(:) acrossright(:) acrossall(:)];
ORNtoboxplot=[acrossleftO(:) acrossrightO(:) acrossallO(:)];

combined=[PNtoboxplot(:); ORNtoboxplot(:)];
groupings=[];
colors=[];
for i=1:3
    groupings=[groupings i*ones(1,flyNum)];
    colors(i,:)=[0 0.8 0];
end
for i=4:6
    groupings=[groupings i*ones(1,flyNumO)];
    colors(i,:)=[0.9 0.0 0.7];
end

figure
boxplot(combined,groupings,'PlotStyle','compact','Colors',colors)
ylabel('Distance in Coding Space')
xlabels{1}='Across Fly (Left Lobe)';
xlabels{2}='Across Fly (Right Lobe)';
xlabels{3}='Across Fly (Both Lobes)';

xlabels{4}='Across Fly (Left Lobe)';
xlabels{5}='Across Fly (Right Lobe)';
xlabels{6}='Across Fly (Both Lobes)';


set(gca,'xtick',1:6,'xticklabel',xlabels,'FontSize',10)
xtickangle(30)
box off
set(gca,'FontSize',15)

% 
% cc=cell(1,6);
% for i=1:3
%     cc{i}=[0 0 0];
% end
% for i=4:6
%     cc{i}=[0.65 0.65 0.65];
% end
% figure
% distributionPlot(combined,'groups',groupings,'histOpt',1,'color',cc,'showMM',0);


% compare within versus across fly 
PNtoboxplot=[withinleft(:) withinright(:) withinacross(:) acrossleft(:) acrossright(:) acrossall(:)];
ORNtoboxplot=[ withinleftO(:) withinrightO(:) withinacrossO(:) acrossleftO(:) acrossrightO(:) acrossallO(:)];

combined=[PNtoboxplot(:); ORNtoboxplot(:)];
groupings=[];
groupings=[groupings 1*ones(1,3*flyNum)];
groupings=[groupings 3*ones(1,3*flyNum)];


groupings=[groupings 2*ones(1,3*flyNumO)];
groupings=[groupings 4*ones(1,3*flyNumO)];

colors=[0 0.8 0; 0.9 0 0.7; 0 0.8 0; 0.9 0 0.7];
figure
boxplot(combined,groupings,'PlotStyle','compact','Colors',colors)
%boxplot(combined,groupings,'BoxStyle','filled','Colors',colors)
ylabel('coding space distance')
xlabels{1}='within';
xlabels{2}='across';
set(gca,'xtick',[1.5 3.5],'xticklabel',xlabels,'FontSize',10)
xtickangle(30)
box off
set(gca,'FontSize',15)

[s PNwithinVacross] = ttest2(combined(groupings==1),combined(groupings==3))
[s ORNwithinVacross] = ttest2(combined(groupings==2),combined(groupings==4))
[s acrossFlyORNvsPN] = ttest2(combined(groupings==3),combined(groupings==4))
[s withinFlyORNvsPN] = ttest2(combined(groupings==1),combined(groupings==2))


figure
errorbar(1,nanmean(combined(groupings==1)),nanstd(combined(groupings==1))/sqrt(flyNum),'s','LineWidth',4,'Color',pncolor,'MarkerSize',10)
hold on
errorbar(2,nanmean(combined(groupings==2)),nanstd(combined(groupings==2))/sqrt(flyNumO),'s','LineWidth',4,'Color',orncolor,'MarkerSize',10)
errorbar(3,nanmean(combined(groupings==3)),nanstd(combined(groupings==3))/sqrt(flyNum),'s','LineWidth',4,'Color',pncolor,'MarkerSize',10)
errorbar(4,nanmean(combined(groupings==4)),nanstd(combined(groupings==4))/sqrt(flyNumO),'s','LineWidth',4,'Color',orncolor,'MarkerSize',10)
ylabel('coding space distance')
xlabels{1}='within';
xlabels{2}='across';
xlim([0 5])
set(gca,'xtick',[1.5 3.5],'xticklabel',xlabels,'FontSize',10)
xtickangle(30)
box off
set(gca,'FontSize',15)



% % compare within versus across fly , one point per fly
% PNtoboxplot=[nanmean([withinleft(:) withinright(:) withinacross(:)],2) nanmean([acrossleft(:) acrossright(:) acrossall(:)],2)];
% ORNtoboxplot=[nanmean([withinleftO(:) withinrightO(:) withinacrossO(:)],2) nanmean([acrossleftO(:) acrossrightO(:) acrossallO(:)],2)];
% 
% combined=[PNtoboxplot(:); ORNtoboxplot(:)];
% groupings=[];
% groupings=[groupings 1*ones(1,flyNum)];
% groupings=[groupings 3*ones(1,flyNum)];
% 
% 
% groupings=[groupings 2*ones(1,flyNumO)];
% groupings=[groupings 4*ones(1,flyNumO)];
% 
% colors=[0 0.8 0; 0.9 0 0.7; 0 0.8 0; 0.9 0 0.7];
% figure
% boxplot(combined,groupings,'PlotStyle','compact','Colors',colors)
% %boxplot(combined,groupings,'BoxStyle','filled','Colors',colors)
% ylabel('coding space distance')
% xlabels{1}='within';
% xlabels{2}='across';
% set(gca,'xtick',[1.5 3.5],'xticklabel',xlabels,'FontSize',10)
% xtickangle(30)
% box off
% set(gca,'FontSize',15)
% 
% figure
% errorbar(1,nanmean(combined(groupings==1)),nanstd(combined(groupings==1))/sqrt(flyNum),'s','LineWidth',4,'Color',pncolor,'MarkerSize',10)
% hold on
% errorbar(2,nanmean(combined(groupings==2)),nanstd(combined(groupings==2))/sqrt(flyNumO),'s','LineWidth',4,'Color',orncolor,'MarkerSize',10)
% errorbar(3,nanmean(combined(groupings==3)),nanstd(combined(groupings==3))/sqrt(flyNum),'s','LineWidth',4,'Color',pncolor,'MarkerSize',10)
% errorbar(4,nanmean(combined(groupings==4)),nanstd(combined(groupings==4))/sqrt(flyNumO),'s','LineWidth',4,'Color',orncolor,'MarkerSize',10)
% ylabel('coding space distance')
% xlabels{1}='within';
% xlabels{2}='across';
% xlim([0 5])
% set(gca,'xtick',[1.5 3.5],'xticklabel',xlabels,'FontSize',10)
% xtickangle(30)
% box off
% set(gca,'FontSize',15)


