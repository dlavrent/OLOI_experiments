% Compare PN and ORN data with single ORN gcamp responses
% load PNs
clear all
close all

load individualORNgcampresponses190423

allColorMaps

pncolor=[0 0.8 0];
orncolor=[0.9 0 0.7];
rng('default')

manualLabelHome='/Users/mattchurgin/Dropbox/flyimaging/analysis/ORNvsPN_analysis/pnsadded181218plus2019';

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

fracIn=0.4;% best results when fracIn is high, ~0.5, only using high confidence glomeruli

responsesNoResponseRemoved=responsesGlomByOdor;
responsesNoResponseRemovedT=responsesTimeCourseGlomByOdor;
responsesTensor=responsesTensorFullyOrthogonal;


% remove fully empty rows
%totallyempty=~any(isfinite(responsesNoResponseRemoved),2);
%responsesNoResponseRemoved(totallyempty,:)=[];

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

manualLabelHome='/Users/mattchurgin/Dropbox/flyimaging/analysis/ORNvsPN_analysis/ornflies';

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

fracIn=0.25; % best results when fracIn is high, ~0.5, only using high confidence glomeruli

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

meanOrco=nanmean(orcoRawData,2);
meanGh146=nanmean(gh146rawdata,2);

meanOrcoReshaped=reshape(meanOrco,nodors,length(gNames));
meanGh146Reshaped=reshape(meanGh146,nodors,length(gNames));

orco_orn_corr = corr(meanOrcoReshaped,meanIndORNresponse);
gh146_orn_corr = corr(meanGh146Reshaped,meanIndORNresponse);

figure;
imagesc(orco_orn_corr,[0 0.6])
set(gca,'ytick',1:length(gNames))
set(gca,'yticklabel',[string(gNames)])
set(gca,'xtick',1:length(indORN_glom))
set(gca,'xticklabel',[string(indORN_glom)])
set(gca,'FontSize',15)
hcb=colorbar;
title(hcb,'r')
xlabel('single OR Ca++ response')
ylabel('Orco Ca++ response')


%% Plot normalized and raw odor activations for DM2 and DC2 in ORCO and Single-Or drivers

% Plot raw odor activations
% plot orco-dm2 and or22a
figure
subplot(2,2,1)
plot(meanOrcoReshaped(:,4),'o','LineWidth',3,'MarkerSize',10)
hold on
plot(meanIndORNresponse(:,3),'x','LineWidth',3,'MarkerSize',10)
ylabel('peak dF/F')
xlabel('odor')
legend('orco DM2','Or22a')
legend boxoff
box off
set(gca,'FontSize',15)

% plot orco-dc2 and or13a
subplot(2,2,2)
plot(meanOrcoReshaped(:,1),'o','LineWidth',3,'MarkerSize',10)
hold on
plot(meanIndORNresponse(:,2),'x','LineWidth',3,'MarkerSize',10)
ylabel('peak dF/F')
xlabel('odor')
legend('orco DC2','Or13a')
legend boxoff
box off
set(gca,'FontSize',15)

% plot orco-dm2 and or13a
subplot(2,2,3)
plot(meanOrcoReshaped(:,4),'o','LineWidth',3,'MarkerSize',10)
hold on
plot(meanIndORNresponse(:,2),'x','LineWidth',3,'MarkerSize',10)
ylabel('peak dF/F')
xlabel('odor')
legend('orco DM2','Or13a')
legend boxoff
box off
set(gca,'FontSize',15)

% plot orco-dc2 and or22a
subplot(2,2,4)
plot(meanOrcoReshaped(:,1),'o','LineWidth',3,'MarkerSize',10)
hold on
plot(meanIndORNresponse(:,3),'x','LineWidth',3,'MarkerSize',10)
ylabel('peak dF/F')
xlabel('odor')
legend('orco DC2','Or22a')
legend boxoff
box off
set(gca,'FontSize',15)

% 
% % % Plot normalized odor activations
% % plot orco-dm2 and or22a
% figure
% subplot(2,2,1)
% plot(meanOrcoReshaped(:,4)/mean(meanOrcoReshaped(:,4)),'o','LineWidth',3,'MarkerSize',10)
% hold on
% plot(meanIndORNresponse(:,3)/mean(meanIndORNresponse(:,3)),'x','LineWidth',3,'MarkerSize',10)
% ylabel('normalized response')
% xlabel('odor')
% legend('orco DM2','Or22a')
% legend boxoff
% box off
% set(gca,'FontSize',15)
% 
% % plot orco-dc2 and or13a
% subplot(2,2,2)
% plot(meanOrcoReshaped(:,1)/mean(meanOrcoReshaped(:,1)),'o','LineWidth',3,'MarkerSize',10)
% hold on
% plot(meanIndORNresponse(:,2)/mean(meanIndORNresponse(:,2)),'x','LineWidth',3,'MarkerSize',10)
% ylabel('normalized response')
% xlabel('odor')
% legend('orco DC2','Or13a')
% legend boxoff
% box off
% set(gca,'FontSize',15)
% 
% % plot orco-dm2 and or13a
% subplot(2,2,3)
% plot(meanOrcoReshaped(:,4)/mean(meanOrcoReshaped(:,4)),'o','LineWidth',3,'MarkerSize',10)
% hold on
% plot(meanIndORNresponse(:,2)/mean(meanIndORNresponse(:,2)),'x','LineWidth',3,'MarkerSize',10)
% ylabel('normalized response')
% xlabel('odor')
% legend('orco DM2','Or13a')
% legend boxoff
% box off
% set(gca,'FontSize',15)
% 
% % plot orco-dc2 and or22a
% subplot(2,2,4)
% plot(meanOrcoReshaped(:,1)/mean(meanOrcoReshaped(:,1)),'o','LineWidth',3,'MarkerSize',10)
% hold on
% plot(meanIndORNresponse(:,3)/mean(meanIndORNresponse(:,3)),'x','LineWidth',3,'MarkerSize',10)
% ylabel('normalized response')
% xlabel('odor')
% legend('orco DC2','Or22a')
% legend boxoff
% box off
% set(gca,'FontSize',15)
%% bootstrap sample individual ORN flies for correlation
iters=1000;
orco_orn_corr_bs=zeros(size(orco_orn_corr));
for i=1:iters
    indORNtemp=zeros(nodors,length(indORN));
   for j=1:length(indORN)
      temp=randperm(size(indORN{j},2));
      indORNtemp(:,j)=indORN{j}(:,temp(1));
   end
   orco_orn_corr_bs = orco_orn_corr_bs+corr(meanOrcoReshaped,indORNtemp);
end
orco_orn_corr_bs=orco_orn_corr_bs/iters;

figure;
imagesc(orco_orn_corr_bs)
set(gca,'ytick',1:length(gNames))
set(gca,'yticklabel',[string(gNames)])
set(gca,'xtick',1:length(indORN_glom))
set(gca,'xticklabel',[string(indORN_glom)])
set(gca,'FontSize',15)
colorbar
xlabel('single OR Ca++ response')
ylabel('Orco Ca++ response')