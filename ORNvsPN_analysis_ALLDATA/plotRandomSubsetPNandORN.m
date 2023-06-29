% Characterize PN and ORN calcium responses

% load PNs
clear all
close all

load ORN_PN_COLORS
allColorMaps

rng('default')

load analysis_dir_path

manualLabelHome=fullfile(analysis_dir_path, 'ORNvsPN_analysis_ALLDATA/pn_alldata');

publishedOdorPath=fullfile(analysis_dir_path, 'utilities/odorPanel_12_DoORData.mat');
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
                        responsesTensorTemp(k,oo,:)=grnResponse(oo,j,odortimesfortimecourse)-nanmedian(grnResponse(oo,j,1:5),3);
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

% perform pca on PN responses

clear responsesNoResponseRemoved

fracIn=0.4;% best results when fracIn is high, ~0.5, only using high confidence glomeruli

responsesNoResponseRemoved=responsesGlomByOdor;
responsesNoResponseRemovedT=responsesTimeCourseGlomByOdor;
responsesTensor=responsesTensorFullyOrthogonal;

% blank odors (when solenoids failed for example)

odorstoremove=[2 6 9 10 12];
blankThresh=[0.3 100 100 0.3 0.3]; % blank odor if 75% prctile response is below this (for odors that sporadically worked)

% hard coded all flies before fly 23.  1:88 !!! Hard coded!!
for i=1:length(odorstoremove)
    flyodorstoblank=find(prctile(responsesNoResponseRemoved(odorstoremove(i):nodors:end,:),75)<=blankThresh(i));
    for j=1:(size(responsesNoResponseRemoved,1)/nodors)
        temp=find(isfinite(responsesNoResponseRemoved(((j-1)*nodors+odorstoremove(i)),[1:88])));
        tofill=intersect(flyodorstoblank,temp);
        if isfinite(nanmean(responsesNoResponseRemoved(((j-1)*nodors+odorstoremove(i)),89:end)))
            responsesNoResponseRemoved(((j-1)*nodors+odorstoremove(i)),tofill)=nanmean(responsesNoResponseRemoved(((j-1)*nodors+odorstoremove(i)),89:end));
        end
    end
end

% remove MCH for flies on 181108
for j=1:(size(responsesNoResponseRemoved,1)/nodors)
    tofill=find(isfinite(responsesNoResponseRemoved(((j-1)*nodors+11),89:100)))+88;
    if isfinite(nanmean(responsesNoResponseRemoved(((j-1)*nodors+11),[1:89 101:end])))
        responsesNoResponseRemoved(((j-1)*nodors+11),tofill)=nanmean(responsesNoResponseRemoved(((j-1)*nodors+11),[1:89 101:end]));
    end
end

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

% load ORNs

manualLabelHome=fullfile(analysis_dir_path, 'ORNvsPN_analysis_ALLDATA/orn_alldata');

publishedOdorPath=fullfile(analysis_dir_path, 'utilities/odorPanel_12_DoORData.mat');
load(publishedOdorPath);
manualLabelledFolders=dir(manualLabelHome);
manualLabelledFolders=manualLabelledFolders(3:end);

labels=cell(1,length(manualLabelledFolders));
glomsannotatedO=zeros(1,length(manualLabelledFolders));
totalPutativeGloms=zeros(1,length(manualLabelledFolders));

nodors=13;
odortimes=[6:9]; % hard-coded specific time interval for summarizing odor response
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
        

        behaviorOccO=[behaviorOccO occ-preocc];  
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
                        responsesTensorTemp(k,oo,:)=grnResponse(oo,j,odortimesfortimecourse)-nanmedian(grnResponse(oo,j,1:5),3);
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


%% plot heat map and response traces for a random subset of flies
nflies=50;
rng(2)
pnflies = randperm(flyNum); pnflies=pnflies(1:nflies);
ornflies = randperm(flyNumO); ornflies=ornflies(1:nflies);

% plot average response for each fly over all trials
gh146flyaverage=zeros(size(gh146rawdata,1),length(flyindices));
orcoflyaverage=zeros(size(orcoRawData,1),length(flyindicesO));
for i=1:length(flyindices)
    gh146flyaverage(:,i)=nanmean(gh146rawdata(:,flyindices{i}),2);
end
for i=1:length(flyindicesO)
    orcoflyaverage(:,i)=nanmean(orcoRawData(:,flyindicesO{i}),2);
end


plotindividualodorpanels=1;
if plotindividualodorpanels
    % resize raw data for plotting
    ghrawplot=NaN*zeros(length(flyindices)*5,13*4);
    orrawplot=NaN*zeros(length(flyindicesO)*5,13*4);
    ii=1;
    for i = 1:length(flyindices)
        jj=1;
        for j = 1:length(flyindices{i})
            for k=0:4
                ghrawplot(ii+k,jj:(jj+12))=transpose(gh146rawdata((13*k+1):(13*k+13),flyindices{i}(j)));
            end
            jj=jj+13;
        end
        for j=(j+1):4
            for k=0:4
                ghrawplot(ii+k,jj:(jj+12))=6.2;
            end
            jj=jj+13;
        end
        ii=ii+5;
    end
    ii=1;
    for i = 1:length(flyindicesO)
        jj=1;
        for j = 1:length(flyindicesO{i})
            for k=0:4
                orrawplot(ii+k,jj:(jj+12))=transpose(orcoRawData((13*k+1):(13*k+13),flyindicesO{i}(j)));
            end
            jj=jj+13;
        end
        for j=(j+1):4
            for k=0:4
                orrawplot(ii+k,jj:(jj+12))=6.2;
            end
            jj=jj+13;
        end
        ii=ii+5;
    end
    % plot raw data for all flies individual odor panels without mean filling
    
    g2=ghrawplot;
    o2=orrawplot;
    minc = min(min(g2(:)),min(o2(:)));
    g2=g2-minc;
    o2=o2-minc;
    g2(isnan(g2))=-1;
    o2(isnan(o2))=-1;
    maxc = max(max(g2(:)),max(o2(:)));
    g2=255*g2/maxc;
    o2=255*o2/maxc;
    g2=g2+1;
    o2=o2+1;
    vectorPixels(g2,hot(256),[0 0 0])
    vectorPixels(o2,hot(256),[0 0 0])
end

% plot raw data for all flies (average response across trials) without mean filling

g2=gh146flyaverage;
o2=orcoflyaverage;
minc = min(min(g2(:)),min(o2(:)));
g2=g2-minc;
o2=o2-minc;
g2(isnan(g2))=0;
o2(isnan(o2))=0;
maxc = max(max(g2(:)),max(o2(:)));
g2=255*g2/maxc;
o2=255*o2/maxc;
g2=g2+1;
o2=o2+1;
% SUP FIG PNs heatmap individual glom-odor responses
vectorPixels(g2,hot(256),[0 0 0])
% SUP FIG ORNs heatmap individual glom-odor responses
vectorPixels(o2,hot(256),[0 0 0])

% mean fill
for i=1:size(gh146flyaverage,1)
    temp=find(isnan(gh146flyaverage(i,:)));
    gh146flyaverage(i,temp)=nanmean(gh146flyaverage(i,:));
end
for i=1:size(orcoflyaverage,1)
    temp=find(isnan(orcoflyaverage(i,:)));
    orcoflyaverage(i,temp)=nanmean(orcoflyaverage(i,:));
end

g2=gh146flyaverage(:,pnflies);
o2=orcoflyaverage(:,ornflies);
minc = min(min(g2(:)),min(o2(:)));
g2=g2-minc;
o2=o2-minc;
maxc = max(max(g2(:)),max(o2(:)));
g2=255*g2/maxc;
o2=255*o2/maxc;
g2=g2+1;
o2=o2+1;
% FIG odor/glom vs individuals PNs
vectorPixels(g2,hot(256),[0 0 0])
% FIG odor/glom vs individuals ORNs
vectorPixels(o2,hot(256),[0 0 0])

%% plot correlation matrix
g2 = corr(gh146flyaverage');
o2 = corr(orcoflyaverage');

minc = min(min(g2(:)),min(o2(:)));
g2=g2-minc;
o2=o2-minc;
maxc = max(max(g2(:)),max(o2(:)));
g2=255*g2/maxc;
o2=255*o2/maxc;
g2=g2+1;
o2=o2+1;
% SUPFIG PNs correlation matrix
vectorPixels(g2,cm.egoalley,[0 0 0])
axis square
% SUPFIG ORNs correlation matrix
vectorPixels(o2,cm.egoalley,[0 0 0])
axis square

figure; %2
imagesc(1:256);colormap(cm.egoalley)
set(gca,'xtick','')
set(gca,'ytick','')

%% plot odor responses for oct, mch, and dc2 and dm2
pnLobeMean=squeeze(nanmean(responsesTensor(pnflies,:,:,:,:),2));
pnMean=squeeze(nanmean(pnLobeMean,1));
ornLobeMean=squeeze(nanmean(responsesTensorO(ornflies,:,:,:,:),2));
ornMean=squeeze(nanmean(ornLobeMean,1));

% FIG df/f vs time DM2/DC2 ORN/PN OCT/MCH
figure %3
subplot(2,2,3)
currglom=1;
currodor=11;
s1=shadedErrorBar(1.2*(odortimesfortimecourse),nanmean(squeeze(ornLobeMean(:,currglom,currodor,:)),1),nanstd(squeeze(ornLobeMean(:,currglom,currodor,:)))/sqrt(size(ornLobeMean,1)),'lineprops','-b','patchSaturation',0.33); 
set(s1.mainLine,'LineWidth',2)
s1.mainLine.Color=ocolor;
s1.patch.FaceColor=ocolor;
hold on
s2=shadedErrorBar(1.2*(odortimesfortimecourse),nanmean(squeeze(pnLobeMean(:,currglom,currodor,:)),1),nanstd(squeeze(pnLobeMean(:,currglom,currodor,:)))/sqrt(size(pnLobeMean,1)),'lineprops','-b','patchSaturation',0.33); 
set(s2.mainLine,'LineWidth',2)
s2.mainLine.Color=pcolor;
s2.patch.FaceColor=pcolor;
text(0,0,'DC2 MCH')
plot(0:25,zeros(1,length(0:25)),'k--','LineWidth',1)
legend('ORNs','PNs')
legend boxoff
xlim([0 23])
ylim([-0.05 .5])
set(gca,'xtick','')
set(gca,'ytick','')
box on
set(gca,'FontSize',15)


subplot(2,2,4)
currglom=4;
currodor=11;
s1=shadedErrorBar(1.2*(odortimesfortimecourse),nanmean(squeeze(ornLobeMean(:,currglom,currodor,:)),1),nanstd(squeeze(ornLobeMean(:,currglom,currodor,:)))/sqrt(size(ornLobeMean,1)),'lineprops','-b','patchSaturation',0.33); 
set(s1.mainLine,'LineWidth',2)
s1.mainLine.Color=ocolor;
s1.patch.FaceColor=ocolor;
hold on
s2=shadedErrorBar(1.2*(odortimesfortimecourse),nanmean(squeeze(pnLobeMean(:,currglom,currodor,:)),1),nanstd(squeeze(pnLobeMean(:,currglom,currodor,:)))/sqrt(size(pnLobeMean,1)),'lineprops','-b','patchSaturation',0.33); 
set(s2.mainLine,'LineWidth',2)
s2.mainLine.Color=pcolor;
s2.patch.FaceColor=pcolor;
text(0,0,'DM2 MCH')
plot(0:25,zeros(1,length(0:25)),'k--','LineWidth',1)
legend('ORNs','PNs')
legend boxoff
xlim([0 23])
ylim([-0.05 0.6])
set(gca,'xtick','')
set(gca,'ytick','')
box on
set(gca,'FontSize',15)




subplot(2,2,1)
currglom=1;
currodor=2;
s1=shadedErrorBar(1.2*(odortimesfortimecourse),nanmean(squeeze(ornLobeMean(:,currglom,currodor,:)),1),nanstd(squeeze(ornLobeMean(:,currglom,currodor,:)))/sqrt(size(ornLobeMean,1)),'lineprops','-b','patchSaturation',0.33); 
set(s1.mainLine,'LineWidth',2)
s1.mainLine.Color=ocolor;
s1.patch.FaceColor=ocolor;
hold on
s2=shadedErrorBar(1.2*(odortimesfortimecourse),nanmean(squeeze(pnLobeMean(:,currglom,currodor,:)),1),nanstd(squeeze(pnLobeMean(:,currglom,currodor,:)))/sqrt(size(pnLobeMean,1)),'lineprops','-b','patchSaturation',0.33); 
set(s2.mainLine,'LineWidth',2)
s2.mainLine.Color=pcolor;
s2.patch.FaceColor=pcolor;
text(0,0,'DC2 OCT')
plot(0:25,zeros(1,length(0:25)),'k--','LineWidth',1)
legend('ORNs','PNs')
legend boxoff
xlim([0 23])
ylim([-0.05 0.9])
set(gca,'xtick','')
set(gca,'ytick','')
box on
set(gca,'FontSize',15)



subplot(2,2,2)
currglom=4;
currodor=2;
s1=shadedErrorBar(1.2*(odortimesfortimecourse),nanmean(squeeze(ornLobeMean(:,currglom,currodor,:)),1),nanstd(squeeze(ornLobeMean(:,currglom,currodor,:)))/sqrt(size(ornLobeMean,1)),'lineprops','-b','patchSaturation',0.33); 
set(s1.mainLine,'LineWidth',2)
s1.mainLine.Color=ocolor;
s1.patch.FaceColor=ocolor;
hold on
s2=shadedErrorBar(1.2*(odortimesfortimecourse),nanmean(squeeze(pnLobeMean(:,currglom,currodor,:)),1),nanstd(squeeze(pnLobeMean(:,currglom,currodor,:)))/sqrt(size(pnLobeMean,1)),'lineprops','-b','patchSaturation',0.33); 
set(s2.mainLine,'LineWidth',2)
s2.mainLine.Color=pcolor;
s2.patch.FaceColor=pcolor;
text(0,0,'DM2 OCT')
plot(0:25,zeros(1,length(0:25)),'k--','LineWidth',1)
legend('ORNs','PNs')
legend boxoff
xlim([0 23])
ylim([-0.05 1.1])
set(gca,'xtick','')
set(gca,'ytick','')
box on
set(gca,'FontSize',15)


%% plot all time - dependent odor traces
pnLobeMean=squeeze(nanmean(responsesTensor(:,:,:,:,:),2));
pnMean=squeeze(nanmean(pnLobeMean,1));
ornLobeMean=squeeze(nanmean(responsesTensorO(:,:,:,:,:),2));
ornMean=squeeze(nanmean(ornLobeMean,1));
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

% SUP FIG glomerulus-odor time-dependent responses
figure %4
k=0;
for j= 1:13
    for i = 1:5
        
        k=k+1;
        h=subplot(13,5,k);
        gcf
        
        
        currglom=i;
        currodor=j;
        s1=shadedErrorBar(1.2*(odortimesfortimecourse),nanmean(squeeze(ornLobeMean(:,currglom,currodor,:)),1),nanstd(squeeze(ornLobeMean(:,currglom,currodor,:)))/sqrt(size(ornLobeMean,1)),'lineprops','-b','patchSaturation',0.33);
        set(s1.mainLine,'LineWidth',2)
        s1.mainLine.Color=ocolor;
        s1.patch.FaceColor=ocolor;
        hold on
        s2=shadedErrorBar(1.2*(odortimesfortimecourse),nanmean(squeeze(pnLobeMean(:,currglom,currodor,:)),1),nanstd(squeeze(pnLobeMean(:,currglom,currodor,:)))/sqrt(size(pnLobeMean,1)),'lineprops','-b','patchSaturation',0.33);
        set(s2.mainLine,'LineWidth',2)
        s2.mainLine.Color=pcolor;
        s2.patch.FaceColor=pcolor;
        plot(0:25,zeros(1,length(0:25)),'k--','LineWidth',1)
        ylim([-0.2 1.5])
        xlim([0 23])
        set(gca,'xtick','')
        set(gca,'ytick','')
        box on
        set(gca,'FontSize',5)
    end
end
% 
%% Plot summary of all fly odor responses

gh=nanmean(gh146rawdata,2);
orc=nanmean(orcoRawData,2);


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



% rescale heatmaps for vectorizing heat map
g2=ghmatrix;
o2=ormatrix;
minc = min(min(g2(:)),min(o2(:)));
g2=g2-minc;
o2=o2-minc;
maxc = max(max(g2(:)),max(o2(:)));
g2=255*g2/maxc;
o2=255*o2/maxc;
g2=g2+1;
o2=o2+1;
% FIG heatmap glom vs odor responses in PNs
vectorPixels(g2,hot(256),[0 0 0])
% FIG heatmap glom vs odor responses in ORNs
vectorPixels(o2,hot(256),[0 0 0])

%[coeffP scoreP latentP tsqP explainedP] = pca(g2);
%[coeffO scoreO latentO tsqO explainedO] = pca(o2);

% for making colorbar
figure; %5
imagesc(1:256)
colormap(hot)

%% plot PN vs ORN odor by odor responses

% plot GH146 vs Orco activation

% SUP FIG PN vs ORN peak calcium response
figure %6
%plot(orc,gh,'k.','LineWidth',3,'MarkerSize',15)
colorms=parula(5);
hold on
for i=1:5
   plot(orc(((i-1)*13+1):(i*13)),gh(((i-1)*13+1):(i*13)),'.','Color',colorms(i,:),'LineWidth',2,'MarkerSize',15)
   %text(0,0,gNames{i},'FontSize',25,'Color',colorms(i,:))
end
xlabel('ORN activation (dF/F)')
ylabel('PN activation (dF/F)')
set(gca,'FontSize',15)
axis([0 1.8 0 1.8])
axis square
legend('DC2','DL5','DM1','DM2','DM3')
legend boxoff
set(gca,'xtick','')
set(gca,'ytick','')
box on
%% plot PN vs ORN activation for each glomerulus
figure %7
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
corrmat=corr(ormatrix',ghmatrix');
% SUP FIG heatmap correlation matrix PN vs ORN
figure %8
imagesc(corrmat)
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
g2=corrmat;
minc = min(corrmat(:));
g2=g2-minc;
maxc = max(g2(:));
g2=255*g2/maxc;
g2=g2+1;
vectorPixels(g2,parula(256),[0 0 0])
axis square

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
% plot(ox,oy/sum(oy),'Color',ocolor,'LineWidth',3)
% hold on
% plot(gx,gy/sum(gy),'Color',pcolor,'LineWidth',3)
% xlabel('activation (dF/F)')
% ylabel('fraction of responses')
% legend('ORNs','PNs')
% legend boxoff
% set(gca,'FontSize',15)
% box off

%% compare odor rank correlations between gh146 and orco with published ORN
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
orcorrShuffled=zeros(1,nshuffles);
ghcorrShuffled=zeros(1,nshuffles);
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

% SUP FIG calcium response vs DoOR response
figure; %9
plot(meanPublishedResponseorco,orc,'.','Color',ocolor,'LineWidth',2,'MarkerSize',15)
hold on
plot(meanPublishedResponsegh146,gh,'.','Color',pcolor,'LineWidth',2,'MarkerSize',15)
xlabel('DoOR ORN data')
ylabel('calcium response (dF/F)')
legend('ORNs','PNs')
legend boxoff
box on
text(0,1,['ORN r = ' num2str(orcorr(1,2),'%2.2f')],'FontSize',15)
text(0,0,['PN r = ' num2str(ghcorr(1,2),'%2.2f')],'FontSize',15)
set(gca,'FontSize',15)
axis([-0.05 0.9 0 1.8])
axis square
set(gca,'xtick','')
set(gca,'ytick','')

% SUP FIG calcium response vs DoOR by glomerulus
figure; %10
for i = 1:length(gNames)
    subplot(1,length(gNames),i)
    mprotemp=transpose(meanPublishedResponseorco(((i-1)*nodors+1):(i)*nodors));
    otemp=orc(((i-1)*nodors+1):(i)*nodors);
    mprptemp=transpose(meanPublishedResponsegh146(((i-1)*nodors+1):(i)*nodors));
    gtemp=gh(((i-1)*nodors+1):(i)*nodors);
    plot(mprotemp,otemp,'.','Color',ocolor,'LineWidth',2,'MarkerSize',15)
    hold on
    plot(mprptemp,gtemp,'.','Color',pcolor,'LineWidth',2,'MarkerSize',15)
   
    box on
    r1 = corrcoef(mprotemp(isfinite(mprotemp)),otemp(isfinite(mprotemp)));
    r2 = corrcoef(mprptemp(isfinite(mprptemp)),gtemp(isfinite(mprptemp)));
    text(0.1,0.3,['ORN r = ' num2str(r1(1,2),'%02.2f')],'FontSize',15)
    text(0.1,0.3,['PN r = ' num2str(r2(1,2),'%02.2f')],'FontSize',15)
    set(gca,'FontSize',15)
    set(gca,'xtick','')
    set(gca,'ytick','')

    axis square
    if i==1
        axis([0 0.6 0 1.5])
        ylabel('calcium response (dF/F)')
        legend('ORNs','PNs')
        legend boxoff
    elseif i==2
        axis([0 0.35 0 0.85])
    elseif i==3
        axis([-.1 0.65 0 1])
        xlabel('DoOR ORN data')
    elseif i==4
        axis([-.1 0.5 0 1.7])
    elseif i==5
        axis([-.1 0.9 0 1.3])
    end
end

%% 
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

figure; %11
imagesc(ghOdorResponseMapDistance,[-1 1])
set(gca,'XTick',[1:nodors])
set(gca,'XTickLabel',string(odornames))
xtickangle(30)
set(gca,'YTick',[1:nodors])
set(gca,'YTickLabel',string(odornames))
ytickangle(30)
title('PN')
set(gca,'FontSize',15)

figure; %12
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
figure; %13
imagesc(ghGlomResponseMapDistance,[-1 1])
set(gca,'YTick',[1:length(gNames)])
set(gca,'YTickLabel',string(gNames))
xtickangle(30)
set(gca,'XTick',[1:length(gNames)])
set(gca,'XTickLabel',string(gNames))
ytickangle(30)
title('PN')
set(gca,'FontSize',15)

figure; %14
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

figure %15
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
figure %16
plot(mean(glomactivation,2),std(glomactivation'),'k*','LineWidth',3,'MarkerSize',10)
hold on
text(mean(glomactivation,2)+0.01,std(glomactivation'),gNames,'FontSize',15)
x=[0.1:0.1:0.6];
xlabel('pn glomerulus df/f \mu')
ylabel('pn glomerulus df/f \sigma')
plot(x,x,'--','Color',[0.65 0.65 0.65],'LineWidth',2)
box off
set(gca,'FontSize',15)

figure %17
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
figure %18
scatter(SCORE(:,1),SCORE(:,2),70,mycolors,'filled')
text(SCORE(:,1),SCORE(:,2),odornames)
title('PNs')
xlabel('PC 1')
ylabel('PC 2')
box off
set(gca,'FontSize',15)

figure %19
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

figure %20
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

figure %21
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

% FIG variance explaned vs. PC
figure; %22
plot(cumsum(EXPLAINEDO),'-','Color',ocolor,'LineWidth',3,'MarkerSize',15)
hold on
plot(cumsum(EXPLAINED),'--','Color',pcolor,'LineWidth',3,'MarkerSize',15)
ylabel('cumulative variance explained (%)')
xlabel('PC #')
legend('ORN','PN')
legend boxoff
box on
axis([0 66 0 100])
set(gca,'xtick','')
set(gca,'ytick','')
set(gca,'FontSize',15)

% SUP FIG PN PC loadings grouped by glomerulus
figure; %23
for i=1:10
    subplot(2,5,i)
    plot(COEFF(:,i),'.','Color',pcolor,'LineWidth',2,'MarkerSize',15)
    hold on
    plot(zeros(1,length(COEFF(:,i))),'k--','LineWidth',3)
    j=1;
    for ii=(nodors+1):nodors:length(COEFF)
        plot((ii-0.5)*ones(1,5), linspace(-.55,.55,5),'k--','LineWidth',1)
        %text(i+floor(nodors/3),min(PCContribution),num2str(glomsFound(j)),'FontSize',15)
        j=j+1;
    end

    axis([0 size(COEFFO,2)+1 -.55 .55])
    set(gca,'xtick',(1:nodors:length(COEFF(:,i)))+floor(nodors/2),'xticklabel',string(gNames),'FontSize',10)
    xtickangle(30)
    title(['PC ' num2str(i) ' (' num2str(EXPLAINED(i),'%2.1f') '%)'])
    box on
    set(gca,'FontSize',15)
    set(gca,'xtick','')
    set(gca,'ytick','')
end

% SUP FIG ORN PC loadings grouped by glomerulus
figure; %24
for i=1:10
    subplot(2,5,i)
    plot(COEFFO(:,i),'.','Color',ocolor,'LineWidth',2,'MarkerSize',15)
    hold on
    plot(zeros(1,length(COEFFO(:,i))),'k--','LineWidth',3)
    j=1;
    for ii=(nodors+1):nodors:length(COEFFO)
        plot((ii-0.5)*ones(1,5), linspace(-.55,.55,5),'k--','LineWidth',1)
        %text(i+floor(nodors/3),min(PCContribution),num2str(glomsFound(j)),'FontSize',15)
        j=j+1;
    end
    axis([0 size(COEFFO,2)+1 -.55 .55])
    set(gca,'xtick',(1:nodors:length(COEFFO(:,i)))+floor(nodors/2),'xticklabel',string(gNamesO),'FontSize',10)
    xtickangle(30)
    title(['PC ' num2str(i) ' (' num2str(EXPLAINEDO(i),'%2.1f') '%)'])
    box on
    set(gca,'FontSize',15)
    set(gca,'xtick','')
    set(gca,'ytick','')
end

% SUP FIG PN PC loadings grouped by odor
figure; %25
for i=1:10
    subplot(2,5,i)
    currpc = COEFF(:,i);
    currpcr = zeros(size(currpc,1),1);
    ngloms=length(gNames);
    for j=1:nodors
        currpcr(((j-1)*ngloms+1):(j*ngloms))=currpc(j:nodors:end);
    end
    plot(currpcr,'.','Color',pcolor,'LineWidth',2,'MarkerSize',15)
    hold on
    plot(zeros(1,length(COEFF(:,i))),'k--','LineWidth',3)
    j=1;
    for ii=(ngloms+1):ngloms:length(currpcr)
        plot((ii-0.5)*ones(1,5), linspace(-.55,.55,5),'k--','LineWidth',1)
        %text(i+floor(nodors/3),min(PCContribution),num2str(glomsFound(j)),'FontSize',15)
        j=j+1;
    end

    axis([0 size(COEFFO,2)+1 -.55 .55])
    %set(gca,'xtick',(1:nodors:length(COEFF(:,i)))+floor(nodors/2),'xticklabel',string(gNames),'FontSize',10)
    xtickangle(30)
    title(['PC ' num2str(i) ' (' num2str(EXPLAINED(i),'%2.1f') '%)'])
    box on
    set(gca,'FontSize',15)
    set(gca,'xtick','')
    set(gca,'ytick','')
end


% SUP FIG ORN PC loadigns grouped by odor
figure; %26
for i=1:10
    subplot(2,5,i)
    currpc = COEFFO(:,i);
    currpcr = zeros(size(currpc,1),1);
    ngloms=length(gNamesO);
    for j=1:nodors
        currpcr(((j-1)*ngloms+1):(j*ngloms))=currpc(j:nodors:end);
    end
    plot(currpcr,'.','Color',ocolor,'LineWidth',2,'MarkerSize',15)
    hold on
    plot(zeros(1,length(COEFFO(:,i))),'k--','LineWidth',3)
    j=1;
    for ii=(ngloms+1):ngloms:length(currpcr)
        plot((ii-0.5)*ones(1,5), linspace(-.55,.55,5),'k--','LineWidth',1)
        %text(i+floor(nodors/3),min(PCContribution),num2str(glomsFound(j)),'FontSize',15)
        j=j+1;
    end

    axis([0 size(COEFFO,2)+1 -.55 .55])
    %set(gca,'xtick',(1:nodors:length(COEFF(:,i)))+floor(nodors/2),'xticklabel',string(gNames),'FontSize',10)
    xtickangle(30)
    title(['PC ' num2str(i) ' (' num2str(EXPLAINEDO(i),'%2.1f') '%)'])
    box on
    set(gca,'FontSize',15)
    set(gca,'xtick','')
    set(gca,'ytick','')
end


%% plot 25 pc loadings

% plot pc loadings  grouped by glomerulus
figure; %27
for i=1:25
    subplot(5,5,i)
    plot(COEFF(:,i),'.','Color',pcolor,'LineWidth',2,'MarkerSize',15)
    hold on
    plot(zeros(1,length(COEFF(:,i))),'k--','LineWidth',3)
    j=1;
    for ii=(nodors+1):nodors:length(COEFF)
        plot((ii-0.5)*ones(1,5), linspace(-.55,.55,5),'k--','LineWidth',1)
        %text(i+floor(nodors/3),min(PCContribution),num2str(glomsFound(j)),'FontSize',15)
        j=j+1;
    end

    axis([0 size(COEFFO,2)+1 -.55 .55])
    set(gca,'xtick',(1:nodors:length(COEFF(:,i)))+floor(nodors/2),'xticklabel',string(gNames),'FontSize',10)
    xtickangle(30)
    title(['PC ' num2str(i) ' (' num2str(EXPLAINED(i),'%2.1f') '%)'])
    box on
    set(gca,'FontSize',15)
    set(gca,'xtick','')
    set(gca,'ytick','')
end

figure; %28
for i=1:25
    subplot(5,5,i)
    plot(COEFFO(:,i),'.','Color',ocolor,'LineWidth',2,'MarkerSize',15)
    hold on
    plot(zeros(1,length(COEFFO(:,i))),'k--','LineWidth',3)
    j=1;
    for ii=(nodors+1):nodors:length(COEFFO)
        plot((ii-0.5)*ones(1,5), linspace(-.55,.55,5),'k--','LineWidth',1)
        %text(i+floor(nodors/3),min(PCContribution),num2str(glomsFound(j)),'FontSize',15)
        j=j+1;
    end
    axis([0 size(COEFFO,2)+1 -.55 .55])
    set(gca,'xtick',(1:nodors:length(COEFFO(:,i)))+floor(nodors/2),'xticklabel',string(gNamesO),'FontSize',10)
    xtickangle(30)
    title(['PC ' num2str(i) ' (' num2str(EXPLAINEDO(i),'%2.1f') '%)'])
    box on
    set(gca,'FontSize',15)
    set(gca,'xtick','')
    set(gca,'ytick','')
end

% plot pc loadings grouped by odor
figure; %29
for i=1:25
    subplot(5,5,i)
    currpc = COEFF(:,i);
    currpcr = zeros(size(currpc,1),1);
    ngloms=length(gNames);
    for j=1:nodors
        currpcr(((j-1)*ngloms+1):(j*ngloms))=currpc(j:nodors:end);
    end
    plot(currpcr,'.','Color',pcolor,'LineWidth',2,'MarkerSize',15)
    hold on
    plot(zeros(1,length(COEFF(:,i))),'k--','LineWidth',3)
    j=1;
    for ii=(ngloms+1):ngloms:length(currpcr)
        plot((ii-0.5)*ones(1,5), linspace(-.55,.55,5),'k--','LineWidth',1)
        %text(i+floor(nodors/3),min(PCContribution),num2str(glomsFound(j)),'FontSize',15)
        j=j+1;
    end

    axis([0 size(COEFFO,2)+1 -.55 .55])
    %set(gca,'xtick',(1:nodors:length(COEFF(:,i)))+floor(nodors/2),'xticklabel',string(gNames),'FontSize',10)
    xtickangle(30)
    title(['PC ' num2str(i) ' (' num2str(EXPLAINED(i),'%2.1f') '%)'])
    box on
    set(gca,'FontSize',15)
    set(gca,'xtick','')
    set(gca,'ytick','')
end


% plot pc loadings grouped by odor
figure; %30
for i=1:25
    subplot(5,5,i)
    currpc = COEFFO(:,i);
    currpcr = zeros(size(currpc,1),1);
    ngloms=length(gNamesO);
    for j=1:nodors
        currpcr(((j-1)*ngloms+1):(j*ngloms))=currpc(j:nodors:end);
    end
    plot(currpcr,'.','Color',ocolor,'LineWidth',2,'MarkerSize',15)
    hold on
    plot(zeros(1,length(COEFFO(:,i))),'k--','LineWidth',3)
    j=1;
    for ii=(ngloms+1):ngloms:length(currpcr)
        plot((ii-0.5)*ones(1,5), linspace(-.55,.55,5),'k--','LineWidth',1)
        %text(i+floor(nodors/3),min(PCContribution),num2str(glomsFound(j)),'FontSize',15)
        j=j+1;
    end

    axis([0 size(COEFFO,2)+1 -.55 .55])
    %set(gca,'xtick',(1:nodors:length(COEFF(:,i)))+floor(nodors/2),'xticklabel',string(gNames),'FontSize',10)
    xtickangle(30)
    title(['PC ' num2str(i) ' (' num2str(EXPLAINEDO(i),'%2.1f') '%)'])
    box on
    set(gca,'FontSize',15)
    set(gca,'xtick','')
    set(gca,'ytick','')
end
%% plot PN and ORN data in PC space

mycmap=hsv(nflies);

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


all_corrL = [ ] ;
all_corrR = [ ] ;
for j=1:(flyNum)
    
    ltemps=find(glombyodorflies==(j));
    
    lcurr=ltemps;
    
    
    rtemps=find(glombyodorflies==(j+leftRightSplitVal));
    
    rcurr=rtemps;
    
    % calculate within-fly distances
    if length(lcurr)>1
        withinLLobe(j)=sqrt(sum((co(lcurr(1),1:pcstouse)-co(lcurr(2),1:pcstouse)).^2));
        corrLtemp = corrcoef(co(lcurr(1),1:pcstouse), co(lcurr(2),1:pcstouse));
        all_corrL = [all_corrL corrLtemp(1,2)];
    end
    if length(rcurr)>1
        withinRLobe(j)=sqrt(sum((co(rcurr(1),1:pcstouse)-co(rcurr(2),1:pcstouse)).^2));
        corrRtemp = corrcoef(co(rcurr(1),1:pcstouse), co(rcurr(2),1:pcstouse));
        all_corrR = [all_corrR corrRtemp(1,2)];
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

figure %31
hold on
i=0;
%for j=1:(flyNum)
for j=pnflies  
    i=i+1;
    ltemps=find(glombyodorflies==(j));
    lcurr=ltemps(1:end);
    rtemps=find(glombyodorflies==(j+1000));
    rcurr=rtemps(1:end);
    
    
    if length(lcurr)>0
        h3=plot(mean(co(lcurr(1:size(lcurr,2)),1)),mean(co(lcurr(1:size(lcurr,2)),2)),'Color',mycmap(i,:),'Marker','*','LineWidth',lsize,'MarkerSize',msize);
        %text(mean(co(lcurr(1:size(lcurr,2)),1)),mean(co(lcurr(1:size(lcurr,2)),2)),[num2str(j)],'FontSize',15)
    end
    if  length(rcurr)>0
        h4=plot(mean(co(rcurr(1:size(rcurr,2)),1)),mean(co(rcurr(1:size(rcurr,2)),2)),'Color',mycmap(i,:),'Marker','o','LineWidth',lsize,'MarkerSize',msize);
    end
    if  length(lcurr)>0 && length(rcurr)>0
        h5=plot([mean(co(lcurr(1:size(lcurr,2)),1))' mean(co(rcurr(1:size(rcurr,2)),1))'],[mean(co(lcurr(1:size(lcurr,2)),2))' mean(co(rcurr(1:size(rcurr,2)),2))'],'Color',mycmap(i,:),'LineStyle','--','LineWidth',lsize,'MarkerSize',1);
    end
    
end
h1=plot(-100,-100,'Marker','o','MarkerSize',msize,'LineStyle','none','LineWidth',lsize,'Color','k');
h2=plot(-100,-100,'Marker','*','MarkerSize',msize,'LineStyle','none','LineWidth',lsize,'Color','k');
legend([h1 h2],{'Left Lobe','Right Lobe'})
legend boxoff
box on
set(gca,'xtick','')
set(gca,'ytick','')
axis([min(co(:,1)) max(co(:,1)) min(co(:,2)) max(co(:,2))])
xlabel(['PC 1 Score (' num2str(EXPLAINED(1),'%0.1f') '%)'])
ylabel(['PC 2 Score (' num2str(EXPLAINED(2),'%0.1f') '%)'])
set(gca,'FontSize',15)


withinleft=withinLLobe;
withinright=withinRLobe;
withinacross=withinDifferentLobe;
acrossleft=acrossLLobe;
acrossright=acrossRLobe;
acrossall=acrossAllLobe;


% plot ORN data

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
mycmapO=mycmap;

figure %32
hold on
i=0;
%for j=1:(flyNumO)
for j=ornflies 
    i=i+1;
    ltemps=find(glombyodorfliesO==(j));
    lcurr=ltemps(1:end);
    rtemps=find(glombyodorfliesO==(j+1000));
    rcurr=rtemps(1:end);
    
    
    if length(lcurr)>0
        h3=plot(mean(co(lcurr(1:size(lcurr,2)),1)),mean(co(lcurr(1:size(lcurr,2)),2)),'Color',mycmapO(i,:),'Marker','*','LineWidth',lsize,'MarkerSize',msize);
        %text(mean(co(lcurr(1:size(lcurr,2)),1)),mean(co(lcurr(1:size(lcurr,2)),2)),[num2str(j)],'FontSize',15)
    end
    if  length(rcurr)>0
        h4=plot(mean(co(rcurr(1:size(rcurr,2)),1)),mean(co(rcurr(1:size(rcurr,2)),2)),'Color',mycmapO(i,:),'Marker','o','LineWidth',lsize,'MarkerSize',msize);
    end
    if  length(lcurr)>0 && length(rcurr)>0
        h5=plot([mean(co(lcurr(1:size(lcurr,2)),1))' mean(co(rcurr(1:size(rcurr,2)),1))'],[mean(co(lcurr(1:size(lcurr,2)),2))' mean(co(rcurr(1:size(rcurr,2)),2))'],'Color',mycmapO(i,:),'LineStyle','--','LineWidth',lsize,'MarkerSize',1);
    end
    
end
h1=plot(-100,-100,'Marker','o','MarkerSize',msize,'LineStyle','none','LineWidth',lsize,'Color','k');
h2=plot(-100,-100,'Marker','*','MarkerSize',msize,'LineStyle','none','LineWidth',lsize,'Color','k');
legend([h1 h2],{'Left Lobe','Right Lobe'})
legend boxoff
box on
set(gca,'xtick','')
set(gca,'ytick','')
axis([min(co(:,1)) max(co(:,1)) min(co(:,2)) max(co(:,2))])
xlabel(['PC 1 Score (' num2str(EXPLAINEDO(1),'%0.1f') '%)'])
ylabel(['PC 2 Score (' num2str(EXPLAINEDO(2),'%0.1f') '%)'])
set(gca,'FontSize',15)

withinleftO=withinLLobe;
withinrightO=withinRLobe;
withinacrossO=withinDifferentLobe;
acrossleftO=acrossLLobe;
acrossrightO=acrossRLobe;
acrossallO=acrossAllLobe;

%%
%% plot PN and ORN data in PC space (random 2 trials)
nfliestoshow=20;

ginds = glombyodorflies; ginds(ginds>1000)=ginds(ginds>1000)-1000;
oinds = glombyodorfliesO; oinds(oinds>1000)=oinds(oinds>1000)-1000;
mycmap=hsv(nfliestoshow);

varianceToKeep=100; % percent of variance to keep

rng(42)
co = SCORE;
totalVarianceExplained=cumsum(EXPLAINED);
pcsWithinVariance=find(totalVarianceExplained<varianceToKeep);
pcstouse=pcsWithinVariance(end);
disp(['Using ' num2str(pcstouse) ' PCs'])

% plot an odor in odor space for each fly and each lobe
allFlies=1:(flyNum);

msize=15;
lsize=1;

% FIG variance explained vs. PC (PNs)
figure %33
hold on
i=0;
%for j=1:(flyNum)
for j=pnflies(1:nfliestoshow)
    i=i+1;
    ltemps=find(ginds==(j));
    lcurr=ltemps(1:end);
    
    if length(lcurr)>1
       temp=randperm(length(lcurr)); 
       h3=plot(co(lcurr(temp(1)),1),co(lcurr(temp(1)),2),'Color',mycmap(i,:),'Marker','*','LineWidth',lsize,'MarkerSize',msize);
       h4=plot(co(lcurr(temp(2)),1),co(lcurr(temp(2)),2),'Color',mycmap(i,:),'Marker','o','LineWidth',lsize,'MarkerSize',msize);
       h5=plot([co(lcurr(temp(1)),1) co(lcurr(temp(2)),1)],[co(lcurr(temp(1)),2) co(lcurr(temp(2)),2)],'Color',mycmap(i,:),'LineStyle','--','LineWidth',lsize,'MarkerSize',1);
    end
end
h1=plot(-100,-100,'Marker','o','MarkerSize',msize,'LineStyle','none','LineWidth',lsize,'Color','k');
h2=plot(-100,-100,'Marker','*','MarkerSize',msize,'LineStyle','none','LineWidth',lsize,'Color','k');
legend([h1 h2],{'trial 1','trial 2'})
legend boxoff
box on
set(gca,'xtick','')
set(gca,'ytick','')
axis square
axis([1.1*min(co(:,1)) 1.1*max(co(:,1)) 1.1*min(co(:,2)) 1.1*max(co(:,2))])
xlabel(['PC 1 Score (' num2str(EXPLAINED(1),'%0.1f') '%)'])
ylabel(['PC 2 Score (' num2str(EXPLAINED(2),'%0.1f') '%)'])
set(gca,'FontSize',15)

% plot ORN data

co = SCOREO;
totalVarianceExplained=cumsum(EXPLAINEDO);
pcsWithinVariance=find(totalVarianceExplained<varianceToKeep);
pcstouse=pcsWithinVariance(end);
disp(['Using ' num2str(pcstouse) ' PCs'])


% plot an odor in odor space for each fly and each lobe
allFlies=1:(flyNumO);


mycmapO=mycmap;

% FIG variance explained vs. PC (ORNs)
figure %34
hold on
i=0;
%for j=1:(flyNumO)
for j=ornflies(1:nfliestoshow)
    i=i+1;
    ltemps=find(oinds==(j));
    lcurr=ltemps(1:end);
    
    if length(lcurr)>1
       temp=randperm(length(lcurr)); 
       h3=plot(co(lcurr(temp(1)),1),co(lcurr(temp(1)),2),'Color',mycmap(i,:),'Marker','*','LineWidth',lsize,'MarkerSize',msize);
       h4=plot(co(lcurr(temp(2)),1),co(lcurr(temp(2)),2),'Color',mycmap(i,:),'Marker','o','LineWidth',lsize,'MarkerSize',msize);
       h5=plot([co(lcurr(temp(1)),1) co(lcurr(temp(2)),1)],[co(lcurr(temp(1)),2) co(lcurr(temp(2)),2)],'Color',mycmap(i,:),'LineStyle','--','LineWidth',lsize,'MarkerSize',1);
    end
end
h1=plot(-100,-100,'Marker','o','MarkerSize',msize,'LineStyle','none','LineWidth',lsize,'Color','k');
h2=plot(-100,-100,'Marker','*','MarkerSize',msize,'LineStyle','none','LineWidth',lsize,'Color','k');
legend([h1 h2],{'trial 1','trial 2'})
legend boxoff
box on
axis square
set(gca,'xtick','')
set(gca,'ytick','')
axis([1.1*min(co(:,1)) 1.1*max(co(:,1)) 1.1*min(co(:,2)) 1.1*max(co(:,2))])
xlabel(['PC 1 Score (' num2str(EXPLAINEDO(1),'%0.1f') '%)'])
ylabel(['PC 2 Score (' num2str(EXPLAINEDO(2),'%0.1f') '%)'])
set(gca,'FontSize',15)


%% plot PN and ORN data in PC space, trial to trial

mycmap=hsv(nflies);
%mycmap=hsv(flyNum);

varianceToKeep=100; % percent of variance to keep


co = SCORE;
totalVarianceExplained=cumsum(EXPLAINED);
pcsWithinVariance=find(totalVarianceExplained<varianceToKeep);
pcstouse=pcsWithinVariance(end);
disp(['Using ' num2str(pcstouse) ' PCs'])

allFlies=1:(flyNum);


msize=10;
lsize=2;
% SUP FIG PN left lobe
figure %35
hold on
i=0;
%for j=1:(flyNum)
for j=pnflies  
    i=i+1;
    ltemps=find(glombyodorflies==(j));
    lcurr=ltemps(1:end);
    rtemps=find(glombyodorflies==(j+1000));
    rcurr=rtemps(1:end);
    
    
    if length(lcurr)>1
        h3=plot(co(lcurr(1),1),co(lcurr(1),2),'Color',mycmap(i,:),'Marker','o','LineWidth',lsize,'MarkerSize',msize);
        %text(mean(co(lcurr(1:size(lcurr,2)),1)),mean(co(lcurr(1:size(lcurr,2)),2)),[num2str(j)],'FontSize',15)
        h4=plot(co(lcurr(2),1),co(lcurr(2),2),'Color',mycmap(i,:),'Marker','*','LineWidth',lsize,'MarkerSize',msize);
        h5=plot(co(lcurr(1:2),1),co(lcurr(1:2),2),'--','Color',mycmap(i,:),'Marker','none','LineWidth',lsize,'MarkerSize',msize);
    end

end
h1=plot(-100,-100,'Marker','o','MarkerSize',msize,'LineStyle','none','LineWidth',lsize,'Color','k');
h2=plot(-100,-100,'Marker','*','MarkerSize',msize,'LineStyle','none','LineWidth',lsize,'Color','k');
legend([h1 h2],{'trial 1','trial 2'})
legend boxoff
box on
axis square
set(gca,'xtick','')
set(gca,'ytick','')
axis([min(co(:,1)) max(co(:,1)) min(co(:,2)) max(co(:,2))])
xlabel(['PC 1 Score (' num2str(EXPLAINED(1),'%0.1f') '%)'])
ylabel(['PC 2 Score (' num2str(EXPLAINED(2),'%0.1f') '%)'])
set(gca,'FontSize',15)

% SUP FIG PN right lobe
figure %36
hold on
i=0;
%for j=1:(flyNum)
for j=pnflies  
    i=i+1;
    ltemps=find(glombyodorflies==(j));
    lcurr=ltemps(1:end);
    rtemps=find(glombyodorflies==(j+1000));
    rcurr=rtemps(1:end);
    
    lcurr=rcurr;
    if length(lcurr)>1
        h3=plot(co(lcurr(1),1),co(lcurr(1),2),'Color',mycmap(i,:),'Marker','o','LineWidth',lsize,'MarkerSize',msize);
        %text(mean(co(lcurr(1:size(lcurr,2)),1)),mean(co(lcurr(1:size(lcurr,2)),2)),[num2str(j)],'FontSize',15)
        h4=plot(co(lcurr(2),1),co(lcurr(2),2),'Color',mycmap(i,:),'Marker','*','LineWidth',lsize,'MarkerSize',msize);
        h5=plot(co(lcurr(1:2),1),co(lcurr(1:2),2),'--','Color',mycmap(i,:),'Marker','none','LineWidth',lsize,'MarkerSize',msize);
    end

end
h1=plot(-100,-100,'Marker','o','MarkerSize',msize,'LineStyle','none','LineWidth',lsize,'Color','k');
h2=plot(-100,-100,'Marker','*','MarkerSize',msize,'LineStyle','none','LineWidth',lsize,'Color','k');
legend([h1 h2],{'trial 1','trial 2'})
legend boxoff
box on
axis square
set(gca,'xtick','')
set(gca,'ytick','')
axis([min(co(:,1)) max(co(:,1)) min(co(:,2)) max(co(:,2))])
xlabel(['PC 1 Score (' num2str(EXPLAINED(1),'%0.1f') '%)'])
ylabel(['PC 2 Score (' num2str(EXPLAINED(2),'%0.1f') '%)'])
set(gca,'FontSize',15)

% plot ORN data

co = SCOREO;
totalVarianceExplained=cumsum(EXPLAINEDO);
pcsWithinVariance=find(totalVarianceExplained<varianceToKeep);
pcstouse=pcsWithinVariance(end);
disp(['Using ' num2str(pcstouse) ' PCs'])

mycmapO=mycmap;
% SUP FIG ORN left lobe
figure %37
hold on
i=0;
%for j=1:(flyNum)
for j=ornflies  
    i=i+1;
    ltemps=find(glombyodorfliesO==(j));
    lcurr=ltemps(1:end);
    rtemps=find(glombyodorfliesO==(j+1000));
    rcurr=rtemps(1:end);
    

    if length(lcurr)>1
        h3=plot(co(lcurr(1),1),co(lcurr(1),2),'Color',mycmapO(i,:),'Marker','o','LineWidth',lsize,'MarkerSize',msize);
        %text(mean(co(lcurr(1:size(lcurr,2)),1)),mean(co(lcurr(1:size(lcurr,2)),2)),[num2str(j)],'FontSize',15)
        h4=plot(co(lcurr(2),1),co(lcurr(2),2),'Color',mycmapO(i,:),'Marker','*','LineWidth',lsize,'MarkerSize',msize);
        h5=plot(co(lcurr(1:2),1),co(lcurr(1:2),2),'--','Color',mycmapO(i,:),'Marker','none','LineWidth',lsize,'MarkerSize',msize);
    end

end
h1=plot(-100,-100,'Marker','o','MarkerSize',msize,'LineStyle','none','LineWidth',lsize,'Color','k');
h2=plot(-100,-100,'Marker','*','MarkerSize',msize,'LineStyle','none','LineWidth',lsize,'Color','k');
legend([h1 h2],{'trial 1','trial 2'})
legend boxoff
box on
axis square
set(gca,'xtick','')
set(gca,'ytick','')
axis([min(co(:,1)) max(co(:,1)) min(co(:,2)) max(co(:,2))])
xlabel(['PC 1 Score (' num2str(EXPLAINEDO(1),'%0.1f') '%)'])
ylabel(['PC 2 Score (' num2str(EXPLAINEDO(2),'%0.1f') '%)'])
set(gca,'FontSize',15)

% SUP FIG ORN right lobe
figure %38
hold on
i=0;
%for j=1:(flyNum)
for j=ornflies  
    i=i+1;
    ltemps=find(glombyodorfliesO==(j));
    lcurr=ltemps(1:end);
    rtemps=find(glombyodorfliesO==(j+1000));
    rcurr=rtemps(1:end);
    
    lcurr=rcurr;
    if length(lcurr)>1
        h3=plot(co(lcurr(1),1),co(lcurr(1),2),'Color',mycmapO(i,:),'Marker','o','LineWidth',lsize,'MarkerSize',msize);
        %text(mean(co(lcurr(1:size(lcurr,2)),1)),mean(co(lcurr(1:size(lcurr,2)),2)),[num2str(j)],'FontSize',15)
        h4=plot(co(lcurr(2),1),co(lcurr(2),2),'Color',mycmapO(i,:),'Marker','*','LineWidth',lsize,'MarkerSize',msize);
        h5=plot(co(lcurr(1:2),1),co(lcurr(1:2),2),'--','Color',mycmapO(i,:),'Marker','none','LineWidth',lsize,'MarkerSize',msize);
    end

end
h1=plot(-100,-100,'Marker','o','MarkerSize',msize,'LineStyle','none','LineWidth',lsize,'Color','k');
h2=plot(-100,-100,'Marker','*','MarkerSize',msize,'LineStyle','none','LineWidth',lsize,'Color','k');
legend([h1 h2],{'trial 1','trial 2'})
legend boxoff
box on
axis square
set(gca,'xtick','')
set(gca,'ytick','')
axis([min(co(:,1)) max(co(:,1)) min(co(:,2)) max(co(:,2))])
xlabel(['PC 1 Score (' num2str(EXPLAINEDO(1),'%0.1f') '%)'])
ylabel(['PC 2 Score (' num2str(EXPLAINEDO(2),'%0.1f') '%)'])
set(gca,'FontSize',15)
%% plots of within and across fly distances for ORN and PN
withinlobe=nanmean([withinleft(:) withinright(:)],2);
withinlobeO=nanmean([withinleftO(:) withinrightO(:)],2);

PNtoboxplot=[withinlobe];
ORNtoboxplot=[withinlobeO];
combined=[ORNtoboxplot(:); PNtoboxplot(:)];
groupings=[];
colors=[];

for i=1
    groupings=[groupings i*ones(1,flyNumO)];
    colors(i,:)=ocolor;
end
for i=2
    groupings=[groupings i*ones(1,flyNum)];
    colors(i,:)=pcolor;
end

% SUP FIG within lobe ORN/PN
figure %39
boxplot(combined,groupings,'PlotStyle','compact','Colors',colors,'symbol','')
ylabel('coding space distance')
xlabels{1}='Within Fly (Left Lobe)';
xlabels{2}='Within Fly (Right Lobe)';
xlabels{3}='Within Fly (Opposite Lobes)';

xlabels{4}='Within Fly (Left Lobe)';
xlabels{5}='Within Fly (Right Lobe)';
xlabels{6}='Within Fly (Opposite Lobes)';

set(gca,'xtick','')
set(gca,'ytick','')
box on
axis([0 3 0 7])
set(gca,'FontSize',15)

%PNtoboxplot=[withinlobe withinacross(:) acrossall(:)];
%ORNtoboxplot=[withinlobeO withinacrossO(:) acrossallO(:)];
PNtoboxplot=[withinacross(:)];
ORNtoboxplot=[withinacrossO(:)];
combined=[ORNtoboxplot(:); PNtoboxplot(:)];
groupings=[];
colors=[];

for i=1
    groupings=[groupings i*ones(1,flyNumO)];
    colors(i,:)=ocolor;
end
for i=2
    groupings=[groupings i*ones(1,flyNum)];
    colors(i,:)=pcolor;
end

% SUP FIG across lobe ORN/PN
figure %40
boxplot(combined,groupings,'PlotStyle','compact','Colors',colors,'symbol','')
ylabel('coding space distance')
xlabels{1}='Within Fly (Left Lobe)';
xlabels{2}='Within Fly (Right Lobe)';
xlabels{3}='Within Fly (Opposite Lobes)';

xlabels{4}='Within Fly (Left Lobe)';
xlabels{5}='Within Fly (Right Lobe)';
xlabels{6}='Within Fly (Opposite Lobes)';

set(gca,'xtick','')
set(gca,'ytick','')
box on
axis([0 3 0 7])
set(gca,'FontSize',15)


%PNtoboxplot=[withinlobe withinacross(:) acrossall(:)];
%ORNtoboxplot=[withinlobeO withinacrossO(:) acrossallO(:)];
PNtoboxplot=[acrossall(:)];
ORNtoboxplot=[acrossallO(:)];
combined=[ORNtoboxplot(:); PNtoboxplot(:)];
groupings=[];
colors=[];

for i=1
    groupings=[groupings i*ones(1,flyNumO)];
    colors(i,:)=ocolor;
end
for i=2
    groupings=[groupings i*ones(1,flyNum)];
    colors(i,:)=pcolor;
end

% SUP FIG across fly ORN/PN
figure %41
boxplot(combined,groupings,'PlotStyle','compact','Colors',colors,'symbol','')
ylabel('coding space distance')
xlabels{1}='Within Fly (Left Lobe)';
xlabels{2}='Within Fly (Right Lobe)';
xlabels{3}='Within Fly (Opposite Lobes)';

xlabels{4}='Within Fly (Left Lobe)';
xlabels{5}='Within Fly (Right Lobe)';
xlabels{6}='Within Fly (Opposite Lobes)';

set(gca,'xtick','')
set(gca,'ytick','')
box on
axis([0 3 0 7])
set(gca,'FontSize',15)


%% plots of within and across fly distances

% plot all data

PNtoboxplot=[withinleft(:) withinright(:) withinacross(:) acrossleft(:) acrossright(:) acrossall(:)];
ORNtoboxplot=[ withinleftO(:) withinrightO(:) withinacrossO(:) acrossleftO(:) acrossrightO(:) acrossallO(:)];

% combined=[PNtoboxplot(:); ORNtoboxplot(:)];
% groupings=[];
% colors=[];
% for i=1:6
%     groupings=[groupings i*ones(1,flyNum)];
%     colors(i,:)=ocolor;
% end
% for i=7:12
%     groupings=[groupings i*ones(1,flyNumO)];
%     colors(i,:)=pcolor;
% end

% 
% figure
% boxplot(combined,groupings,'PlotStyle','compact','Colors',colors)
% ylabel('coding space distance')
% xlabels{1}='within (left)';
% xlabels{2}='within (right)';
% xlabels{3}='within (opposite)';
% 
% xlabels{4}='across (left)';
% xlabels{5}='across (right)';
% xlabels{6}='across (both)';
% 
% xlabels{7}='within (left)';
% xlabels{8}='within (right)';
% xlabels{9}='within (opposite)';
% 
% xlabels{10}='across (left)';
% xlabels{11}='across (right)';
% xlabels{12}='across (both)';
% 
% set(gca,'xtick',1:12,'xticklabel',xlabels,'FontSize',10)
% xtickangle(30)
% box off
% set(gca,'FontSize',15)

% Figure 42
violinPlot(ORNtoboxplot,ocolor)
ylabel('coding space distance')
xlabels{1}='within (left)';
xlabels{2}='within (right)';
xlabels{3}='within (opposite)';
xlabels{4}='across (left)';
xlabels{5}='across (right)';
xlabels{6}='across (both)';
set(gca,'xtick',1:6,'xticklabel',xlabels,'FontSize',10)
xtickangle(30)
box off
set(gca,'FontSize',15)

% Figure 43
violinPlot(PNtoboxplot,pcolor)
ylabel('coding space distance')
xlabels{1}='within (left)';
xlabels{2}='within (right)';
xlabels{3}='within (opposite)';
xlabels{4}='across (left)';
xlabels{5}='across (right)';
xlabels{6}='across (both)';
set(gca,'xtick',1:6,'xticklabel',xlabels,'FontSize',10)
xtickangle(30)
box off
set(gca,'FontSize',15)



PNtoboxplot=[withinleft(:) withinright(:) withinacross(:)];
ORNtoboxplot=[ withinleftO(:) withinrightO(:) withinacrossO(:)];

% combined=[PNtoboxplot(:); ORNtoboxplot(:)];
% groupings=[];
% colors=[];
% for i=1:3
%     groupings=[groupings i*ones(1,flyNum)];
%     colors(i,:)=ocolor;
% end
% for i=4:6
%     groupings=[groupings i*ones(1,flyNumO)];
%     colors(i,:)=pcolor;
% end

% figure
% boxplot(combined,groupings,'PlotStyle','compact','Colors',colors)
% ylabel('Distance in Coding Space')
% xlabels{1}='Within Fly (Left Lobe)';
% xlabels{2}='Within Fly (Right Lobe)';
% xlabels{3}='Within Fly (Opposite Lobes)';
% 
% xlabels{4}='Within Fly (Left Lobe)';
% xlabels{5}='Within Fly (Right Lobe)';
% xlabels{6}='Within Fly (Opposite Lobes)';
% 
% set(gca,'xtick',1:6,'xticklabel',xlabels,'FontSize',10)
% xtickangle(30)
% box off
% set(gca,'FontSize',15)

% Figure 44
violinPlot(ORNtoboxplot,ocolor)
ylabel('coding space distance')
xlabels{1}='within (left)';
xlabels{2}='within (right)';
xlabels{3}='within (opposite)';
set(gca,'xtick',1:3,'xticklabel',xlabels,'FontSize',10)
xtickangle(30)
box off
set(gca,'FontSize',15)

% Figure 45
violinPlot(PNtoboxplot,pcolor)
ylabel('coding space distance')
xlabels{1}='within (left)';
xlabels{2}='within (right)';
xlabels{3}='within (opposite)';
set(gca,'xtick',1:3,'xticklabel',xlabels,'FontSize',10)
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
% 
% combined=[PNtoboxplot(:); ORNtoboxplot(:)];
% groupings=[];
% colors=[];
% for i=1:3
%     groupings=[groupings i*ones(1,flyNum)];
%     colors(i,:)=ocolor;
% end
% for i=4:6
%     groupings=[groupings i*ones(1,flyNumO)];
%     colors(i,:)=pcolor;
% end
% 
% figure
% boxplot(combined,groupings,'PlotStyle','compact','Colors',colors)
% ylabel('Distance in Coding Space')
% xlabels{1}='Across Fly (Left Lobe)';
% xlabels{2}='Across Fly (Right Lobe)';
% xlabels{3}='Across Fly (Both Lobes)';
% 
% xlabels{4}='Across Fly (Left Lobe)';
% xlabels{5}='Across Fly (Right Lobe)';
% xlabels{6}='Across Fly (Both Lobes)';
% 
% 
% set(gca,'xtick',1:6,'xticklabel',xlabels,'FontSize',10)
% xtickangle(30)
% box off
% set(gca,'FontSize',15)

% Figure 46
violinPlot(ORNtoboxplot,ocolor)
ylabel('coding space distance')
xlabels{1}='across (left)';
xlabels{2}='across (right)';
xlabels{3}='across (both)';
set(gca,'xtick',1:3,'xticklabel',xlabels,'FontSize',10)
xtickangle(30)
box off
set(gca,'FontSize',15)

% Figure 47
violinPlot(PNtoboxplot,pcolor)
ylabel('coding space distance')
xlabels{1}='across (left)';
xlabels{2}='across (right)';
xlabels{3}='across (both)';
set(gca,'xtick',1:3,'xticklabel',xlabels,'FontSize',10)
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
PNw=[withinleft(:) withinright(:) withinacross(:)]; PNw=PNw(:);
PNa=[acrossleft(:) acrossright(:) acrossall(:)]; PNa=PNa(:);
ORNw=[withinleftO(:) withinrightO(:) withinacrossO(:)]; ORNw=ORNw(:);
ORNa=[acrossleftO(:) acrossrightO(:) acrossallO(:)]; ORNa=ORNa(:);


PNtoboxplot=[nanmean([withinleft(:) withinright(:) withinacross(:)],2) nanmean([acrossleft(:) acrossright(:) acrossall(:)],2)];
ORNtoboxplot=[nanmean([withinleftO(:) withinrightO(:) withinacrossO(:)],2) nanmean([acrossleftO(:) acrossrightO(:) acrossallO(:)],2)];

combined=[ORNtoboxplot(:); PNtoboxplot(:)];

groupings=[];
groupings=[groupings 1*ones(1,flyNumO)];
groupings=[groupings 3*ones(1,flyNumO)];


groupings=[groupings 2*ones(1,flyNum)];
groupings=[groupings 4*ones(1,flyNum)];


colors=[ocolor; pcolor; ocolor; pcolor];
% FIG Euclidean distance in full coding space
figure %48
boxplot(combined,groupings,'PlotStyle','compact','Colors',colors,'symbol','')
%boxplot(combined,groupings,'BoxStyle','filled','Colors',colors)
ylabel('coding space distance')
xlabels{1}='within';
xlabels{2}='across';
set(gca,'xtick','')
set(gca,'ytick','')
axis([0 5 0 7])
xtickangle(30)
box on
set(gca,'FontSize',15)

d = cell(1,4);
d{1}=ORNw;
d{2}=PNw;
d{3}=ORNa;
d{4}=PNa;
figure
violinPlot(d,[ocolor;pcolor;ocolor;pcolor])
%boxplot(combined,groupings,'BoxStyle','filled','Colors',colors)
ylabel('coding space distance')
xlabels{1}='within';
xlabels{2}='across';
set(gca,'xtick',[1.5 3.5],'xticklabel',xlabels,'FontSize',10)
xtickangle(30)
box off
set(gca,'FontSize',15)

%[s PNwithinVacross] = ttest2(combined(groupings==1),combined(groupings==3))
%[s ORNwithinVacross] = ttest2(combined(groupings==2),combined(groupings==4))
%[s acrossFlyORNvsPN] = ttest2(combined(groupings==3),combined(groupings==4))
%[s withinFlyORNvsPN] = ttest2(combined(groupings==1),combined(groupings==2))



