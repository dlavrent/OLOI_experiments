
clear all
close all
rng('default')

load ORN_PN_colors
load analysis_dir_path
manualLabelHome=fullfile(analysis_dir_path, 'PN_analysis/yokedControlExperiment_yokedControlFlies');

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
responsesTensorFullyOrthogonal=NaN*zeros(27,2,39,nodors,length(odortimesfortimecourse)); % flynum HARD CODED
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
        currFlyNum=currname((underscores(4)+4):(underscores(5)-1));
        
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
        %behaviorOcc=[behaviorOcc occ];
        behaviorpreOcc=[behaviorpreOcc preocc];

        manualClusterLabels=clusterLabels;
        totalPutativeGloms(i)=length(manualClusterLabels);
        
        
        gs=median(grnResponse(:,:,odortimes),3); % use median
        gs=zeros(13,size(grnResponse,2));
        for j=1:size(grnResponse,2)
            temp = squeeze(grnResponse(:,j,odortimes)-nanmedian(grnResponse(:,j,1:5),3));
            gs(:,j)= transpose(max(temp'));
        end
        %gs=prctile(grnResponse(:,:,odortimes),75); % use percentile
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


% shuffle data points that correspond to each fly
flyindicesShuffled=cell(1,flyNum);
j=1;
for i=randperm(flyNum)
    [temp temp2]=find(glombyodorflies==i);
    [tem3 temp4]=find(glombyodorflies==(i+leftRightSplitVal));
    flyindicesShuffled{j}=[temp2 temp4];
    j=j+1;
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

%% load experimental flies

manualLabelHome=fullfile(analysis_dir_path, 'PN_analysis/yokedControlExperiment_experimentalFlies');

publishedOdorPath=fullfile(analysis_dir_path, 'utilities/odorPanel_12_DoORData.mat');
load(publishedOdorPath);

manualLabelledFolders=dir(manualLabelHome);
manualLabelledFolders=manualLabelledFolders(3:end);

labels=cell(1,length(manualLabelledFolders));
glomsannotatedE=zeros(1,length(manualLabelledFolders));
totalPutativeGloms=zeros(1,length(manualLabelledFolders));

nodors=13;
odortimes=[6:9]; % hard-coded specific time interval for summarizing odor response
odortimesfortimecourse=[1:18];
responsesGlomByOdorE=[];
responsesTimeCourseGlomByOdorE=[];
flyodorsE=[];
flyglomsE=[];
glombyodorfliesE=[];
behaviorOccE=[];
behaviorpreOccE=[];
leftRightSplitVal=1000;
firstSecondPanelSplitVal=500;
flyNumE=0;
glomfoundE=zeros(1,39);

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
            flyNumE=flyNumE+1;
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

        behaviorOccE=[behaviorOccE occ-preocc];
        %behaviorOcc=[behaviorOcc occ];
        behaviorpreOccE=[behaviorpreOccE preocc];

        manualClusterLabels=clusterLabels;
        totalPutativeGloms(i)=length(manualClusterLabels);
        
        
        gs=median(grnResponse(:,:,odortimes),3); % use median
        gs=zeros(13,size(grnResponse,2));
        for j=1:size(grnResponse,2)
            temp = squeeze(grnResponse(:,j,odortimes)-nanmedian(grnResponse(:,j,1:5),3));
            gs(:,j)= transpose(max(temp'));
        end
        %gs=prctile(grnResponse(:,:,odortimes),75); % use percentile
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
                    glomsannotatedE(i)=glomsannotatedE(i)+1;
                    glomfoundE(k)=glomfoundE(k)+1;
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
        
        responsesGlomByOdorE=[responsesGlomByOdorE responseTempT(:)];
        responsesTimeCourseGlomByOdorE=[responsesTimeCourseGlomByOdorE responseTempTimeCourseT(:)];
        glombyodorfliesE=[glombyodorfliesE (flyNumE+(cL-1)*leftRightSplitVal)];
        
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
flyindicesE=cell(1,flyNumE);
flyindicesEL=cell(1,flyNumE);
flyindicesER=cell(1,flyNumE);
for i=1:flyNumE
    [temp temp2]=find(glombyodorfliesE==i);
    [tem3 temp4]=find(glombyodorfliesE==(i+leftRightSplitVal));
    flyindicesE{i}=[temp2 temp4];
    flyindicesEL{i}=temp2;
    flyindicesER{i}=temp4;
end


% shuffle data points that correspond to each fly
flyindicesShuffledE=cell(1,flyNumE);
j=1;
for i=randperm(flyNumE)
    [temp temp2]=find(glombyodorfliesE==i);
    [tem3 temp4]=find(glombyodorfliesE==(i+leftRightSplitVal));
    flyindicesShuffledE{j}=[temp2 temp4];
    j=j+1;
end



%% perform pca on yoked controls

clear responsesNoResponseRemoved

fracIn=0.01; % best results when fracIn is high, ~0.5, only using high confidence glomeruli


medianResponseOrTimeCourse=1; % 1 for median response only, 0 for full time course

if medianResponseOrTimeCourse
    responsesNoResponseRemoved=responsesGlomByOdor;
else
    responsesNoResponseRemoved=responsesTimeCourseGlomByOdor;
end    


% remove fully empty rows 
%totallyempty=~any(isfinite(responsesNoResponseRemoved),2);
%responsesNoResponseRemoved(totallyempty,:)=[];

gNames=publishedOR.gh146glomerulusNames;
glomsFound=glomfound;
numFinite=sum(isfinite(responsesNoResponseRemoved),2);
toRemove=find(numFinite/size(responsesNoResponseRemoved,2)<=fracIn);
responsesNoResponseRemoved(toRemove,:)=[];

% remove air and ethanol
%responsesNoResponseRemoved(1:13:end)=0;
%responsesNoResponseRemoved(8:13:end)=0;

if medianResponseOrTimeCourse
    temp=nodors;
    fp=toRemove(find(mod(toRemove,temp)==1));
    glomsremoved=((fp-1)/temp)+1;
    gNames(glomsremoved)=[];
    glomsFound(glomsremoved)=[];
else
    temp=length(odortimesfortimecourse)*nodors;
    fp=toRemove(find(mod(toRemove,temp)==1));
    glomsremoved=((fp-1)/temp)+1;
    gNames(glomsremoved)=[];
end
gNames(1)=[];
responsesNoResponseRemoved(1:13,:)=[];
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

responsesNoResponseRemoved=responsesNoResponseRemoved';
responsesNoResponseRemoved=(responsesNoResponseRemoved-mean(responsesNoResponseRemoved));
responsesNoResponseRemoved=responsesNoResponseRemoved';


opt = statset('pca');
opt.Display='iter';
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(responsesNoResponseRemoved','Options',opt);


figure; %1
plot(cumsum(EXPLAINED),'o-','LineWidth',3)
ylabel('Variance Explained (%)')
xlabel('PC #')
box off
set(gca,'FontSize',15)


%% perform pca on experimental flies controls

clear responsesNoResponseRemovedE

fracIn=0.2; % best results when fracIn is high, ~0.5, only using high confidence glomeruli


medianResponseOrTimeCourse=1; % 1 for median response only, 0 for full time course

if medianResponseOrTimeCourse
    responsesNoResponseRemovedE=responsesGlomByOdorE;
else
    responsesNoResponseRemovedE=responsesTimeCourseGlomByOdorE;
end    


% remove fully empty rows 
%totallyempty=~any(isfinite(responsesNoResponseRemoved),2);
%responsesNoResponseRemoved(totallyempty,:)=[];

gNamesE=publishedOR.gh146glomerulusNames;
glomsFoundE=glomfoundE;
numFinite=sum(isfinite(responsesNoResponseRemovedE),2);
toRemove=find(numFinite/size(responsesNoResponseRemovedE,2)<=fracIn);
responsesNoResponseRemovedE(toRemove,:)=[];

% remove air and ethanol
%responsesNoResponseRemoved(1:13:end)=0;
%responsesNoResponseRemoved(8:13:end)=0;

if medianResponseOrTimeCourse
    temp=nodors;
    fp=toRemove(find(mod(toRemove,temp)==1));
    glomsremovedE=((fp-1)/temp)+1;
    gNamesE(glomsremovedE)=[];
    glomsFoundE(glomsremovedE)=[];
end

% % fill nans with mean
for i=1:size(responsesNoResponseRemovedE,1)
    for j=1:size(responsesNoResponseRemovedE,2)
        if isnan(responsesNoResponseRemovedE(i,j))
            responsesNoResponseRemovedE(i,j)=nanmean(responsesNoResponseRemovedE(i,:));
        end
    end
end

% remove D
responsesNoResponseRemovedE(1:nodors,:)=[];
gNamesE(1)=[];


responsesNoResponseRemovedE=responsesNoResponseRemovedE';
responsesNoResponseRemovedE=(responsesNoResponseRemovedE-mean(responsesNoResponseRemovedE));
responsesNoResponseRemovedE=responsesNoResponseRemovedE';


opt = statset('pca');
opt.Display='iter';
[COEFFE, SCOREE, LATENTE, TSQUAREDE, EXPLAINEDE] = pca(responsesNoResponseRemovedE','Options',opt);


figure; %2
plot(cumsum(EXPLAINEDE),'o-','LineWidth',3)
ylabel('Variance Explained (%)')
xlabel('PC #')
box off
set(gca,'FontSize',15)


%% apply training data predictor to yoked control data
trainedModel = load('trainDataModel.mat');

matchedResponses=responsesNoResponseRemoved;
COEFFmatched=COEFF;
trainedPC=trainedModel.mypc;
trainedGNames=trainedModel.gNames;
% match available gloms in test data with gloms from trainedPC
testGlomInTrainGlom=zeros(1,length(gNames));
for i=1:length(gNames)
    for j=1:length(trainedGNames)
        if strcmp(gNames{i},trainedGNames{j})
            testGlomInTrainGlom(i)=1;
            break
        end
    end
end
% remove data for test gloms not in train set
todeleteTest=[];
for i=1:length(gNames)
   if ~testGlomInTrainGlom(i)
       todeleteTest=[todeleteTest ((nodors*(i-1)+1):(nodors*i))];
   end
end
matchedResponses(todeleteTest,:)=[];


trainGlomInTestGlom=zeros(1,length(trainedGNames));
for i=1:length(trainedGNames)
    for j=1:length(gNames)
        if strcmp(gNames{j},trainedGNames{i})
            trainGlomInTestGlom(i)=1;
            break
        end
    end
end
% remove data for train gloms PC not in test set
todeleteTrain=[];
for i=1:length(trainedGNames)
   if ~trainGlomInTestGlom(i)
       todeleteTrain=[todeleteTrain ((nodors*(i-1)+1):(nodors*i))];
   end
end
trainedPC(todeleteTrain)=[];

COEFFmatched(todeleteTest,:)=[];
% measure correlation between the trained pc and test data PC
pccorr=zeros(1,10);
for i=1:10
    tempr=corr(trainedPC,COEFFmatched(:,i));
    pccorr(i)=tempr;
end

behaviorprediction=matchedResponses'*trainedPC;

flyTruePref=zeros(1,flyNum);
flyPredictedPref=zeros(1,flyNum);
ally=behaviorOcc';

nactivity=zeros(flyNum,1);
for i=1:flyNum
    flyTruePref(i)=mean(ally(flyindices{i}));
    flyPredictedPref(i)=mean(mean(behaviorprediction(flyindices{i},:)));
    nactivity(i,:)=mean(behaviorprediction(flyindices{i},:));
end
linmodel=trainedModel.linmodelPrecorrected;
myprediction=predict(linmodel,nactivity);
figure %3
plot(myprediction,flyTruePref,'.','Color',pcolor,'LineWidth',3)
hold on
xlabel('Predicted Preference')
ylabel('Measured Preference')
box on
axis([-.6 .3 -.6 .3])
axis square
set(gca,'FontSize',15)


corrcoef(myprediction,flyTruePref)

% FIG 2k
figure %4
hold on;
xVals = (myprediction-mean(myprediction))/(std(myprediction));
yVals = (flyTruePref-mean(flyTruePref))/(std(flyTruePref)); 
linreg = linearRegressionCI2(xVals, yVals.', 1, 0, -2.6, 2.2);

areaBar(linreg.xVals,polyval(linreg.pOverall,linreg.xVals),2*std(linreg.fits),[0 0 0],[0.9 0.9 0.9])
plot(xVals,yVals,'.','Color',pcolor, 'LineWidth',3,'MarkerSize',20)
for i=1:flyNum
   hold on
   %text(myprediction(i)+0.01,flyTruePref(i),num2str(i),'FontSize',15)
end
[r p]=corrcoef(xVals,yVals);
text(-.25,-.25,['r = ' num2str(r(1,2),'%2.2f')],'FontSize',15)
text(-.3,-.25,['p = ' num2str(p(1,2),'%2.2f')],'FontSize',15)
set(gca,'FontSize',15)
box on
axis([-2.6 2.2 -2.6 2.2])
axis square
xlabel('predicted preference (z-scored)')
ylabel('measured preference (z-scored)')
set(gca,'xtick','')
set(gca,'ytick','')

%% apply trained preditor to experimental flies

matchedResponses=responsesNoResponseRemovedE;
COEFFmatched=COEFFE;
trainedPC=trainedModel.mypc;
trainedGNames=trainedModel.gNames;
% match available gloms in test data with gloms from trainedPC
testGlomInTrainGlom=zeros(1,length(gNames));
for i=1:length(gNamesE)
    for j=1:length(trainedGNames)
        if strcmp(gNamesE{i},trainedGNames{j})
            testGlomInTrainGlom(i)=1;
            break
        end
    end
end
% remove data for test gloms not in train set
todeleteTest=[];
for i=1:length(gNamesE)
   if ~testGlomInTrainGlom(i)
       todeleteTest=[todeleteTest ((nodors*(i-1)+1):(nodors*i))];
   end
end
matchedResponses(todeleteTest,:)=[];


trainGlomInTestGlom=zeros(1,length(trainedGNames));
for i=1:length(trainedGNames)
    for j=1:length(gNamesE)
        if strcmp(gNamesE{j},trainedGNames{i})
            trainGlomInTestGlom(i)=1;
            break
        end
    end
end
% remove data for train gloms PC not in test set
todeleteTrain=[];
for i=1:length(trainedGNames)
   if ~trainGlomInTestGlom(i)
       todeleteTrain=[todeleteTrain ((nodors*(i-1)+1):(nodors*i))];
   end
end
trainedPC(todeleteTrain)=[];

COEFFmatched(todeleteTest,:)=[];
% measure correlation between the trained pc and test data PC
pccorr=zeros(1,10);
for i=1:10
    tempr=corr(trainedPC,COEFFmatched(:,i));
    pccorr(i)=tempr;
end

behaviorprediction=matchedResponses'*trainedPC;

flyTruePref=zeros(1,flyNumE);
flyPredictedPref=zeros(1,flyNumE);
ally=behaviorOccE';

nactivity=zeros(flyNumE,1);
for i=1:flyNumE
    flyTruePref(i)=mean(ally(flyindicesE{i}));
    flyPredictedPref(i)=mean(mean(behaviorprediction(flyindicesE{i},:)));
    nactivity(i,:)=mean(behaviorprediction(flyindicesE{i},:));
end
linmodel=trainedModel.linmodelPrecorrected;
myprediction=predict(linmodel,nactivity);
figure %5
plot(myprediction,flyTruePref,'.','Color',pcolor,'LineWidth',3)
hold on
xlabel('Predicted Preference')
ylabel('Measured Preference')
box on
axis([-.6 .3 -.6 .3])
axis square
set(gca,'FontSize',15)


corrcoef(myprediction,flyTruePref)


figure %6
plot((myprediction-mean(myprediction))/(std(myprediction)),(flyTruePref-mean(flyTruePref))/(std(flyTruePref)),'.','Color',pcolor, 'LineWidth',3)
for i=1:flyNum
   hold on
   %text(myprediction(i)+0.01,flyTruePref(i),num2str(i),'FontSize',15)
end
[r p]=corrcoef((myprediction-mean(myprediction))/(std(myprediction)),(flyTruePref-mean(flyTruePref))/(std(flyTruePref)));
text(-.25,-.25,['r = ' num2str(r(1,2),'%2.2f')],'FontSize',15)
text(-.3,-.25,['p = ' num2str(p(1,2),'%2.2f')],'FontSize',15)
set(gca,'FontSize',15)
box on
%axis([-2.1 2.1 -2.1 2.1])
axis square
xlabel('predicted preference (z-scored)')
ylabel('measured preference (z-scored)')
%set(gca,'xtick','')
%set(gca,'ytick','')






