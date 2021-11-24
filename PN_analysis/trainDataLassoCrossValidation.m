% look at manually labelled gh146 flies in odor space
% load all data
% only use 80% of data to a) determine principal components and b) fit a
% linear prediction to behavior
% use the remaining 20% of data, apply the predictors fit with training
% data, and measure how well the predictors predict behavior
clear all
close all
rng('default')

manualLabelHome='/Users/mattchurgin/Dropbox/flyimaging/analysis/PN_analysis/train';

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
        %behaviorOcc=[behaviorOcc occ];
        behaviorpreOcc=[behaviorpreOcc preocc];
        
        
        manualClusterLabels=clusterLabels;
        totalPutativeGloms(i)=length(manualClusterLabels);
        
        
        gs=median(grnResponse(:,:,odortimes),3); % use median
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



%% run train/test many times and save replicate data
warning('off')

iters=100;
yesZscore=1; % zscore data matrix before performing PCA (0=no, 1=yes)
lasso_alpha=1E-5;

fracIn=0.5; % best results when fracIn is high, ~0.5, only using high confidence glomeruli
medianResponseOrTimeCourse=1; % 1 for median response only, 0 for full time course

if medianResponseOrTimeCourse
    responsesNoResponseRemoved=responsesGlomByOdor;
else
    responsesNoResponseRemoved=responsesTimeCourseGlomByOdor;
end

% blank odors (when solenoids failed for example)
if medianResponseOrTimeCourse
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
else
    
end
gNames=publishedOR.gh146glomerulusNames;
glomsFound=glomfound;
numFinite=sum(isfinite(responsesNoResponseRemoved),2);
toRemove=find(numFinite/size(responsesNoResponseRemoved,2)<=fracIn);
responsesNoResponseRemoved(toRemove,:)=[];

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

predictors=zeros(length(gNames)*nodors,iters);
mycorr=zeros(1,iters);
mycorrshuffled=zeros(1,iters);
for jjj=1:iters
    
    % Split data into randomly assigned train/test sets
    testsize=40;
    trainsize=100-testsize;
    randomizedOrder=randperm(flyNum);
    holdoutflies=randomizedOrder(1:round((length(randomizedOrder)*testsize/100)));
    trainflies=setxor(1:flyNum,holdoutflies);
    traintoremove=[];
    for i=1:length(holdoutflies)
        temp=find(glombyodorflies==(holdoutflies(i)));
        temp2=find(glombyodorflies==(holdoutflies(i)+leftRightSplitVal));
        
        traintoremove=[traintoremove temp temp2];
    end
    
    testtoremove=setxor(1:length(glombyodorflies),traintoremove);
    
    flyindicestrain=cell(1,length(trainflies));
    glombyodorfliestrain=glombyodorflies;
    glombyodorfliestrain(traintoremove)=[];
    for i=1:length(trainflies)
        [temp temp2]=find(glombyodorfliestrain==trainflies(i));
        [tem3 temp4]=find(glombyodorfliestrain==(trainflies(i)+leftRightSplitVal));
        flyindicestrain{i}=[temp2 temp4];
    end
    
    flyindicestest=cell(1,length(holdoutflies));
    glombyodorfliestest=glombyodorflies;
    glombyodorfliestest(testtoremove)=[];
    for i=1:length(holdoutflies)
        [temp temp2]=find(glombyodorfliestest==holdoutflies(i));
        [tem3 temp4]=find(glombyodorfliestest==(holdoutflies(i)+leftRightSplitVal));
        flyindicestest{i}=[temp2 temp4];
    end
    
    
    % fill nans with mean
    
    for i=1:size(responsesNoResponseRemoved,1)
        for j=1:size(responsesNoResponseRemoved,2)
            if isnan(responsesNoResponseRemoved(i,j))
                responsesNoResponseRemoved(i,j)=nanmean(responsesNoResponseRemoved(i,:));
            end
        end
    end
    
    testdata=responsesNoResponseRemoved;
    testdata(:,testtoremove)=[];
    
    traindata=responsesNoResponseRemoved;
    traindata(:,traintoremove)=[];
    
    behaviorData=behaviorOcc;
    trainBehavior=behaviorData;
    trainBehavior(traintoremove)=[];
    testBehavior=behaviorData;
    testBehavior(testtoremove)=[];
    
    if yesZscore
        % z-score train and test data
        
        ztrain=traindata';
        ztrain=(ztrain-mean(ztrain))./(std(ztrain));
        ztrain=ztrain';
        
        ztest=testdata';
        ztest=(ztest-mean(ztest))./(std(ztest));
        ztest=ztest';
        
        traindata=ztrain;
        testdata=ztest;
    end
    
    meanTrainData=zeros(size(traindata,1),length(flyindicestrain));
    meanTestData=zeros(size(traindata,1),length(flyindicestest));
    flyTrainBehavior=zeros(1,length(flyindicestrain));
    flyTestBehavior=zeros(1,length(flyindicestest));

    for i=1:length(flyindicestrain)
       meanTrainData(:,i)=mean(traindata(:,flyindicestrain{i}),2);
       flyTrainBehavior(i)=mean(trainBehavior(flyindicestrain{i}));
    end
    
    for i=1:length(flyindicestest)
       meanTestData(:,i)=mean(testdata(:,flyindicestest{i}),2);
       flyTestBehavior(i)=mean(testBehavior(flyindicestest{i}));
    end
    
    
    B = lasso(meanTrainData',flyTrainBehavior,'Alpha',lasso_alpha,'CV',10);
    mypredictor=B(:,1);
    
    myprediction=meanTestData'*mypredictor;
    [r p]=corrcoef(myprediction,flyTestBehavior);
    
    [rshuffled pshuffled]=corrcoef(myprediction,flyTestBehavior(randperm(length(flyTestBehavior))));
    
    
    mycorr(jjj)=r(1,2);
    mycorrshuffled(jjj)=rshuffled(1,2);
    predictors(:,jjj)=mypredictor;
    if mod(jjj,10)==0
        disp(['iteration ' num2str(jjj)])
    end

end
rsq=mycorr.^2;
rsqshuffled=mycorrshuffled.^2;

figure;
distributionPlot([rsq' rsqshuffled'],'histOpt',1,'colormap',1-gray(64),'showMM',0);
ylabel('R^2')
set(gca,'xticklabel',[{'unshuffled'} {'shuffled'}])
xtickangle(30)
box off
set(gca,'FontSize',15)


mypredictor=mean(predictors,2);

figure;
plot(mypredictor,'*','LineWidth',2,'MarkerSize',8)
hold on
plot(zeros(1,length(mypredictor)),'k--','LineWidth',3)
j=1;
for i=1:nodors:length(mypredictor)
   plot((i-0.5)*ones(1,5), linspace(min(mypredictor),max(mypredictor),5),'k--','LineWidth',2)
   
   j=j+1;
end
ylabel('lasso coefficient')
set(gca,'xtick',(1:nodors:length(mypredictor))+floor(nodors/2),'xticklabel',string(gNames),'FontSize',10)
xtickangle(30)
box off
set(gca,'FontSize',15)


% save lasso coefficients
%save(['190220_lassoCoeffs.mat'],'mypredictor','gNames')