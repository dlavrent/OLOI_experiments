% keep or remove single glomeruli to determine importance for model
clear all
close all
rng('default')
allColorMaps
manualLabelHome='/Users/mattchurgin/Dropbox/flyimaging/analysis/PN_analysis/added181218plus2019';

publishedOdorPath='/Users/mattchurgin/Desktop/dblab/mattFunctions/odorpanelResponsesDoOR/odorPanel_12/odorPanel_12_DoORData.mat';
load(publishedOdorPath);
odornames=[{'air'},{'3-octanol'},{'1-hexanol'},{'ethyl lactate'},{'citronella'},{'2-heptanone'},{'1-pentanol'},{'ethanol'},{'geranyl acetate'},{'hexyl acetate'},{'4-methylcyclohexanol'},{'pentyl acetate'},{'1-butanol'}];
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



%% keep single glomerulus
warning('off')

iters=100;
highestPCtouse=5;
yesZscore=0; % zscore data matrix before performing PCA (0=no, 1=yes)
normalizeResponse=0;

fracIn=0.4; % best results when fracIn is high, ~0.5, only using high confidence glomeruli
medianResponseOrTimeCourse=1; % 1 for median response only, 0 for full time course

nGloms=6;

testR=zeros(highestPCtouse,iters,nGloms);
testRshuffled=zeros(highestPCtouse,iters,nGloms);

for glomKeep=1:nGloms
    disp(['keeping only glomerulus ' num2str(glomKeep)])
    if medianResponseOrTimeCourse
        responsesNoResponseRemoved=responsesGlomByOdor;
    else
        responsesNoResponseRemoved=responsesTimeCourseGlomByOdor;
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
    
    % remove all but one glomerulus
    for removeGloms=1:nGloms
        if removeGloms~=glomKeep
            responsesNoResponseRemoved(((removeGloms-1)*nodors+1):((removeGloms)*nodors),:)=0;
        end
    end
    
    
    % normalize responsesNoResponseRemoved
    if normalizeResponse
        for i=1:size(responsesNoResponseRemoved,2)
            responsesNoResponseRemoved(:,i)=responsesNoResponseRemoved(:,i)/nanmean(abs(responsesNoResponseRemoved(:,i)));
        end
    end
    
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
        
        
        % % fill nans with mean
        
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
        
        [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(traindata');
        co = SCORE;
        
        if jjj==1
           rawPC{glomKeep}=zeros(nodors,iters,highestPCtouse);
        end

        for pcstouse=1:highestPCtouse
            currpc=COEFF(:,pcstouse);
            
            rawPC{glomKeep}(:,jjj,pcstouse)=currpc(((glomKeep-1)*nodors+1):((glomKeep)*nodors))';
            %rawPC{glomKeep}(:,jjj,pcstouse)=currpc;
            
            % generate linear model based on current PC
            behaviorprediction=(traindata'*currpc);
            flyTruePref=zeros(1,length(trainflies));
            flyTruePrefShuffled=zeros(1,length(trainflies));
            ally=behaviorOcc';
            
            ally(traintoremove)=[];
            
            nactivity=zeros(length(trainflies),length(pcstouse));
            shuffledindtrain=randperm(length(trainflies)); % for comparing to shuffled null model
            for i=1:length(trainflies)
                flyTruePref(i)=mean(ally(flyindicestrain{i}));
                flyTruePrefShuffled(i)=mean(ally(flyindicestrain{shuffledindtrain(i)}));
                nactivity(i,:)=mean(behaviorprediction(flyindicestrain{i},:));
            end
            linmodel=fitlm(nactivity,flyTruePref);
            
            beta=linmodel.Coefficients.Estimate;
            
            PCContribution=currpc*beta(2:end);
            
            % apply predictor to testdata and evaluate how held out points fit in
            % model
            
            %behaviorprediction=testdata'*PCContribution;
            behaviorprediction=testdata'*currpc;
            testflyTruePref=zeros(1,length(holdoutflies));
            testflyTruePrefShuffled=zeros(1,length(holdoutflies));
            flyPredictedPref=zeros(1,length(holdoutflies));
            ally=behaviorOcc';
            
            ally(testtoremove)=[];
            
            testnactivity=zeros(length(holdoutflies),1);
            shuffledindtest=randperm(length(holdoutflies)); % for comparing to shuffled null model
            for i=1:length(holdoutflies)
                testflyTruePref(i)=mean(ally(flyindicestest{i}));
                testflyTruePrefShuffled(i)=mean(ally(flyindicestest{shuffledindtest(i)}));
                testnactivity(i,:)=mean(behaviorprediction(flyindicestest{i},:));
            end
            
            mytestprediction=predict(linmodel,testnactivity);
            
            myr=corrcoef(mytestprediction,testflyTruePref);
            myrshuffled=corrcoef(mytestprediction,testflyTruePref(randperm(length(testflyTruePref))));
            
            testR(pcstouse,jjj,glomKeep)=myr(1,2);
            testRshuffled(pcstouse,jjj,glomKeep)=myrshuffled(1,2);
            
        end
        
        
        if mod(jjj,10)==0
            disp(['iteration ' num2str(jjj)])
        end
    end
end


myr2=(testR.^2).*sign(testR);
myr2shuffled=(testRshuffled.^2).*sign(testRshuffled);

for i=1:highestPCtouse
    pc2=squeeze(myr2(i,:,:));
    pc2shuffled=squeeze(myr2shuffled(i,:,:));
    
    figure
    subplot(1,2,1)
    distributionPlot(pc2,'histOpt',1,'colormap',1-gray(64),'showMM',0)
    xlabel('glomerulus used for linear model')
    ylabel('Unshuffled R^2')
    set(gca,'XTickLabel',gNames)
    xtickangle(30)
    axis([0 length(gNames)+1 -0.8 0.8])
    box off
    set(gca,'FontSize',15)
    subplot(1,2,2)
    distributionPlot(pc2shuffled,'histOpt',1,'colormap',1-gray(64),'showMM',0)
    xlabel('odor used for linear model')
    ylabel('Shuffled R^2')
    set(gca,'XTickLabel',gNames)
    xtickangle(30)
    axis([0 length(gNames)+1 -0.8 0.8])
    box off
    set(gca,'FontSize',15)
    title(['PC ' num2str(i)])
end



figure
imagesc(squeeze(mean(myr2,2)),[0 0.3])
set(gca,'XTick',[1:nodors])
set(gca,'XTickLabel',gNames)
xtickangle(30)
ylabel('PC #')
xlabel('glomerulus used for model')
set(gca,'YTick',[1:highestPCtouse])
h=colorbar;
title(h,'R^2')
set(gca,'FontSize',15)

% find best pc for each odor
averageR2=squeeze(mean(myr2,2));
[bestR2 bestPCforodor]=max(averageR2,[],1);

for j=1:nGloms
pc2(:,j)=squeeze(myr2(bestPCforodor(j),:,j));
pc2shuffled(:,j)=squeeze(myr2shuffled(bestPCforodor(j),:,j));
end

figure
subplot(1,2,1)
distributionPlot(pc2,'histOpt',1,'colormap',1-gray(64),'showMM',0)
xlabel('glom. used for linear model')
ylabel('Unshuffled R^2')
set(gca,'XTickLabel',gNames)
xtickangle(30)
axis([0 length(gNames)+1 -0.9 0.9])
box off
set(gca,'FontSize',15)
subplot(1,2,2)
distributionPlot(pc2shuffled,'histOpt',1,'colormap',1-gray(64),'showMM',0)
xlabel('glom. used for linear model')
ylabel('Shuffled R^2')
set(gca,'XTickLabel',gNames)
xtickangle(30)
axis([0 length(gNames)+1 -0.9 0.9])
box off
set(gca,'FontSize',15)

%visualize predictive pcs for each odor
bestPCmeanloading=zeros(length(gNames),size(rawPC{1},1));
for j=1:length(gNames)
   %bestPCmeanloading(j,:)=mean(rawPC{j}(:,:,bestPCforodor(j)),2);
   
   % use single iteration of loadings rather than average
   [bestv besti]=max(myr2(bestPCforodor(j),:,j));
   bestPCmeanloading(j,:)=rawPC{j}(:,besti,bestPCforodor(j));
end
figure;
imagesc(bestPCmeanloading)
colormap(cm.egoalley)
set(gca,'YTick',[1:length(gNames)])
set(gca,'YTickLabel',string(gNames))
ylabel('glomerulus kept')
xtickangle(30)
set(gca,'XTick',[1:nodors])
set(gca,'XTickLabel',string(odornames))
ytickangle(30)
h=colorbar;
title(h,'PC loading')
set(gca,'FontSize',15)

figure
for i=1:nGloms
    subplot(2,3,i)
    plot(bestPCmeanloading(i,:),'k*','LineWidth',3)
    hold on
    plot(0:(nodors+1),zeros(nodors+2),'k--')
    axis([0 nodors+1 -0.6 0.6])
    set(gca,'XTick',[1:nodors])
    set(gca,'XTickLabel',string(odornames))
    xtickangle(40)
    ylabel('PC loading')
    title(gNames(i))
    box off
end
%% omit single odor
warning('off')

iters=100;
highestPCtouse=5;
yesZscore=0; % zscore data matrix before performing PCA (0=no, 1=yes)
normalizeResponse=0;

fracIn=0.4; % best results when fracIn is high, ~0.5, only using high confidence glomeruli
medianResponseOrTimeCourse=1; % 1 for median response only, 0 for full time course

nGloms=6;

testR=zeros(highestPCtouse,iters,nGloms);
testRshuffled=zeros(highestPCtouse,iters,nGloms);

for glomKeep=1:nGloms
    disp(['omitting glomerulus ' num2str(glomKeep)])
    if medianResponseOrTimeCourse
        responsesNoResponseRemoved=responsesGlomByOdor;
    else
        responsesNoResponseRemoved=responsesTimeCourseGlomByOdor;
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
    
    % remove one glomerulus
    responsesNoResponseRemoved(((glomKeep-1)*nodors+1):((glomKeep)*nodors),:)=0;
    
    % normalize responsesNoResponseRemoved
    if normalizeResponse
        for i=1:size(responsesNoResponseRemoved,2)
            responsesNoResponseRemoved(:,i)=responsesNoResponseRemoved(:,i)/nanmean(abs(responsesNoResponseRemoved(:,i)));
        end
    end
    
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
        
        
        % % fill nans with mean
        
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
        
        [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(traindata');
        co = SCORE;
        
        if jjj==1
            rawPC{glomKeep}=zeros(size(COEFF,1),iters,highestPCtouse);
        end

        for pcstouse=1:highestPCtouse
            currpc=COEFF(:,pcstouse);
            
            rawPC{glomKeep}(:,jjj,pcstouse)=currpc;
            
            % generate linear model based on current PC
            behaviorprediction=(traindata'*currpc);
            flyTruePref=zeros(1,length(trainflies));
            flyTruePrefShuffled=zeros(1,length(trainflies));
            ally=behaviorOcc';
            
            ally(traintoremove)=[];
            
            nactivity=zeros(length(trainflies),length(pcstouse));
            shuffledindtrain=randperm(length(trainflies)); % for comparing to shuffled null model
            for i=1:length(trainflies)
                flyTruePref(i)=mean(ally(flyindicestrain{i}));
                flyTruePrefShuffled(i)=mean(ally(flyindicestrain{shuffledindtrain(i)}));
                nactivity(i,:)=mean(behaviorprediction(flyindicestrain{i},:));
            end
            linmodel=fitlm(nactivity,flyTruePref);
            
            beta=linmodel.Coefficients.Estimate;
            
            PCContribution=currpc*beta(2:end);
            
            % apply predictor to testdata and evaluate how held out points fit in
            % model
            
            %behaviorprediction=testdata'*PCContribution;
            behaviorprediction=testdata'*currpc;
            testflyTruePref=zeros(1,length(holdoutflies));
            testflyTruePrefShuffled=zeros(1,length(holdoutflies));
            flyPredictedPref=zeros(1,length(holdoutflies));
            ally=behaviorOcc';
            
            ally(testtoremove)=[];
            
            testnactivity=zeros(length(holdoutflies),1);
            shuffledindtest=randperm(length(holdoutflies)); % for comparing to shuffled null model
            for i=1:length(holdoutflies)
                testflyTruePref(i)=mean(ally(flyindicestest{i}));
                testflyTruePrefShuffled(i)=mean(ally(flyindicestest{shuffledindtest(i)}));
                testnactivity(i,:)=mean(behaviorprediction(flyindicestest{i},:));
            end
            
            mytestprediction=predict(linmodel,testnactivity);
            
            myr=corrcoef(mytestprediction,testflyTruePref);
            myrshuffled=corrcoef(mytestprediction,testflyTruePref(randperm(length(testflyTruePref))));
            
            testR(pcstouse,jjj,glomKeep)=myr(1,2);
            testRshuffled(pcstouse,jjj,glomKeep)=myrshuffled(1,2);
            
        end
        
        
        if mod(jjj,10)==0
            disp(['iteration ' num2str(jjj)])
        end
    end
end


myr2=(testR.^2).*sign(testR);
myr2shuffled=(testRshuffled.^2).*sign(testRshuffled);

for i=1:highestPCtouse
    pc2=squeeze(myr2(i,:,:));
    pc2shuffled=squeeze(myr2shuffled(i,:,:));
    
    figure
    subplot(1,2,1)
    distributionPlot(pc2,'histOpt',1,'colormap',1-gray(64),'showMM',0)
    xlabel('glomerulus omitted')
    ylabel('Unshuffled R^2')
    set(gca,'XTickLabel',gNames)
    xtickangle(30)
    axis([0 length(gNames)+1 -0.8 0.8])
    box off
    set(gca,'FontSize',15)
    subplot(1,2,2)
    distributionPlot(pc2shuffled,'histOpt',1,'colormap',1-gray(64),'showMM',0)
    xlabel('glomerulus omitted')
    ylabel('Shuffled R^2')
    set(gca,'XTickLabel',gNames)
    xtickangle(30)
    axis([0 length(gNames)+1 -0.8 0.8])
    box off
    set(gca,'FontSize',15)
    title(['PC ' num2str(i)])
end



figure
imagesc(squeeze(mean(myr2,2)),[0 0.3])
set(gca,'XTick',[1:length(gNames)])
set(gca,'XTickLabel',gNames)
xtickangle(30)
ylabel('PC #')
xlabel('glomerulus omitted')
set(gca,'YTick',[1:highestPCtouse])
h=colorbar;
title(h,'R^2')
set(gca,'FontSize',15)

% find best pc for each odor
averageR2=squeeze(mean(myr2,2));
[bestR2 bestPCforodor]=max(averageR2,[],1);

for j=1:length(gNames)
pc2(:,j)=squeeze(myr2(bestPCforodor(j),:,j));
pc2shuffled(:,j)=squeeze(myr2shuffled(bestPCforodor(j),:,j));
end

figure
subplot(1,2,1)
distributionPlot(pc2,'histOpt',1,'colormap',1-gray(64),'showMM',0)
xlabel('glomerulus omitted')
ylabel('Unshuffled R^2')
set(gca,'XTickLabel',gNames)
xtickangle(30)
axis([0 length(gNames)+1 -0.8 0.8])
box off
set(gca,'FontSize',15)
subplot(1,2,2)
distributionPlot(pc2shuffled,'histOpt',1,'colormap',1-gray(64),'showMM',0)
xlabel('glomerulus omitted')
ylabel('Shuffled R^2')
set(gca,'XTickLabel',gNames)
xtickangle(30)
axis([0 length(gNames)+1 -0.8 0.8])
box off
set(gca,'FontSize',15)

% visualize predictive pcs for each odor
bestPCmeanloading=zeros(nGloms,size(rawPC{1},1));
for j=1:nGloms
  % bestPCmeanloading(j,:)=mean(rawPC{j}(:,:,bestPCforodor(j)),2);
   
  % use single iteration of loadings rather than average
   [bestv besti]=max(myr2(bestPCforodor(j),:,j));
   bestPCmeanloading(j,:)=rawPC{j}(:,besti,bestPCforodor(j));
end
figure;
imagesc(bestPCmeanloading)
colormap(cm.egoalley)
set(gca,'XTick',nodors/2+[1:nodors:(nodors*length(gNames))])
set(gca,'XTickLabel',string(gNames))
xtickangle(30)
set(gca,'YTick',[1:nGloms])
set(gca,'YTickLabel',string(gNames))
ylabel('glomerulus removed')
ytickangle(30)
h=colorbar;
title(h,'PC loading')
set(gca,'FontSize',15)

figure
for i=1:nGloms
    subplot(2,3,i)
    plot(bestPCmeanloading(i,:),'k*','LineWidth',3)
    hold on
    plot(0:(length(gNames)+1),zeros(length(gNames)+2),'k--')
    axis([0 nodors*length(gNames)+1 -0.5 .5])
    set(gca,'XTick',[1:nodors:(nodors*length(gNames))])
    set(gca,'XTickLabel',string(gNames))
    xtickangle(40)
    ylabel('PC loading')
    title(gNames(i))
    box off
end

