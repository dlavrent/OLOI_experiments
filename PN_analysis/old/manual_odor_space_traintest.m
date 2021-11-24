% look at manually labelled gh146 flies in odor space
% load all data
% only use 80% of data to a) determine principal components and b) fit a
% linear prediction to behavior
% use the remaining 20% of data, apply the predictors fit with training
% data, and measure how well the predictors predict behavior
clear all
close all
rng('default')

manualLabelHome='/Users/mattchurgin/Dropbox/flyimaging/analysis/PN_analysis/trainingData_through181210';

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

%% perform pca on responses

clear responsesNoResponseRemoved

fracIn=0.5;% best results when fracIn is high, ~0.5, only using high confidence glomeruli


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
    
    %  odorstoremove=[1 8]; % delete air and ethanol
    %  for i=1:length(odorstoremove)
    %      for j=1:nodors
    %      responsesNoResponseRemoved(((j-1)*nodors+odorstoremove(i)),:)=[];
    %      end
    %  end
else
    
end


% remove fully empty rows
%totallyempty=~any(isfinite(responsesNoResponseRemoved),2);
%responsesNoResponseRemoved(totallyempty,:)=[];

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

opt = statset('pca');
opt.Display='iter';
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(responsesNoResponseRemoved','Options',opt);


figure;
plot(cumsum(EXPLAINED),'o-','LineWidth',3)
ylabel('Variance Explained (%)')
xlabel('PC #')
box off
set(gca,'FontSize',15)

% calculate glom contribution to each pc
glomcontributionAbsolute=zeros(length(gNames),size(COEFF,2));
glomcontributionMean=zeros(length(gNames),size(COEFF,2));
odorcontributionAbsolute=zeros(nodors,size(COEFF,2));
odorcontributionMean=zeros(nodors,size(COEFF,2));
for i=1:length(gNames)
    if medianResponseOrTimeCourse
        temp=nodors;
    else
        temp=length(odortimesfortimecourse)*nodors;
    end
    
    currind=((i-1)*temp)+1;
    for j=1:size(COEFF,2)
        glomcontributionMean(i,j)=mean(COEFF(currind:(currind+temp-1),j));
        glomcontributionAbsolute(i,j)=mean(abs(COEFF(currind:(currind+temp-1),j)));
    end
end
if medianResponseOrTimeCourse
    for i=1:nodors
        for j=1:size(COEFF,2)
            odorcontributionMean(i,j)=mean(COEFF(i:nodors:end,j));
            odorcontributionAbsolute(i,j)=mean(abs(COEFF(i:nodors:end,j)));
        end
    end
end
figure;
imagesc(glomcontributionMean)
set(gca,'ytick',1:length(gNames),'yticklabel',string(gNames),'FontSize',10)
ytickangle(30)
xlabel('PC #')
set(gca,'FontSize',15)
%% find predictive power for each PC
warning('off')

iters=100;

trsize=30;

ally=behaviorOcc';
%ally(traintoremove)=[];

RON=zeros(size(co,2),iters);
PON=zeros(size(co,2),iters);
ROFF=zeros(size(co,2),iters);
POFF=zeros(size(co,2),iters);

for pcstouse=1:size(co,2)
    predictor=(co(:,pcstouse));
    predictor=predictor';
    
    bont=zeros(2,iters);
    bofft=zeros(2,iters);
    totaly=[];
    totalyp=[];
    totalprediction=[];
    totalpredictionp=[];
    
    for ii=1:iters
        
        temp=randperm(length(trainflies));
        trflies=temp(1:floor(length(temp)*trsize/100));
        validationflies=temp(ceil(length(temp)*trsize/100):end);
        
        trainset=[];
        for j=1:length(trflies)
            trainset=[trainset find(glombyodorfliestrain==trflies(j))];
            trainset=[trainset find(glombyodorfliestrain==(trflies(j)+leftRightSplitVal))];
        end
        
        validationset=[];
        for j=1:length(validationflies)
            validationset=[validationset find(glombyodorfliestrain==validationflies(j))];
            validationset=[validationset find(glombyodorfliestrain==(validationflies(j)+leftRightSplitVal))];
        end
        
        % Predict behavioral preference during odor
        % split data into train and test set
        bon=zeros(1,1);
        boff=zeros(1,1);
        
        % Predict behavioral preference during odor
        ytrain=ally(trainset);
        Xtr=predictor(:,trainset)';
        Xtrain=[ones(size(Xtr,1),1) Xtr];
        
        % run regression
        [bon,bint,r,rint,stats] = regress(ytrain,Xtrain);
        
        % Predict behavioral preference during preodor
        ytrainpre=behaviorpreOcc(trainset)';
        [boff,bintpre,rpre,rintpre,statspre] = regress(ytrainpre,Xtrain);
        
        % evaluate regression on held-out observations
        y=ally(validationset);
        ypre=behaviorpreOcc(validationset)';
        Xt=predictor(:,validationset)';
        X=[ones(size(Xt,1),1) Xt];
        
        totaly=[totaly y'];
        totalyp=[totalyp ypre'];
        totalprediction=[totalprediction (X*bon)'];
        totalpredictionp=[totalpredictionp (X*boff)'];
        
        [ron pon]=corrcoef(y,X*bon);
        [roff poff]=corrcoef(ypre,X*boff);
        RON(pcstouse,ii)=ron(1,2);
        PON(pcstouse,ii)=pon(1,2);
        ROFF(pcstouse,ii)=roff(1,2);
        POFF(pcstouse,ii)=poff(1,2);
        
        bont(:,ii)=bon;
        bofft(:,ii)=boff;
    end
    
    bonAverage=mean(bont,2);
    boffAverage=mean(bofft,2);
    if mod(pcstouse,10)==0
        disp(['calculated fit using ' num2str(pcstouse) ' PCs of ' num2str(size(co,2))])
    end
end

left_color=[0 0 0];
right_color=[0.5 0 0.9];

fig=figure;
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(nanmean(RON,2),'-','LineWidth',2)
ylabel('correlation')
hold on
yyaxis right
ylabel('p-value')
plot(nanmean(PON,2),'o-','LineWidth',2)
xlabel('# PCs')
set(gca,'FontSize',15)

fig=figure;
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(nanmean(ROFF,2),'-','LineWidth',2)
ylabel('correlation')
hold on
yyaxis right
ylabel('p-value')
plot(nanmean(POFF,2),'o-','LineWidth',2)
xlabel('# PCs')
set(gca,'FontSize',15)


%% use fitlm
pncolor=[0 0.8 0];
pcstouse=[2];

behaviorprediction=(SCORE(:,pcstouse));
flyTruePref=zeros(1,flyNum);
flyPredictedPref=zeros(1,flyNum);
ally=behaviorOcc';
linmodel=fitlm(behaviorprediction,ally);
myprediction=predict(linmodel,behaviorprediction);
figure
plot(myprediction,ally,'o','LineWidth',3)
xlabel('Predicted Preference')
ylabel('Measured Preference')
set(gca,'FontSize',15)
box off
linmodel

nactivity=zeros(flyNum,length(pcstouse));
for i=1:flyNum
    flyTruePref(i)=mean(ally(flyindices{i}));
    flyPredictedPref(i)=mean(mean(behaviorprediction(flyindices{i},:)));
    nactivity(i,:)=mean(behaviorprediction(flyindices{i},:));
end
linmodel=fitlm(nactivity,flyTruePref);
myprediction=predict(linmodel,nactivity);
figure
plot(myprediction,flyTruePref,'o','Color',pncolor,'LineWidth',3)
for i=1:flyNum
    hold on
    %text(myprediction(i)+0.01,flyTruePref(i),num2str(i),'FontSize',15)
end
xlabel('Predicted Preference')
ylabel('Measured Preference')
set(gca,'FontSize',15)
box off
linmodel

beta=linmodel.Coefficients.Estimate;

PCContribution=COEFF(:,pcstouse)*beta(2:end);
figure;
plot(PCContribution,'*','LineWidth',2,'MarkerSize',8)
hold on
plot(zeros(1,length(PCContribution(:,1))),'k--','LineWidth',3)
j=1;
for i=1:nodors:length(PCContribution)
    plot((i-0.5)*ones(1,5), linspace(min(PCContribution),max(PCContribution),5),'k--','LineWidth',2)
    %text(i+floor(nodors/3),min(PCContribution),num2str(glomsFound(j)),'FontSize',15)
    j=j+1;
end

set(gca,'xtick',(1:nodors:length(PCContribution))+floor(nodors/2),'xticklabel',string(gNames),'FontSize',10)
xtickangle(30)
ylabel('PC 4 loadings')
box off
set(gca,'FontSize',15)

currpc=COEFF(:,pcstouse);

% get odor valences
for i=1:nodors
    odorValence(i)=mean(PCContribution(i:nodors:end));
end


%% apply predictor to testdata

behaviorprediction=testdata'*PCContribution;
flyTruePref=zeros(1,length(holdoutflies));
flyPredictedPref=zeros(1,length(holdoutflies));
ally=behaviorOcc';
ally(testtoremove)=[];
linmodel=fitlm(behaviorprediction,ally);
myprediction=predict(linmodel,behaviorprediction);
figure
plot(myprediction,ally,'o','LineWidth',3)
xlabel('Predicted Preference')
ylabel('Measured Preference')
set(gca,'FontSize',15)
box off
linmodel

nactivity=zeros(length(holdoutflies),1);
for i=1:length(holdoutflies)
    flyTruePref(i)=mean(ally(flyindicestest{i}));
    flyPredictedPref(i)=mean(mean(behaviorprediction(flyindicestest{i},:)));
    nactivity(i,:)=mean(behaviorprediction(flyindicestest{i},:));
end
linmodel=fitlm(nactivity,flyTruePref);
myprediction=predict(linmodel,nactivity);
figure
plot(myprediction,flyTruePref,'o','LineWidth',3)
% for i=1:flyNum
%    hold on
%    text(myprediction(i)+0.01,flyTruePref(i),num2str(i),'FontSize',15)
% end
xlabel('Predicted Preference')
ylabel('Measured Preference')
set(gca,'FontSize',15)
box off
linmodel

beta=linmodel.Coefficients.Estimate;

%% use elastic net regression


flyTruePref=zeros(1,length(trainflies));
flyTruePreftest=zeros(1,length(holdoutflies));
ally=behaviorOcc';
ally(traintoremove)=[];
testbehavior=behaviorOcc';
testbehavior(testtoremove)=[];

nactivity=zeros(size(traindata,1),length(trainflies));
for i=1:length(trainflies)
    flyTruePref(i)=mean(ally(flyindicestrain{i}));
    
    nactivity(:,i)=mean(traindata(:,flyindicestrain{i}),2);
end

nactivitytest=zeros(size(testdata,1),length(holdoutflies));
for i=1:length(holdoutflies)
    flyTruePreftest(i)=mean(testbehavior(flyindicestest{i}));
    
    nactivitytest(:,i)=mean(testdata(:,flyindicestest{i}),2);
end


[b fitinfo]=lasso(nactivity',flyTruePref,'Alpha',0.7,'CV',30);
figure;
lassoPlot(b,fitinfo,'PlotType','CV');
legend('show') % Show legend

coef0=fitinfo.Intercept(1);
beta=b(:,1);

myprediction=coef0 + nactivity'*beta;
figure
plot(myprediction,flyTruePref,'o','LineWidth',3)
xlabel('Predicted Preference')
ylabel('Measured Preference')
set(gca,'FontSize',15)
box off

myprediction=coef0 + nactivitytest'*beta;
figure
plot(myprediction,flyTruePreftest,'o','LineWidth',3)
xlabel('Predicted Preference')
ylabel('Measured Preference')
set(gca,'FontSize',15)
box off

figure;
plot(b(:,1),'*','LineWidth',3,'MarkerSize',10)
hold on
plot(zeros(1,length(b(:,1))),'k--','LineWidth',3)
j=1;
for i=1:nodors:length(b(:,1))
    plot((i-0.5)*ones(1,5), linspace(min(b(:,1)),max(b(:,1)),5),'k--','LineWidth',2)
    text(i+floor(nodors/3),min(b(:,1)),num2str(glomsFound(j)),'FontSize',15)
    j=j+1;
end

set(gca,'xtick',(1:nodors:length(b(:,1)))+floor(nodors/2),'xticklabel',string(gNames),'FontSize',10)
xtickangle(30)
ylabel('Weight')
box off
set(gca,'FontSize',15)

%% run train/test many times
warning('off')

iters=250;
highestPCtouse=80;
testR=zeros(highestPCtouse,iters);
averagePredictor=cell(highestPCtouse,1);
vExplained=cell(1,highestPCtouse);

fracIn=0.05;
medianResponseOrTimeCourse=1; % 1 for median response only, 0 for full time course

% fill missing values
responsesMissingOdorsFixed=responsesGlomByOdor;
odorstoremove=[2 6 9 10 12];
blankThresh=[0.3 100 100 0.3 0.3]; % blank odor if 75% prctile response is below this (for odors that sporadically worked)

% mean fill time points when odors failed
for i=1:length(odorstoremove)
    flyodorstoblank=find(prctile(responsesMissingOdorsFixed(odorstoremove(i):nodors:end,:),75)<=blankThresh(i));
    for j=1:(size(responsesMissingOdorsFixed,1)/nodors)
        temp=find(isfinite(responsesMissingOdorsFixed(((j-1)*nodors+odorstoremove(i)),[1:79])));
        temp2=find(isfinite(responsesMissingOdorsFixed(((j-1)*nodors+odorstoremove(i)),[115:end])));
        tofill=intersect(flyodorstoblank,temp);
        tofill2=intersect(flyodorstoblank,temp2)+114;
        if isfinite(nanmean(responsesMissingOdorsFixed(((j-1)*nodors+odorstoremove(i)),80:114)))
            responsesMissingOdorsFixed(((j-1)*nodors+odorstoremove(i)),tofill)=nanmean(responsesMissingOdorsFixed(((j-1)*nodors+odorstoremove(i)),80:114));
            responsesMissingOdorsFixed(((j-1)*nodors+odorstoremove(i)),tofill2)=nanmean(responsesMissingOdorsFixed(((j-1)*nodors+odorstoremove(i)),80:114));
        end
    end
end

% remove MCH for flies on 181108
for j=1:(size(responsesMissingOdorsFixed,1)/nodors)
    tofill=find(isfinite(responsesMissingOdorsFixed(((j-1)*nodors+11),80:92)))+79;
    if isfinite(nanmean(responsesMissingOdorsFixed(((j-1)*nodors+11),[1:79 93:end])))
        responsesMissingOdorsFixed(((j-1)*nodors+11),tofill)=nanmean(responsesMissingOdorsFixed(((j-1)*nodors+11),[1:79 93:end]));
    end
end


for jjj=1:iters
    
    clear responsesNoResponseRemoved
    
    responsesNoResponseRemoved=responsesMissingOdorsFixed;
    
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
    
    
    % Split data into randomly assigned train/test sets
    testsize=20;
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
    
    traindata=responsesNoResponseRemoved;
    traindata(:,traintoremove)=[];
    
    % % fill nans with mean
    for i=1:size(traindata,1)
        for j=1:size(traindata,2)
            if isnan(traindata(i,j))
                traindata(i,j)=nanmean(traindata(i,:));
            end
        end
    end
    
    for i=1:size(responsesNoResponseRemoved,1)
        for j=1:size(responsesNoResponseRemoved,2)
            if isnan(responsesNoResponseRemoved(i,j))
                responsesNoResponseRemoved(i,j)=nanmean(responsesNoResponseRemoved(i,:));
            end
        end
    end
    testdata=responsesNoResponseRemoved;
    testdata(:,testtoremove)=[];
    
    %     % z-score train and test data
    %
    %     ztrain=traindata';
    %     ztrain=(ztrain-mean(ztrain))./(std(ztrain));
    %     ztrain=ztrain';
    %
    %     ztest=testdata';
    %     ztest=(ztest-mean(ztest))./(std(ztest));
    %     ztest=ztest';
    %
    %     traindata=ztrain;
    %     testdata=ztest;
    
    %     opt = statset('pca');
    %     opt.Display='iter';
    %[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(traindata','Options',opt);
    [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(traindata');
    
    vExplained{jjj}=EXPLAINED;
    co = SCORE;
    
    for pcstouse=1:highestPCtouse
        currpc=COEFF(:,pcstouse);
        
        % generate linear model based on current PC
        behaviorprediction=(traindata'*currpc);
        flyTruePref=zeros(1,length(trainflies));
        ally=behaviorOcc';
        ally(traintoremove)=[];
        
        nactivity=zeros(length(trainflies),length(pcstouse));
        for i=1:length(trainflies)
            flyTruePref(i)=mean(ally(flyindicestrain{i}));
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
        flyPredictedPref=zeros(1,length(holdoutflies));
        ally=behaviorOcc';
        ally(testtoremove)=[];
        
        testnactivity=zeros(length(holdoutflies),1);
        for i=1:length(holdoutflies)
            testflyTruePref(i)=mean(ally(flyindicestest{i}));
            testnactivity(i,:)=mean(behaviorprediction(flyindicestest{i},:));
        end
        
        mytestprediction=predict(linmodel,testnactivity);
        
        myr=corrcoef(mytestprediction,testflyTruePref);
        
        if jjj==1
            averagePredictor{pcstouse}=PCContribution';
        else
            averagePredictor{pcstouse}=averagePredictor{pcstouse}+PCContribution';
        end
        
        testR(pcstouse,jjj)=myr(1,2);
        
    end
    if mod(jjj,1)==0
        disp(['iteration ' num2str(jjj)])
    end
    if mod(jjj,20)==0
        boxplot(testR(:,1:jjj)')
        xlabel('PC #')
        ylabel('Test Data r value')
        box off
        set(gca,'FontSize',15)
        drawnow
    end
end

figure
boxplot(testR')
xlabel('PC used for linear model')
ylabel('Correlation between predicted vs. true preference (test data)')
box off
set(gca,'FontSize',15)

figure;
plot(averagePredictor{4},'*','LineWidth',3,'MarkerSize',10)
hold on
plot(averagePredictor{5},'*','LineWidth',3,'MarkerSize',10)

plot(zeros(1,length(averagePredictor{4})),'k--','LineWidth',3)
j=1;
for i=1:nodors:length(averagePredictor{4})
    plot((i-0.5)*ones(1,100), linspace(min(averagePredictor{4}),max(averagePredictor{4})),'k--','LineWidth',2)
    text(i+floor(nodors/3),min(averagePredictor{4}),num2str(glomsFound(j)),'FontSize',15)
    j=j+1;
end
set(gca,'xtick',(1:nodors:length(averagePredictor{4}))+floor(nodors/2),'xticklabel',string(gNames),'FontSize',10)
xtickangle(30)
ylabel('Loadings')
box off
legend('PC 4','PC 5')
legend boxoff
set(gca,'FontSize',15)

% plot variance explained for first 80 PCs, averaged over all replicates
vExplainedCropped=zeros(80,1);
for i=1:iters
    vExplainedCropped=vExplainedCropped+vExplained{i}(1:80);
end
vExplainedCropped=vExplainedCropped/iters;

figure;
subplot(2,1,1)
plot(vExplainedCropped,'LineWidth',2)
xlabel('PC #')
ylabel('Variance Explained')
set(gca,'FontSize',15)
box off
subplot(2,1,2)
plot(log10(vExplainedCropped),'LineWidth',2)
xlabel('PC #')
ylabel('log10(Variance Explained)')
set(gca,'FontSize',15)
box off

figure
npairs=64;
for i=1:npairs
    subplot(sqrt(npairs),sqrt(npairs),i)
    plot(averagePredictor{i},averagePredictor{i+1},'*','LineWidth',2)
    box off
    set(gca,'FontSize',15)
end

figure
npairs=64;
for i=1:npairs
    subplot(sqrt(npairs),sqrt(npairs),i)
    plot(averagePredictor{i},'*','LineWidth',1)
    box off
    set(gca,'FontSize',10)
end


%% run train/test many times and save replicate data
warning('off')

iters=100;
highestPCtouse=10;
yesZscore=0; % zscore data matrix before performing PCA (0=no, 1=yes)


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

% remove air and ethanol
responsesNoResponseRemoved(1:nodors:end,:)=0;
responsesNoResponseRemoved(8:nodors:end,:)=0;

testR=zeros(highestPCtouse,iters);
testNRSS=zeros(highestPCtouse,iters); % normalized residual sum of squares
averagePredictor=cell(highestPCtouse,1);
bestPredictorLoading=cell(iters,1); % loadings for the best predictor
bestPredictor=zeros(iters,1); % PC which is the best predictor
bestPredictorR=zeros(iters,1);
pcPredictorRank=zeros(iters,highestPCtouse);
vExplained=cell(1,highestPCtouse);
relativeOctMchactivation=zeros(iters,highestPCtouse);
rawRelativeOctactivation=zeros(iters,highestPCtouse,length(gNames));
rawRelativeMchactivation=zeros(iters,highestPCtouse,length(gNames));
rawPC=cell(1,highestPCtouse);

bestPCOctMchactivation=zeros(1,iters);
bestPCOctMchcorrelation=zeros(1,iters);
octmchloadingcorr=zeros(iters,highestPCtouse); % correlation between oct and mch loading for each pc

for jjj=1:iters
    
    % Split data into randomly assigned train/test sets
    testsize=20;
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
    
    vExplained{jjj}=EXPLAINED;
    
    co = SCORE;
    
    
    rawActivation{jjj}=zeros(size(COEFF,1),highestPCtouse);
    
    
    for pcstouse=1:highestPCtouse
        currpc=COEFF(:,pcstouse);
        
        if jjj==1
            rawPC{pcstouse}=zeros(size(COEFF,1),iters);
        end
        rawPC{pcstouse}(:,jjj)=currpc;
        
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
        
        if jjj==1
            %averagePredictor{pcstouse}=PCContribution';
            averagePredictor{pcstouse}=currpc';
        else
            %averagePredictor{pcstouse}=averagePredictor{pcstouse}+PCContribution';
            averagePredictor{pcstouse}=averagePredictor{pcstouse}+currpc';
        end
        
        testR(pcstouse,jjj)=myr(1,2);
        testNRSS(pcstouse,jjj)=nanmean((((testflyTruePref-mytestprediction')).^2));
        
        relativeOctMchactivation(jjj,pcstouse)=mean((currpc(2:13:end)-currpc(11:13:end))./(1+(abs(currpc(2:13:end))+abs(currpc(11:13:end)))));
        rawRelativeOctactivation(jjj,pcstouse,:)=currpc(2:13:end);
        rawRelativeMchactivation(jjj,pcstouse,:)=currpc(11:13:end);
        
        [octmchcorr asd]=corr(currpc(2:13:end),currpc(11:13:end));
        octmchloadingcorr(jjj,pcstouse)=octmchcorr;
        
    end
    
    % find best predictor
    [val ind]=max(testR(:,jjj));
    bestPredictor(jjj)=ind;
    bestPredictorLoading{jjj}=COEFF(:,ind);
    bestPredictorR(jjj)=val;
    [vals inds]=sort(testR(:,jjj),'descend');
    pcPredictorRank(jjj,:)=inds;
    
    bestpc=COEFF(:,ind);
    bestPCOctMchactivation(jjj)=mean((bestpc(2:13:end)-bestpc(11:13:end))./(1+(abs(bestpc(2:13:end))+abs(bestpc(11:13:end)))));
    bestPCOctMchcorrelation(jjj)=octmchloadingcorr(jjj,ind);
    
    if mod(jjj,1)==0
        disp(['iteration ' num2str(jjj)])
    end
    %     if mod(jjj,20)==0
    %         boxplot(testR(:,1:jjj)')
    %         xlabel('PC #')
    %         ylabel('Test Data r value')
    %         box off
    %         set(gca,'FontSize',15)
    %         drawnow
    %     end
    %catch
    %end
end


for i=1:highestPCtouse
    averagePredictor{i}=averagePredictor{i}/iters;
end


%rawOctMchdifference=(rawRelativeOctactivation-rawRelativeMchactivation)./(1+abs(rawRelativeOctactivation)+abs(rawRelativeMchactivation));
rawOctMchdifference=(rawRelativeOctactivation-rawRelativeMchactivation);
totalOctMchdifference=zeros(highestPCtouse,iters*length(gNames));
for j=1:highestPCtouse
    totalOctMchdifference(j,:)=reshape(rawOctMchdifference(:,j,:),iters*length(gNames),1);
    
    % plot histogram of oct-mch loadings
    %     nbins=linspace(-0.55,0.55,50);
    %     [ay ax]=hist(totalOctMchdifference(j,:),nbins);
    %     plot(ax,log10(ay),'k','LineWidth',2)
    %     axis([-0.5 0.5 0 4])
    %     pause
end

meanAbsoluteOctMchdifference=zeros(iters,highestPCtouse);
for j=1:iters
    meanAbsoluteOctMchdifference(j,:)=mean(rawOctMchdifference(j,:,:),3);
end

figure
distributionPlot(meanAbsoluteOctMchdifference(:,:),'histOpt',1,'colormap',1-gray(64),'showMM',0);
xlabel('PC #')
ylabel('Average OCT-MCH Loading')
set(gca,'FontSize',15)

figure;
boxplot(octmchloadingcorr)
box off
xlabel('PC #')
ylabel('Correlation (OCT vs. MCH loading)')
set(gca,'FontSize',15)

figure
%subplot(2,1,1)
boxplot(testR')
xlabel('PC used for linear model')
ylabel('Correlation (true vs. predicted preference)')
box off
set(gca,'FontSize',15)
% subplot(2,1,2)
% boxplot(testNRSS')
% xlabel('PC used for linear model')
% ylabel('Normalized residual sum of squares')
% box off
% set(gca,'FontSize',15)


figure
plot(averagePredictor{4},'*','LineWidth',2,'MarkerSize',8)
hold on
plot(zeros(1,size(COEFF,1)),'k--','LineWidth',3)
j=1;
for i=1:nodors:length(averagePredictor{4})
    plot((i-0.5)*ones(1,100), linspace(min(averagePredictor{4}),max(averagePredictor{4})),'k--','LineWidth',2)
    %text(i+floor(nodors/3),min(averagePredictor{4}),num2str(glomsFound(j)),'FontSize',15)
    j=j+1;
end
set(gca,'xtick',(1:nodors:length(averagePredictor{4}))+floor(nodors/2),'xticklabel',string(gNames),'FontSize',10)
xtickangle(30)
ylabel('PC 4 loadings')
box off
set(gca,'FontSize',15)


figure
subplot(2,2,1)
boxplot(testR')
xlabel('PC used for linear model')
ylabel('Correlation between predicted vs. true preference (test data)')
box off
set(gca,'FontSize',15)

subplot(2,2,2)
hist(bestPredictor)
xlabel('Best predictor PC #')

% boxplot(pcPredictorRank)
% xlabel('PC #')
% ylabel('PC Predictor Rank')

subplot(2,2,3)
hist(bestPredictorR)
xlabel('Best predictor''s r-value')

bploadings=zeros(1,size(COEFF,1));
allbpl=zeros(iters,size(COEFF,1));
for j=1:iters
    bploadings=bploadings+bestPredictorLoading{j}';
    allbpl(j,:)=bestPredictorLoading{j}';
end
bploadings=bploadings/iters;


subplot(2,2,4)
plot(bploadings,'*','LineWidth',3,'MarkerSize',10)
hold on
plot(zeros(1,size(COEFF,1)),'k--','LineWidth',3)
j=1;
for i=1:nodors:length(averagePredictor{4})
    plot((i-0.5)*ones(1,100), linspace(min(bploadings),max(bploadings)),'k--','LineWidth',2)
    text(i+floor(nodors/3),min(bploadings),num2str(glomsFound(j)),'FontSize',15)
    j=j+1;
end
set(gca,'xtick',(1:nodors:length(bploadings))+floor(nodors/2),'xticklabel',string(gNames),'FontSize',10)
xtickangle(30)
ylabel('Best predictor''s loadings')
box off
set(gca,'FontSize',15)

% figure
% boxplot(relativeOctMchactivation)
% xlabel('PC #')
% ylabel('Abs(Oct-MCH) Loading')
% set(gca,'FontSize',15)
%
% figure
% hist(relativeOctMchactivation(:,4),50)
% xlabel('PC #')
% ylabel('Abs(Oct-MCH) Loading')
% set(gca,'FontSize',15)
%
% figure
% distributionPlot(relativeOctMchactivation,'histOpt',1,'colormap',1-gray(64));
% xlabel('PC #')
% ylabel('Average OCT-MCH Loading')
% set(gca,'FontSize',15)

% show best predictor loadings for each iteration
bpls=zeros(length(bestPredictorLoading{1}),iters);
oa=zeros(1,iters);
for j=1:iters
    bpls(:,j)=bestPredictorLoading{j}';
    oa(j)=corr(bestPredictorLoading{j}(2:13:end),bestPredictorLoading{j}(11:13:end));
end
figure;
imagesc(bpls)
xlabel('iteration')
ylabel('dimension')

figure; hist(oa)

ex=zeros(1,50);
for j=1:iters
    ex=ex+vExplained{j}(1:50)';
end
ex=ex/iters;
figure
plot(log10(ex),'LineWidth',3)
xlabel('PC #')
ylabel('log10(Variance Explained)')
set(gca,'FontSize',20)

% run kmeans on allbpl (how many "best pcs" are there)
[IDX, bestpredictorC, SUMD]= kmeans(allbpl, 20);
figure;
imagesc(bestpredictorC)
xlabel('Dimensions')
ylabel('Predictor Class (from kmeans)')
title('Best PC predictor loadings (all)')
set(gca,'FontSize',15)

% make a "best predictor matrix"
bestpredictors=bestpredictorC(1,:);
for i=2:size(bestpredictorC,1)
    for j=1:size(bestpredictors,1)
        temp=corrcoef(bestpredictorC(i,:),bestpredictors(j,:));
        temp2=temp(1,2);
        if abs(temp2)>0.4
            if temp2>0
                bestpredictors(j,:)=(bestpredictors(j,:)+bestpredictorC(i,:))/2;
            else
                bestpredictors(j,:)=(bestpredictors(j,:)-bestpredictorC(i,:))/2;
            end
            break
        elseif j==size(bestpredictors,1)
            bestpredictors=[bestpredictors; bestpredictorC(i,:)];
        end
        
    end
end

figure;
imagesc(bestpredictors)
xlabel('Dimensions')
ylabel('Predictor Class (from kmeans)')
title('Best PC predictor loadings (clustered)')
set(gca,'FontSize',15)


% run kmeans on each rawPC sweep to find possibilities
maxK=20;
for i=1:highestPCtouse
    [IDX, C, SUMD]= kmeans(rawPC{i}', maxK);
    kidx{i}=IDX;
    kc{i}=C;
end
figure
for i=1:highestPCtouse
    subplot(2,5,i)
    boxplot(testR(i,:),kidx{i})
    ylabel('Measured-Predicted Fit Correlation')
    title(['PC ' num2str(i)])
end

figure
for i=1:highestPCtouse
    subplot(2,5,i)
    
    boxplot(octmchloadingcorr(:,i),kidx{i})
    ylabel('OCT-MCH correlation')
    title(['PC ' num2str(i)])
end

% go through all PCs and select each one which led to a behavior prediction
predcorr=0.2; % select PCs with behavior prediction correlation higher than this number
bestOfPc=[];
for i=1:highestPCtouse
    for j=1:maxK
        if median(testR(i,kidx{i}==j))>predcorr
            if bestOfPc
                bestOfPc=kc{i}(j,:);
            else
                bestOfPc=[bestOfPc; kc{i}(j,:)];
            end
            disp(['saving PC# ' num2str(i) ' cluster ' num2str(j)])
        end
    end
end

% make a "best of predictor matrix"
bestpredictorsAll=bestOfPc(1,:);
for i=2:size(bestOfPc,1)
    for j=1:size(bestpredictorsAll,1)
        temp=corrcoef(bestOfPc(i,:),bestpredictorsAll(j,:));
        temp2=temp(1,2);
        if abs(temp2)>0.6
            if temp2>0
                bestpredictorsAll(j,:)=(bestpredictorsAll(j,:)+bestOfPc(i,:))/2;
            else
                bestpredictorsAll(j,:)=(bestpredictorsAll(j,:)-bestOfPc(i,:))/2;
            end
            break
        elseif j==size(bestpredictorsAll,1)
            bestpredictorsAll=[bestpredictorsAll; bestOfPc(i,:)];
        end
        
    end
end

figure;
imagesc(bestpredictorsAll)
xlabel('Dimensions')
ylabel('Predictor Class (from kmeans)')
title('Best PC predictor loadings (clustered)')
set(gca,'FontSize',15)
