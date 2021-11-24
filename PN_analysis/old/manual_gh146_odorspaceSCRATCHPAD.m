% look at manually labelled gh146 flies in odor space
clear all
close all

manualLabelHome='/Users/mattchurgin/Dropbox/flyimaging/analysis/PN_analysis/manuallyLabelled';

publishedOdorPath='/Users/mattchurgin/Desktop/dblab/mattFunctions/odorpanelResponsesDoOR/odorPanel_12/odorPanel_12_DoORData.mat';
load(publishedOdorPath);

manualLabelledFolders=dir(manualLabelHome);
manualLabelledFolders=manualLabelledFolders(3:end);

labels=cell(1,length(manualLabelledFolders));
glomsannotated=zeros(1,length(manualLabelledFolders));
totalPutativeGloms=zeros(1,length(manualLabelledFolders));

nodors=13;
odortimes=[6:10]; % hard-coded specific time interval for summarizing odor response
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
    currDate=currname(1:6);
    underscores=strfind(currname,'_');
    currFlyNum=currname((underscores(1)+4):(underscores(2)-1));
    
    if ~strcmp(currDate,oldDate) || ~strcmp(currFlyNum,oldFlyNum)
        flyNum=flyNum+1;
    end
    
    currFiles=dir([manualLabelHome '/' manualLabelledFolders(i).name]);
    currFiles=currFiles(3:end);
    for j=1:length(currFiles)
        % load manual label file, response data, and behavior data
        load([manualLabelHome '/' manualLabelledFolders(i).name '/' currFiles(j).name])
    end
    
    %behaviorOcc=[behaviorOcc occ_zeroed];
    %behaviorpreOcc=[behaviorpreOcc preocc_zeroed];
    behaviorOcc=[behaviorOcc occ];
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

%% perform pca on responses

clear responsesNoResponseRemoved

fracIn=0.01;


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
    %odorstoremove=[];
    
    % hard coded all flies before fly 23.  1:81 !!! Hard coded!!
    for i=1:length(odorstoremove)
        flyodorstoblank=find(prctile(responsesNoResponseRemoved(odorstoremove(i):nodors:end,:),75)<=blankThresh(i));
        for j=1:(size(responsesNoResponseRemoved,1)/nodors)
            temp=find(isfinite(responsesNoResponseRemoved(((j-1)*nodors+odorstoremove(i)),1:81)));
            tofill=intersect(flyodorstoblank,temp);
            if isfinite(nanmean(responsesNoResponseRemoved(((j-1)*nodors+odorstoremove(i)),82:end)))
                responsesNoResponseRemoved(((j-1)*nodors+odorstoremove(i)),tofill)=nanmean(responsesNoResponseRemoved(((j-1)*nodors+odorstoremove(i)),82:end));
            end
        end
    end
    
    % remove MCH for flies on 181108
    for j=1:(size(responsesNoResponseRemoved,1)/nodors)
        tofill=find(isfinite(responsesNoResponseRemoved(((j-1)*nodors+11),82:94)))+81;
        if isfinite(nanmean(responsesNoResponseRemoved(((j-1)*nodors+11),[1:81 95:end])))
        responsesNoResponseRemoved(((j-1)*nodors+11),tofill)=nanmean(responsesNoResponseRemoved(((j-1)*nodors+11),[1:81 95:end]));
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
%%
varianceToKeep=50; % percent of variance to keep


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
    
    
    if length(lcurr)>0 && length(rcurr)>0
        plot(co(lcurr(1:size(lcurr,2)),1),co(lcurr(1:size(lcurr,2)),2),'Color',mycmap(j,:),'Marker','*','LineWidth',lsize,'MarkerSize',msize);
        
        
        plot(co(rcurr(1:size(rcurr,2)),1),co(rcurr(1:size(rcurr,2)),2),'Color',mycmap(j,:),'Marker','o','LineWidth',lsize,'MarkerSize',msize);
        
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
    
    
    if length(lcurr)>0 && length(rcurr)>0
        h3=plot(mean(co(lcurr(1:size(lcurr,2)),1)),mean(co(lcurr(1:size(lcurr,2)),2)),'Color',mycmap(j,:),'Marker','*','LineWidth',lsize,'MarkerSize',msize);
        %text(mean(co(lcurr(1:size(lcurr,2)),1)),mean(co(lcurr(1:size(lcurr,2)),2)),[num2str(j)],'FontSize',15)
        h4=plot(mean(co(rcurr(1:size(rcurr,2)),1)),mean(co(rcurr(1:size(rcurr,2)),2)),'Color',mycmap(j,:),'Marker','o','LineWidth',lsize,'MarkerSize',msize);
        
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

%% use a single PC to directly predict behavior
behaviorprediction=SCORE(:,1)';
ally=behaviorOcc';

%collapse all single fly points into one data point
flyTruePref=zeros(1,flyNum);
flyPredictedPref=zeros(1,flyNum);
for i=1:flyNum
   flyTruePref(i)=mean(ally(flyindices{i}));
   flyPredictedPref(i)=mean(behaviorprediction(1,flyindices{i}));
end

[rPredicted pPredicted]=corrcoef(flyTruePref,flyPredictedPref);

figure
plot(flyPredictedPref,flyTruePref,'*','LineWidth',3)
xlabel('PC Score 1')
ylabel('Measured Preference')
box off
text(0.8*max(flyPredictedPref),1.1*min(flyTruePref),['r = ' num2str(rPredicted(1,2))],'FontSize',20)
set(gca,'FontSize',20)


% get prediction for left and right lobes separately
flyPredictedPrefL=zeros(1,flyNum);
flyPredictedPrefR=zeros(1,flyNum);
for i=1:flyNum
   flyPredictedPrefL(i)=mean(behaviorprediction(1,flyindicesL{i}));
   flyPredictedPrefR(i)=mean(behaviorprediction(1,flyindicesR{i}));
end
[rPredictedL pPredictedL]=corrcoef(flyTruePref(isfinite(flyPredictedPrefL)),flyPredictedPrefL(isfinite(flyPredictedPrefL)));
[rPredictedR pPredictedR]=corrcoef(flyTruePref(isfinite(flyPredictedPrefR)),flyPredictedPrefR(isfinite(flyPredictedPrefR)));
figure
plot(flyPredictedPrefL,flyTruePref,'*','LineWidth',3)
hold on
plot(flyPredictedPrefR,flyTruePref,'o','LineWidth',3)
xlabel('PC Score 1')
ylabel('Measured Preference')
legend('Left Lobe','Right Lobe')
legend boxoff
box off
text(0.8*max(flyPredictedPrefL),1.1*min(flyTruePref),['r = ' num2str(rPredictedL(1,2))],'FontSize',20)
text(0.8*max(flyPredictedPrefL),0.9*min(flyTruePref),['r = ' num2str(rPredictedR(1,2))],'FontSize',20)
set(gca,'FontSize',20)

figure
for i=1:length(flyTruePref)
plot(flyPredictedPrefL(i),flyTruePref(i),'*','LineWidth',3)
hold on
plot(flyPredictedPrefR(i),flyTruePref(i),'o','LineWidth',3)
text(flyPredictedPrefL(i),flyTruePref(i),[num2str(i)],'FontSize',15)
text(flyPredictedPrefR(i),flyTruePref(i),[num2str(i)],'FontSize',15)
end
xlabel('PC Score 1')
ylabel('Measured Preference')
legend('Left Lobe','Right Lobe')
legend boxoff
box off
text(0.8*max(flyPredictedPrefL),1.1*min(flyTruePref),['r = ' num2str(rPredictedL(1,2))],'FontSize',20)
text(0.8*max(flyPredictedPrefL),0.9*min(flyTruePref),['r = ' num2str(rPredictedR(1,2))],'FontSize',20)
set(gca,'FontSize',20)


% % fit linear model with averages
 linModel=fitlm(flyPredictedPref,flyTruePref);
 neuralPrediction=predict(linModel,flyPredictedPref');

figure
plot(neuralPrediction,flyTruePref,'*','LineWidth',3)
ylabel('True Preference')
xlabel('Predicted Preference')
box off
set(gca,'FontSize',20)

%% use fitlm

pcstouse=[3 4 23];

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

PCContribution=COEFF(:,pcstouse)*beta(2:end);
figure;
plot(PCContribution,'*','LineWidth',3,'MarkerSize',10)
hold on
plot(zeros(1,length(PCContribution(:,1))),'k--','LineWidth',3)
j=1;
for i=1:nodors:length(PCContribution)
   plot((i-0.5)*ones(1,5), linspace(min(PCContribution),max(PCContribution),5),'k--','LineWidth',2)
   text(i+floor(nodors/3),min(PCContribution),num2str(glomsFound(j)),'FontSize',15)
   j=j+1;
end

set(gca,'xtick',(1:nodors:length(PCContribution))+floor(nodors/2),'xticklabel',string(gNames),'FontSize',10)
xtickangle(30)
ylabel('Weight')
box off
set(gca,'FontSize',15)

%% use fitlm
% use only particular gloms


currpc=PCContribution;

% blank any coefficient with small contribution
%currpc([27:143 183:364 378:end])=0;
%currpc([27:143 183:246 274:364 378:end])=0;

behaviorprediction=responsesNoResponseRemoved'*currpc;
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

%% train behavior prediction from individual flies neuronal response in coding space
varianceToKeep=30; % percent of variance to keep
totalVarianceExplained=cumsum(EXPLAINED);
pcsWithinVariance=find(totalVarianceExplained<varianceToKeep);
pcstouse=pcsWithinVariance(end);
disp(['Using ' num2str(pcstouse) ' PCs'])

predictor=(co(:,1:pcstouse)-mean(co(:,1:pcstouse),1))./std(co(:,1:pcstouse));
predictor=predictor';

totaly=[];
totalyp=[];
totalprediction=[];
totalpredictionp=[];

iters=1000;

trainsize=50;

RON=zeros(1,iters);
PON=zeros(1,iters);
ROFF=zeros(1,iters);
POFF=zeros(1,iters);
bont=zeros(pcstouse+1,iters);
bofft=zeros(pcstouse+1,iters);

for ii=1:iters
    
    temp=randperm(flyNum-1);
    trainflies=temp(1:floor(length(temp)*trainsize/100));
    validationflies=temp(ceil(length(temp)*trainsize/100):end);
    
    trainset=[];
    for j=1:length(trainflies)
        trainset=[trainset find(glombyodorflies==trainflies(j))];
        trainset=[trainset find(glombyodorflies==(trainflies(j)+leftRightSplitVal))];
    end
    
    validationset=[];
    for j=1:length(validationflies)
        validationset=[validationset find(glombyodorflies==validationflies(j))];
        validationset=[validationset find(glombyodorflies==(validationflies(j)+leftRightSplitVal))];
    end
    
    % Predict behavioral preference during odor
    % split data into train and test set
    bon=zeros(pcstouse+1,1);
    boff=zeros(pcstouse+1,1);
    
    % Predict behavioral preference during odor
    ytrain=behaviorOcc(trainset)';
    Xtr=predictor(:,trainset)';
    Xtrain=[ones(size(Xtr,1),1) Xtr];
    
    % run regression
    [bon,bint,r,rint,stats] = regress(ytrain,Xtrain);
    
    % Predict behavioral preference during preodor
    ytrainpre=behaviorpreOcc(trainset)';
    [boff,bintpre,rpre,rintpre,statspre] = regress(ytrainpre,Xtrain);
    
    % evaluate regression on held-out observations
    y=behaviorOcc(validationset)';
    ypre=behaviorpreOcc(validationset)';
    Xt=predictor(:,validationset)';
    X=[ones(size(Xt,1),1) Xt];
    
    totaly=[totaly y'];
    totalyp=[totalyp ypre'];
    totalprediction=[totalprediction (X*bon)'];
    totalpredictionp=[totalpredictionp (X*boff)'];
    
    [ron pon]=corrcoef(y,X*bon);
    [roff poff]=corrcoef(ypre,X*boff);
    RON(ii)=ron(1,2);
    PON(ii)=pon(1,2);
    ROFF(ii)=roff(1,2);
    POFF(ii)=poff(1,2);
    
    bont(:,ii)=bon;
    bofft(:,ii)=boff;
    
    if mod(ii,1000)==0
        disp(['iteration ' num2str(ii) ' of ' num2str(iters)])
    end
end

bonAverage=mean(bont,2);
boffAverage=mean(bofft,2);

% show prediction using average coefficients
% not a fully unbiased measure, since this data was used to fit the
% coefficients
% but it gives an idea
ally=behaviorOcc';
allypre=behaviorpreOcc';
Xt=predictor';
X=[ones(size(Xt,1),1) Xt];

%collapse all single fly points into one data point
flyTruePref=zeros(1,flyNum);
flyPredictedPref=zeros(1,flyNum);
flyTruePrefPre=zeros(1,flyNum);
flyPredictedPrefPre=zeros(1,flyNum);
for i=1:flyNum
   flyTruePref(i)=mean(ally(flyindices{i}));
   flyTruePrefPre(i)=mean(allypre(flyindices{i}));
   flyPredictedPref(i)=mean(X(flyindices{i},:)*bonAverage);
   flyPredictedPrefPre(i)=mean(X(flyindices{i},:)*boffAverage);
end

[rPredicted pPredicted]=corrcoef(flyTruePref,flyPredictedPref);
[rPredictedPre pPredictedPre]=corrcoef(flyTruePrefPre,flyPredictedPrefPre);

figure;
plot(flyTruePref,flyPredictedPref,'o','LineWidth',3)
hold on
plot(flyTruePrefPre,flyPredictedPrefPre,'*','LineWidth',3)
legend('odor','pre-odor')
ylabel('Predicted Preference')
xlabel('Measured Preference')
legend boxoff
box off
%text(1.1*min(flyTruePref),0.9*max(flyPredictedPref),['r = ' num2str(rPredicted(1,2))],'FontSize',20)
%text(1.1*min(flyTruePref),0.8*max(flyPredictedPref),['r = ' num2str(rPredictedPre(1,2))],'FontSize',20)
set(gca,'FontSize',20)


% %validation points predicted
% figure
% plot(totaly,totalprediction,'o','LineWidth',3)
% hold on
% plot(totalyp,totalpredictionp,'*','LineWidth',3)
% legend('odor','pre-odor')
% legend boxoff
% box off
% ylabel('Predicted Preference')
% xlabel('Measured Preference')
% set(gca,'FontSize',15)


% correlation of held out observation predictions against true preference
nbins=30;
[yh xh]=hist(RON,nbins); yh=yh/sum(yh);
[yho xho]=hist(ROFF,nbins); yho=yho/sum(yho);

figure
plot(xh,yh,'LineWidth',2)
hold on
plot(xho,yho,'LineWidth',2)
legend('odor preference','pre-odor preference')
legend boxoff
box off
xlabel('correlation')
ylabel('frequency')
set(gca,'FontSize',20)

% show range of intercepts
[yb0 xb0]=hist(bont(1,:),nbins);yb0=yb0/sum(yb0);
[yb1 xb1]=hist(bont(2,:),nbins);yb1=yb1/sum(yb1);
[yb2 xb2]=hist(bont(3,:),nbins);yb2=yb2/sum(yb2);

[yob0 xob0]=hist(bofft(1,:),nbins);yob0=yob0/sum(yob0);
[yob1 xob1]=hist(bofft(2,:),nbins);yob1=yob1/sum(yob1);
[yob2 xob2]=hist(bofft(3,:),nbins);yob2=yob2/sum(yob2);

figure
subplot(1,2,1)
plot(xb0,yb0,'LineWidth',2)
hold on
plot(xb1,yb1,'LineWidth',2)
plot(xb2,yb2,'LineWidth',2)
legend('b0','b1','b2')
legend boxoff
box off
xlabel('weight')
ylabel('frequency')
title('odor preference')
set(gca,'FontSize',20)

subplot(1,2,2)
plot(xob0,yob0,'LineWidth',2)
hold on
plot(xob1,yob1,'LineWidth',2)
plot(xob2,yob2,'LineWidth',2)
legend('b0','b1','b2')
legend boxoff
box off
xlabel('weight')
ylabel('frequency')
title('pre-odor preference')
set(gca,'FontSize',20)

%% train prediction using differing numbers of pcs
warning('off')

iters=100;

trainsize=60;

RON=zeros(size(co,2),iters);
PON=zeros(size(co,2),iters);
ROFF=zeros(size(co,2),iters);
POFF=zeros(size(co,2),iters);

for pcstouse=1:size(co,2)
    predictor=(co(:,1:pcstouse)-mean(co(:,1:pcstouse),1))./std(co(:,1:pcstouse));
    predictor=predictor';
    
    bont=zeros(pcstouse+1,iters);
    bofft=zeros(pcstouse+1,iters);
    totaly=[];
    totalyp=[];
    totalprediction=[];
    totalpredictionp=[];
    
    for ii=1:iters
        
        temp=randperm(flyNum-1);
        trainflies=temp(1:floor(length(temp)*trainsize/100));
        validationflies=temp(ceil(length(temp)*trainsize/100):end);
        
        trainset=[];
        for j=1:length(trainflies)
            trainset=[trainset find(glombyodorflies==trainflies(j))];
            trainset=[trainset find(glombyodorflies==(trainflies(j)+leftRightSplitVal))];
        end
        
        validationset=[];
        for j=1:length(validationflies)
            validationset=[validationset find(glombyodorflies==validationflies(j))];
            validationset=[validationset find(glombyodorflies==(validationflies(j)+leftRightSplitVal))];
        end
        
        % Predict behavioral preference during odor
        % split data into train and test set
        bon=zeros(pcstouse+1,1);
        boff=zeros(pcstouse+1,1);
        
        % Predict behavioral preference during odor
        ytrain=behaviorOcc(trainset)';
        Xtr=predictor(:,trainset)';
        Xtrain=[ones(size(Xtr,1),1) Xtr];
        
        % run regression
        [bon,bint,r,rint,stats] = regress(ytrain,Xtrain);
        
        % Predict behavioral preference during preodor
        ytrainpre=behaviorpreOcc(trainset)';
        [boff,bintpre,rpre,rintpre,statspre] = regress(ytrainpre,Xtrain);
        
        % evaluate regression on held-out observations
        y=behaviorOcc(validationset)';
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
plot(mean(RON,2),'-','LineWidth',2)
ylabel('correlation')
hold on
yyaxis right
ylabel('p-value')
plot(mean(PON,2),'o-','LineWidth',2)
xlabel('# PCs')
set(gca,'FontSize',15)

fig=figure;
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(mean(ROFF,2),'-','LineWidth',2)
ylabel('correlation')
hold on
yyaxis right
ylabel('p-value')
plot(mean(POFF,2),'o-','LineWidth',2)
xlabel('# PCs')
set(gca,'FontSize',15)

%% train prediction using each pc
warning('off')

iters=100;

trainsize=30;

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
        
        temp=randperm(flyNum);
        trainflies=temp(1:floor(length(temp)*trainsize/100));
        validationflies=temp(ceil(length(temp)*trainsize/100):end);
        
        trainset=[];
        for j=1:length(trainflies)
            trainset=[trainset find(glombyodorflies==trainflies(j))];
            trainset=[trainset find(glombyodorflies==(trainflies(j)+leftRightSplitVal))];
        end
        
        validationset=[];
        for j=1:length(validationflies)
            validationset=[validationset find(glombyodorflies==validationflies(j))];
            validationset=[validationset find(glombyodorflies==(validationflies(j)+leftRightSplitVal))];
        end
        
        % Predict behavioral preference during odor
        % split data into train and test set
        bon=zeros(1,1);
        boff=zeros(1,1);
        
        % Predict behavioral preference during odor
        ytrain=behaviorOcc(trainset)';
        Xtr=predictor(:,trainset)';
        Xtrain=[ones(size(Xtr,1),1) Xtr];
        
        % run regression
        [bon,bint,r,rint,stats] = regress(ytrain,Xtrain);
        
        % Predict behavioral preference during preodor
        ytrainpre=behaviorpreOcc(trainset)';
        [boff,bintpre,rpre,rintpre,statspre] = regress(ytrainpre,Xtrain);
        
        % evaluate regression on held-out observations
        y=behaviorOcc(validationset)';
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
plot(mean(RON,2),'-','LineWidth',2)
ylabel('correlation')
hold on
yyaxis right
ylabel('p-value')
plot(mean(PON,2),'o-','LineWidth',2)
xlabel('# PCs')
set(gca,'FontSize',15)

fig=figure;
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(mean(ROFF,2),'-','LineWidth',2)
ylabel('correlation')
hold on
yyaxis right
ylabel('p-value')
plot(mean(POFF,2),'o-','LineWidth',2)
xlabel('# PCs')
set(gca,'FontSize',15)

%% run tensor analysis on tensor of responses
% Fit CP Tensor Decomposition

% these commands require that you download Sandia Labs' tensor toolbox:
% http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.6.html
responsesTensorOrig=responsesTensorFullyOrthogonal;
gNames=publishedOR.gh146glomerulusNames;
fullyorthogonal=1; % if using a fully unwrapped tensor (all variables along independent axes)

% convert data to a tensor object
responsesNoResponseRemoved=responsesTensorOrig;
fracIn=0.025;
numFinite=zeros(1,size(responsesTensorFullyOrthogonal,3));
for j=1:size(responsesTensorFullyOrthogonal,3)
    
numFinite(j)=sum(sum(isfinite(responsesTensorOrig(:,:,j,1,1))));
end
toRemove=find(numFinite/size(responsesTensorOrig,2)<fracIn);

if fullyorthogonal
    responsesNoResponseRemoved=responsesTensorFullyOrthogonal;
end

responsesNoResponseRemoved(:,:,toRemove,:,:)=[];

finalGlomNames=gNames;
finalGlomNames(toRemove)=[];

% set NaNs equal to mean of non NaN values for the same glomerulus/time
% point
if ~fullyorthogonal
    for i=1:size(responsesNoResponseRemoved,1)
        for j=1:size(responsesNoResponseRemoved,2)
            for k=1:size(responsesNoResponseRemoved,3)
                if isnan(responsesNoResponseRemoved(i,j,k))
                    responsesNoResponseRemoved(i,j,k)=nanmean(responsesNoResponseRemoved(i,:,k));
                end
            end
        end
    end
else
    for i=1:size(responsesNoResponseRemoved,1)
        for j=1:size(responsesNoResponseRemoved,2)
            for k=1:size(responsesNoResponseRemoved,3)
                for l=1:size(responsesNoResponseRemoved,4)
                    for m=1:size(responsesNoResponseRemoved,5)
                        if isnan(responsesNoResponseRemoved(i,j,k,l,m))
                            temp=responsesNoResponseRemoved(:,:,k,l,m);
                            responsesNoResponseRemoved(i,j,k,l,m)=nanmean(temp(:));
                        end
                    end
                end
            end
        end
    end
end

data = tensor(responsesNoResponseRemoved);

% plot the ground truth
% true_factors = ktensor(lam, A, B, C);
% true_err = norm(full(true_factors) - data)/norm(true_factors);
% viz_ktensor(true_factors, ...
%             'Plottype', {'bar', 'line', 'scatter'}, ...
%             'Modetitles', {'neurons', 'time', 'trials'})
% set(gcf, 'Name', 'true factors')

% fit the cp decomposition from random initial guesses
n_fits = 10;
R=5;
tcafits=cell(n_fits,1);
err = zeros(n_fits,length(R));
for n = 1:n_fits
    for i=1:length(R)
        disp(['fitting R = ' num2str(R(i)) ', nfit = ' num2str(n) ' of ' num2str(n_fits)])
        % fit model
        est_factors = cp_als(data,R(i),'printitn',0);
        
        tcafits{n}{i}=est_factors;
        % store error
        err(n,i) = norm(full(est_factors) - data)/norm(data);
        
        % visualize fit for first several fits
        if n == 1
            % score aligns the cp decompositions
            %[sc, est_factors] = score(est_factors, true_factors);
            
            % plot the estimated factors
            if ~fullyorthogonal
                viz_ktensor(est_factors, ...
                    'Plottype', {'bar', 'line', 'scatter'}, ...
                    'Modetitles', {'glomeruli', 'flyxodor', 'time'})
                set(gcf, 'Name', ['estimated factors - fit #' num2str(n)])
            else
                viz_ktensor(est_factors, ...
                    'Plottype', {'bar','bar','bar', 'line', 'scatter'}, ...
                    'Modetitles', {'fly #', 'L/R lobe', 'glomerulus','odor','time'})
                set(gcf, 'Name', ['estimated factors - fit #' num2str(n)])
            end
        end
    end
end

%est_factors.lambda
%est_factors.U