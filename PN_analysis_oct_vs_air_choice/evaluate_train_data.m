
clear all
%close all
rng('default')

load ORN_PN_colors

load analysis_dir_path

manualLabelHome=fullfile(analysis_dir_path, 'PN_analysis_oct_vs_air_choice/train');

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


%% perform pca on responses

clear responsesNoResponseRemoved

fracIn=0.1; % best results when fracIn is high, ~0.5, only using high confidence glomeruli


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
finiteFrac = numFinite/size(responsesNoResponseRemoved,2);
glomFiniteFrac = finiteFrac(1:13:end);
responsesNoResponseRemoved(toRemove,:)=[];

% remove air and ethanol
responsesNoResponseRemoved(1:13:end)=0;
responsesNoResponseRemoved(8:13:end)=0;

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

% REMOVE D
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

yesZscore=0;
if yesZscore
    responsesNoResponseRemoved=responsesNoResponseRemoved';
    responsesNoResponseRemoved=(responsesNoResponseRemoved-mean(responsesNoResponseRemoved))./std(responsesNoResponseRemoved);
    responsesNoResponseRemoved=responsesNoResponseRemoved';
end

opt = statset('pca');
opt.Display='iter';
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(responsesNoResponseRemoved','Options',opt);


figure; %1
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
figure; %2
imagesc(glomcontributionMean)
set(gca,'ytick',1:length(gNames),'yticklabel',string(gNames),'FontSize',10)
ytickangle(30)
xlabel('PC #')
set(gca,'FontSize',15)

%%
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
figure %3
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

figure %4
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

figure %5
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

%% use fitlm

pcstouse=[1];

behaviorprediction=(SCORE(:,pcstouse));
flyTruePref=zeros(1,flyNum);
flyPredictedPref=zeros(1,flyNum);
ally=behaviorOcc';
linmodel=fitlm(behaviorprediction,ally);
myprediction=predict(linmodel,behaviorprediction);
figure %6
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
figure %7
plot(myprediction,flyTruePref,'.','Color',pcolor, 'LineWidth',3)
for i=1:flyNum
   hold on
   %text(myprediction(i)+0.01,flyTruePref(i),num2str(i),'FontSize',15)
end
[r p]=corrcoef(myprediction,flyTruePref);
text(-.25,-.25,['r = ' num2str(r(1,2),'%2.2f')],'FontSize',15)
text(-.3,-.25,['p = ' num2str(p(1,2),'%2.2f')],'FontSize',15)
xlabel('Predicted Preference')
ylabel('Measured Preference')
set(gca,'FontSize',15)
box on
axis([-.6 .3 -.6 .3])
axis square
linmodel

% SUP FIG 9a
figure %8
hold on;
xVals = (myprediction-mean(myprediction))/(std(myprediction));
yVals = (flyTruePref-mean(flyTruePref))/(std(flyTruePref));
linreg = linearRegressionCI2(xVals, yVals.', 1, 0, -2, 3);

areaBar(linreg.xVals,polyval(linreg.pOverall,linreg.xVals),2*std(linreg.fits),[0 0 0],[0.9 0.9 0.9])
plot(xVals,yVals,'.','Color',pcolor, 'LineWidth',3)
for i=1:flyNum
   hold on
   %text(myprediction(i)+0.01,flyTruePref(i),num2str(i),'FontSize',15)
end
[r p]=corrcoef(xVals,yVals);
text(-.25,-.25,['r = ' num2str(r(1,2),'%2.2f')],'FontSize',15)
text(-.3,-.25,['p = ' num2str(p(1,2),'%2.2f')],'FontSize',15)
set(gca,'FontSize',15)
box on
axis([-2 3 -2 3])
xlabel('predicted preference (z-scored)')
ylabel('measured preference (z-scored)')
set(gca,'xtick','')
set(gca,'ytick','')
axis square







beta=linmodel.Coefficients.Estimate;

PCContribution=COEFF(:,pcstouse);
figure; %9
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
ylabel('PC loadings')
box off
set(gca,'FontSize',15)

currpc=PCContribution;

%save trainpc currpc gNames
%save trainDataModel linmodel currpc gNames

%% use average across all dimensions
behaviorprediction=transpose(mean(responsesNoResponseRemoved));
flyTruePref=zeros(1,flyNum);
flyPredictedPref=zeros(1,flyNum);
ally=behaviorOcc';
linmodel=fitlm(behaviorprediction,ally);
myprediction=predict(linmodel,behaviorprediction);
figure %10
plot(myprediction,ally,'o','LineWidth',3)
xlabel('Predicted Preference')
ylabel('Measured Preference')
set(gca,'FontSize',15)
box off
linmodel

nactivity=zeros(flyNum,1);
for i=1:flyNum
   flyTruePref(i)=mean(ally(flyindices{i}));
   flyPredictedPref(i)=mean(mean(behaviorprediction(flyindices{i},:)));
   nactivity(i,:)=mean(behaviorprediction(flyindices{i},:));
end
linmodel=fitlm(nactivity,flyTruePref);
myprediction=predict(linmodel,nactivity);
figure %11
plot(myprediction,flyTruePref,'o','Color',pcolor, 'LineWidth',3)
for i=1:flyNum
   hold on
   %text(myprediction(i)+0.01,flyTruePref(i),num2str(i),'FontSize',15)
end
[r p]=corrcoef(myprediction,flyTruePref);
text(-.25,-.25,['r = ' num2str(r(1,2),'%2.2f')],'FontSize',15)
text(-.3,-.25,['p = ' num2str(p(1,2),'%2.2f')],'FontSize',15)
xlabel('Predicted Preference')
ylabel('Measured Preference')
set(gca,'FontSize',15)
box off
linmodel

beta=linmodel.Coefficients.Estimate;


%save trainpc currpc gNames
%save trainDataModel_averageAllDimensions linmodel gNames

%% bootstrap
warning('off')

iters=250;
bootstrapsamples=flyNum;
highestPCtouse=5;
yesZscore=0; % zscore data matrix before performing PCA (0=no, 1=yes)

fracIn=0.1; % best results when fracIn is high, ~0.5, only using high confidence glomeruli
medianResponseOrTimeCourse=1; % 1 for median response only, 0 for full time course

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

% remove D

gNames(1)=[];
responsesNoResponseRemoved(1:13,:)=[];
% remove air and ethanol
%responsesNoResponseRemoved(1:nodors:end,:)=0;
%responsesNoResponseRemoved(8:nodors:end,:)=0;

% % fill nans with mean

for i=1:size(responsesNoResponseRemoved,1)
    for j=1:size(responsesNoResponseRemoved,2)
        if isnan(responsesNoResponseRemoved(i,j))
            responsesNoResponseRemoved(i,j)=nanmean(responsesNoResponseRemoved(i,:));
        end
    end
end

testR=zeros(highestPCtouse,iters);
testRshuffled=zeros(highestPCtouse,iters);
testR2=zeros(highestPCtouse,iters);
testR2shuffled=zeros(highestPCtouse,iters);
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
    
    % select bootstrap sample flies
    trainflies=zeros(1,bootstrapsamples);
    for i=1:bootstrapsamples
        trainflies(i)=1+round((flyNum-1)*rand(1));
    end
 
    flyindicestrain=cell(1,length(trainflies));
    datapts=cell(1,length(trainflies));
    lastdatapt=0;
    for i=1:length(trainflies)
        [temp temp2]=find(glombyodorflies==trainflies(i));
        [tem3 temp4]=find(glombyodorflies==(trainflies(i)+leftRightSplitVal));
        flyindicestrain{i}=[temp2 temp4];
        datapts{i}=[(lastdatapt+1):(lastdatapt+length(flyindicestrain{i}))];
        lastdatapt=(lastdatapt+length(flyindicestrain{i}));
    end
    
    traindata=[];
    for i=1:length(trainflies)
        for j=1:length(flyindicestrain{i})
            traindata=[traindata responsesNoResponseRemoved(:,flyindicestrain{i}(j))];
        end
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
        ally=behaviorOcc;
        
        nactivity=zeros(length(trainflies),length(pcstouse));
        shuffledindtrain=randperm(length(trainflies)); % for comparing to shuffled null model
        for i=1:length(trainflies)
            flyTruePref(i)=mean(ally(flyindicestrain{i}));
            flyTruePrefShuffled(i)=mean(ally(flyindicestrain{shuffledindtrain(i)}));
            nactivity(i,:)=mean(behaviorprediction(datapts{i},:));
        end
        linmodel=fitlm(nactivity,flyTruePref);
        linmodelShuffled=fitlm(nactivity,flyTruePrefShuffled);
        
        beta=linmodel.Coefficients.Estimate;
        
        PCContribution=currpc*beta(2:end);
        
        mycorr=corrcoef(nactivity,flyTruePref);
        mycorrShuffled=corrcoef(nactivity,flyTruePrefShuffled);
        
        myr2=linmodel.Rsquared.Ordinary;
        myr2shuffled=linmodelShuffled.Rsquared.Ordinary;
        
        if jjj==1
            %averagePredictor{pcstouse}=PCContribution';
            averagePredictor{pcstouse}=currpc';
        else
            %averagePredictor{pcstouse}=averagePredictor{pcstouse}+PCContribution';
            averagePredictor{pcstouse}=averagePredictor{pcstouse}+currpc';
        end
        
        testR(pcstouse,jjj)=mycorr(1,2);
        testRshuffled(pcstouse,jjj)=mycorrShuffled(1,2);
        testR2(pcstouse,jjj)=myr2;
        testR2shuffled(pcstouse,jjj)=myr2shuffled;
        
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
    
    if mod(jjj,10)==0
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
% 
% figure
% subplot(2,1,1)
% distributionPlot(testR','histOpt',1,'colormap',1-gray(64),'showMM',0)
% xlabel('PC used for linear model')
% ylabel('Unshuffled R')
% box off
% set(gca,'FontSize',15)
% subplot(2,1,2)
% distributionPlot(testRshuffled','histOpt',1,'colormap',1-gray(64),'showMM',0)
% xlabel('PC used for linear model')
% ylabel('Shuffled R')
% box off
% set(gca,'FontSize',15)
% 
% figure
% subplot(1,2,1)
% distributionPlot(testR2','histOpt',1,'colormap',1-gray(64),'showMM',0)
% xlabel('PC used for linear model')
% ylabel('Unshuffled R^2')
% axis([0 highestPCtouse+1 0 0.8])
% box off
% set(gca,'FontSize',15)
% subplot(1,2,2)
% distributionPlot(testR2shuffled','histOpt',1,'colormap',1-gray(64),'showMM',0)
% xlabel('PC used for linear model')
% ylabel('Shuffled R^2')
% axis([0 highestPCtouse+1 0 0.8])
% box off
% set(gca,'FontSize',15)

% figure 12
violinPlot(testR2',pcolor)
xlabel('PC used for linear model')
ylabel('Unshuffled R^2')
axis([0 highestPCtouse+1 0 1])
box off
set(gca,'FontSize',15)
% figure 13
violinPlot(testR2shuffled',pcolor)
xlabel('PC used for linear model')
ylabel('Shuffled R^2')
axis([0 highestPCtouse+1 0 1])
box off
set(gca,'FontSize',15)


testR2t=testR2';
testR2shuffledt=testR2shuffled';
labs=[ones(1,iters) 2*ones(1,iters) 3*ones(1,iters) 4*ones(1,iters) 5*ones(1,iters)];
% FIG 1k PN OCT-AIR preference prediction
figure %14
boxplot(testR2t(:),labs,'plotstyle','compact','BoxStyle','filled','Colors',pcolor,'medianstyle','target','symbol','','outliersize',1)
xlabel('PC used for linear model')
ylabel('Unshuffled R^2')
axis([0 highestPCtouse+1 0 0.75])
set(gca,'FontSize',15)
set(gca,'xtick','')
set(gca,'ytick','')
figure %15
boxplot(testR2shuffledt(:),labs,'plotstyle','compact','BoxStyle','filled','Colors',pcolor,'medianstyle','target','symbol','','outliersize',1)
xlabel('PC used for linear model')
ylabel('Shuffled R^2')
axis([0 highestPCtouse+1 0 0.75])
set(gca,'FontSize',15)
set(gca,'xtick','')
set(gca,'ytick','')

%% run train/test many times and save replicate data
warning('off')

iters=50;
highestPCtouse=5;
yesZscore=0; % zscore data matrix before performing PCA (0=no, 1=yes)
normalizeResponse=0;

fracIn=0.4; % best results when fracIn is high, ~0.5, only using high confidence glomeruli
medianResponseOrTimeCourse=1; % 1 for median response only, 0 for full time course

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

% remove air and ethanol
%responsesNoResponseRemoved(1:nodors:end,:)=0;
%responsesNoResponseRemoved(8:nodors:end,:)=0;


% normalize responsesNoResponseRemoved
if normalizeResponse
    for i=1:size(responsesNoResponseRemoved,2)
        responsesNoResponseRemoved(:,i)=responsesNoResponseRemoved(:,i)/nanmean(abs(responsesNoResponseRemoved(:,i)));
    end
end

testR=zeros(highestPCtouse,iters);
testRshuffled=zeros(highestPCtouse,iters);
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
    testsize=50;
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
        myrshuffled=corrcoef(mytestprediction,testflyTruePref(randperm(length(testflyTruePref))));
        
        if jjj==1
            %averagePredictor{pcstouse}=PCContribution';
            averagePredictor{pcstouse}=currpc';
        else
            %averagePredictor{pcstouse}=averagePredictor{pcstouse}+PCContribution';
            averagePredictor{pcstouse}=averagePredictor{pcstouse}+currpc';
        end
        
        testR(pcstouse,jjj)=myr(1,2);
        testRshuffled(pcstouse,jjj)=myrshuffled(1,2);
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
    
    if mod(jjj,10)==0
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


myr2=(testR.^2').*sign(testR');
myr2shuffled=(testRshuffled.^2').*sign(testRshuffled');

% figure
% subplot(1,2,1)
% distributionPlot(myr2,'histOpt',1,'colormap',1-gray(64),'showMM',0)
% xlabel('PC used for linear model')
% ylabel('Unshuffled R^2')
% axis([0 highestPCtouse+1 -0.95 0.95])
% box off
% set(gca,'FontSize',15)
% subplot(1,2,2)
% distributionPlot(myr2shuffled,'histOpt',1,'colormap',1-gray(64),'showMM',0)
% xlabel('PC used for linear model')
% ylabel('Shuffled R^2')
% axis([0 highestPCtouse+1 -0.95 0.95])
% box off
% set(gca,'FontSize',15)
% 
% 
figure %16
distributionPlot(myr2,'histOpt',1,'colormap',pcolor,'showMM',0)
xlabel('PC used for linear model')
ylabel('Unshuffled R^2')
axis([0 highestPCtouse+1 -1 1])
box off
set(gca,'FontSize',15)
figure %17
distributionPlot(myr2shuffled,'histOpt',1,'colormap',pcolor,'showMM',0)
xlabel('PC used for linear model')
ylabel('Shuffled R^2')
axis([0 highestPCtouse+1 -1 1])
box off
set(gca,'FontSize',15)

% figure 18
violinPlot(myr2,pcolor)
xlabel('PC used for linear model')
ylabel('Unshuffled R^2')
axis([0 highestPCtouse+1 -1 1])
box off
set(gca,'FontSize',15)
% figure 19
violinPlot(myr2shuffled,pcolor)
xlabel('PC used for linear model')
ylabel('Shuffled R^2')
axis([0 highestPCtouse+1 -1 1])
box off
set(gca,'FontSize',15)
% 
% figure
% subplot(2,2,1)
% boxplot(testR')
% xlabel('PC used for linear model')
% ylabel('Correlation between predicted vs. true preference (test data)')
% box off
% set(gca,'FontSize',15)
% 
% subplot(2,2,2)
% hist(bestPredictor)
% xlabel('Best predictor PC #')

% boxplot(pcPredictorRank)
% xlabel('PC #')
% ylabel('PC Predictor Rank')
% 
% subplot(2,2,3)
% hist(bestPredictorR)
% xlabel('Best predictor''s r-value')
% 
% bploadings=zeros(1,size(COEFF,1));
% allbpl=zeros(iters,size(COEFF,1));
% for j=1:iters
%     bploadings=bploadings+bestPredictorLoading{j}';
%     allbpl(j,:)=bestPredictorLoading{j}';
% end
% bploadings=bploadings/iters;
% 
% 
% subplot(2,2,4)
% plot(bploadings,'*','LineWidth',3,'MarkerSize',10)
% hold on
% plot(zeros(1,size(COEFF,1)),'k--','LineWidth',3)
% j=1;
% for i=1:nodors:length(averagePredictor{4})
%     plot((i-0.5)*ones(1,100), linspace(min(bploadings),max(bploadings)),'k--','LineWidth',2)
%     text(i+floor(nodors/3),min(bploadings),num2str(glomsFound(j)),'FontSize',15)
%     j=j+1;
% end
% set(gca,'xtick',(1:nodors:length(bploadings))+floor(nodors/2),'xticklabel',string(gNames),'FontSize',10)
% xtickangle(30)
% ylabel('Best predictor''s loadings')
% box off
% set(gca,'FontSize',15)
% % 
% % figure
% % boxplot(relativeOctMchactivation)
% % xlabel('PC #')
% % ylabel('Abs(Oct-MCH) Loading')
% % set(gca,'FontSize',15)
% %
% % figure
% % hist(relativeOctMchactivation(:,4),50)
% % xlabel('PC #')
% % ylabel('Abs(Oct-MCH) Loading')
% % set(gca,'FontSize',15)
% %
% % figure
% % distributionPlot(relativeOctMchactivation,'histOpt',1,'colormap',1-gray(64));
% % xlabel('PC #')
% % ylabel('Average OCT-MCH Loading')
% % set(gca,'FontSize',15)
% 
% % show best predictor loadings for each iteration
% bpls=zeros(length(bestPredictorLoading{1}),iters);
% oa=zeros(1,iters);
% for j=1:iters
%     bpls(:,j)=bestPredictorLoading{j}';
%     oa(j)=corr(bestPredictorLoading{j}(2:13:end),bestPredictorLoading{j}(11:13:end));
% end
% figure;
% imagesc(bpls)
% xlabel('iteration')
% ylabel('dimension')
% 
% figure; hist(oa)
% 
% ex=zeros(1,30);
% for j=1:iters
%     ex=ex+vExplained{j}(1:30)';
% end
% ex=ex/iters;
% figure
% plot(log10(ex),'LineWidth',3)
% xlabel('PC #')
% ylabel('log10(Variance Explained)')
% set(gca,'FontSize',20)
% 
% % run kmeans on allbpl (how many "best pcs" are there)
% [IDX, bestpredictorC, SUMD]= kmeans(allbpl, 20);
% figure;
% imagesc(bestpredictorC)
% xlabel('Dimensions')
% ylabel('Predictor Class (from kmeans)')
% title('Best PC predictor loadings (all)')
% set(gca,'FontSize',15)
% 
% % make a "best predictor matrix"
% bestpredictors=bestpredictorC(1,:);
% for i=2:size(bestpredictorC,1)
%     for j=1:size(bestpredictors,1)
%         temp=corrcoef(bestpredictorC(i,:),bestpredictors(j,:));
%         temp2=temp(1,2);
%         if abs(temp2)>0.4
%             if temp2>0
%                 bestpredictors(j,:)=(bestpredictors(j,:)+bestpredictorC(i,:))/2;
%             else
%                 bestpredictors(j,:)=(bestpredictors(j,:)-bestpredictorC(i,:))/2;
%             end
%             break
%         elseif j==size(bestpredictors,1)
%             bestpredictors=[bestpredictors; bestpredictorC(i,:)];
%         end
%         
%     end
% end
% 
% figure;
% imagesc(bestpredictors)
% xlabel('Dimensions')
% ylabel('Predictor Class (from kmeans)')
% title('Best PC predictor loadings (clustered)')
% set(gca,'FontSize',15)
% 
% 
% % run kmeans on each rawPC sweep to find possibilities
% maxK=20;
% for i=1:highestPCtouse
%     [IDX, C, SUMD]= kmeans(rawPC{i}', maxK);
%     kidx{i}=IDX;
%     kc{i}=C;
% end
% figure
% for i=1:highestPCtouse
%     subplot(2,5,i)
%     boxplot(testR(i,:),kidx{i})
%     ylabel('Measured-Predicted Fit Correlation')
%     title(['PC ' num2str(i)])
% end
% 
% figure
% for i=1:highestPCtouse
%     subplot(2,5,i)
%     
%     boxplot(octmchloadingcorr(:,i),kidx{i})
%     ylabel('OCT-MCH correlation')
%     title(['PC ' num2str(i)])
% end
% 
% % go through all PCs and select each one which led to a behavior prediction
% predcorr=0.2; % select PCs with behavior prediction correlation higher than this number
% bestOfPc=[];
% for i=1:highestPCtouse
%     for j=1:maxK
%         if median(testR(i,kidx{i}==j))>predcorr
%             if bestOfPc
%                 bestOfPc=kc{i}(j,:);
%             else
%                 bestOfPc=[bestOfPc; kc{i}(j,:)];
%             end
%             disp(['saving PC# ' num2str(i) ' cluster ' num2str(j)])
%         end
%     end
% end
% 
% % make a "best of predictor matrix"
% bestpredictorsAll=bestOfPc(1,:);
% for i=2:size(bestOfPc,1)
%     for j=1:size(bestpredictorsAll,1)
%         temp=corrcoef(bestOfPc(i,:),bestpredictorsAll(j,:));
%         temp2=temp(1,2);
%         if abs(temp2)>0.6
%             if temp2>0
%                 bestpredictorsAll(j,:)=(bestpredictorsAll(j,:)+bestOfPc(i,:))/2;
%             else
%                 bestpredictorsAll(j,:)=(bestpredictorsAll(j,:)-bestOfPc(i,:))/2;
%             end
%             break
%         elseif j==size(bestpredictorsAll,1)
%             bestpredictorsAll=[bestpredictorsAll; bestOfPc(i,:)];
%         end
%         
%     end
% end
% 
% figure;
% imagesc(bestpredictorsAll)
% xlabel('Dimensions')
% ylabel('Predictor Class (from kmeans)')
% title('Best PC predictor loadings (clustered)')
% set(gca,'FontSize',15)


