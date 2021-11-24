
clear all
%close all
rng('default')

load ORN_PN_colors
manualLabelHome='/Users/mattchurgin/Dropbox (Harvard University)/flyimaging/analysis/PN_analysis_oct_vs_air_choice/alldata';

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

%% use fitlm

pcstouse=[1];

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
   
   flypc(i,:)=mean(SCORE(flyindices{i},:),1);
end
linmodel=fitlm(nactivity,flyTruePref);
myprediction=predict(linmodel,nactivity);
figure
plot((myprediction-mean(myprediction))/std(myprediction),(flyTruePref-mean(flyTruePref))/std(flyTruePref),'.','Color',pcolor,'LineWidth',2)
xlabel('Predicted Preference')
ylabel('Measured Preference')
set(gca,'FontSize',15)
box on
axis square
axis([-5 4 -5 4])
set(gca,'xtick','')
set(gca,'ytick','')
linmodel

beta=linmodel.Coefficients.Estimate;

PCContribution=COEFF(:,pcstouse);
figure;
plot(PCContribution,'.','Color',pcolor,'LineWidth',2,'MarkerSize',5)
hold on
plot(zeros(1,length(PCContribution(:,1))),'k--','LineWidth',3)
j=1;
for i=nodors:nodors:(length(PCContribution)-1)
    plot((i+0.5)*ones(1,5), linspace(2*min(PCContribution),2*max(PCContribution),5),'--','Color',[0.7 0.7 0.7],'LineWidth',2)
    %text(i+floor(nodors/3),min(PCContribution),num2str(glomsFound(j)),'FontSize',15)
    j=j+1;
end

set(gca,'xtick',(1:nodors:length(PCContribution))+floor(nodors/2),'xticklabel',string(gNames),'FontSize',10)
xtickangle(30)
ylabel('PC 1 loadings')
axis([0 66 -.01 .32])
set(gca,'FontSize',15)

figure;
bar(PCContribution,'FaceColor',pcolor)
hold on
plot(zeros(1,length(PCContribution(:,1))),'k--','LineWidth',3)
j=1;
for i=nodors:nodors:(length(PCContribution)-1)
    plot((i+0.5)*ones(1,5), linspace(2*min(PCContribution),2*max(PCContribution),5),'--','Color',[0.7 0.7 0.7],'LineWidth',2)
    %text(i+floor(nodors/3),min(PCContribution),num2str(glomsFound(j)),'FontSize',15)
    j=j+1;
end

set(gca,'xtick',(1:nodors:length(PCContribution))+floor(nodors/2),'xticklabel',string(gNames),'FontSize',10)
set(gca,'xtick','')
set(gca,'ytick','')
xtickangle(30)
ylabel('PC 1 loadings')
axis([0 66 -.05 .32])
set(gca,'FontSize',15)

figure;
bar(ones(1,65),'FaceColor',pcolor)
hold on
plot(zeros(1,length(PCContribution(:,1))),'k--','LineWidth',3)
j=1;
for i=nodors:nodors:(length(PCContribution)-1)
    plot((i+0.5)*ones(1,5), linspace(-.05,1.2,5),'--','Color',[0.7 0.7 0.7],'LineWidth',2)
    %text(i+floor(nodors/3),min(PCContribution),num2str(glomsFound(j)),'FontSize',15)
    j=j+1;
end

set(gca,'xtick',(1:nodors:length(PCContribution))+floor(nodors/2),'xticklabel',string(gNames),'FontSize',10)
set(gca,'xtick','')
set(gca,'ytick','')
xtickangle(30)
ylabel('PC 1 loadings')
axis([0 66 -.05 1.2])
set(gca,'FontSize',15)

currpc=PCContribution;

gbcm=interp1([1 256],[0 0 0; 0 1 0],1:256);
figure
scatter(flypc(:,1),flypc(:,2),50,flyTruePref,'filled')
colormap(gbcm)
box on
xlabel(['PC 1 (' num2str(EXPLAINED(1),'%2.1f') '%)'])
ylabel(['PC 2 (' num2str(EXPLAINED(2),'%2.1f') '%)'])
colorbar
axis square
set(gca,'FontSize',15)


%save trainpc currpc gNames
%save trainDataModel linmodel
%% use average response across all dimensions


behaviorprediction=transpose(mean(responsesNoResponseRemoved));
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

nactivity=zeros(flyNum,1);
for i=1:flyNum
   flyTruePref(i)=mean(ally(flyindices{i}));
   flyPredictedPref(i)=mean(mean(behaviorprediction(flyindices{i},:)));
   nactivity(i,:)=mean(behaviorprediction(flyindices{i},:));
end




% plot histogram
figure
histogram(nactivity,10)
ylabel('# flies')
xlabel('average df/f')
axis square

% plot raw values
figure
plot(nactivity,flyTruePref,'.','Color',pcolor, 'LineWidth',3,'MarkerSize',20)
xlabel('average df/f')
ylabel('measured preference')
axis square

linmodel=fitlm(nactivity,flyTruePref);
myprediction=predict(linmodel,nactivity);
figure
plot(myprediction,flyTruePref,'o','Color',[0 0.7 0],'LineWidth',3)
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



figure
plot((myprediction-mean(myprediction))/(std(myprediction)),(flyTruePref-mean(flyTruePref))/(std(flyTruePref)),'.','Color',pcolor, 'LineWidth',3,'MarkerSize',20)
for i=1:flyNum
   hold on
   %text(myprediction(i)+0.01,flyTruePref(i),num2str(i),'FontSize',15)
end
[r p]=corrcoef((myprediction-mean(myprediction))/(std(myprediction)),(flyTruePref-mean(flyTruePref))/(std(flyTruePref)));
text(-.25,-.25,['r = ' num2str(r(1,2),'%2.2f')],'FontSize',15)
text(-.3,-.25,['p = ' num2str(p(1,2),'%2.2f')],'FontSize',15)
set(gca,'FontSize',15)
box on
axis([-4.5 3.5 -4.5 3.5])
xlabel('predicted preference (z-scored)')
ylabel('measured preference (z-scored)')
set(gca,'xtick','')
set(gca,'ytick','')
axis square

%save trainpc currpc gNames
%save trainDataModel linmodel

