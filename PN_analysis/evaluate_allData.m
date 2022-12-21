% look at manually labelled gh146 flies in odor space
clear all
%close all

load ORN_PN_colors

load analysis_dir_path

manualLabelHome=fullfile(analysis_dir_path, 'PN_analysis/alldata');

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

datNames=[];
for i=1:length(manualLabelledFolders)
    currname=manualLabelledFolders(i).name;
    datNames = [datNames, '\n', currname];
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
finiteFrac = numFinite/size(responsesNoResponseRemoved,2);
glomFiniteFrac = finiteFrac(1:13:end);
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
% % fill nans with mean
for i=1:size(responsesNoResponseRemoved,1)
    for j=1:size(responsesNoResponseRemoved,2)
        if isnan(responsesNoResponseRemoved(i,j))
            responsesNoResponseRemoved(i,j)=nanmean(responsesNoResponseRemoved(i,:));
        end
    end
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
%% plot individual flies (left v right lobe) in PC space
varianceToKeep=80; % percent of variance to keep


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


%% build linear model from PC 2 to predict behavior

pcstouse=[2];

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
    
    flypc(i,:)=mean(SCORE(flyindices{i},:),1);
end
linmodel=fitlm(nactivity,flyTruePref);
myprediction=predict(linmodel,nactivity);

% FIG 2g inset
figure %7
hold on;
xVals = (myprediction-mean(myprediction))/std(myprediction);
yVals = (flyTruePref-mean(flyTruePref))/std(flyTruePref);
linreg = linearRegressionCI2(xVals, yVals.', 1, 0, -3, 2.3);
areaBar(linreg.xVals,polyval(linreg.pOverall,linreg.xVals),2*std(linreg.fits),[0 0 0],[0.9 0.9 0.9])
plot(xVals,yVals,'.','Color',pcolor,'LineWidth',2)
xlabel('Predicted Preference')
ylabel('Measured Preference')
set(gca,'FontSize',15)
box on
axis square
axis([-3 2.3 -3 2.3])
set(gca,'xtick','')
set(gca,'ytick','')
linmodel

beta=linmodel.Coefficients.Estimate;

PCContribution=COEFF(:,pcstouse);
figure; %8
plot(PCContribution,'.','Color',pcolor,'LineWidth',2,'MarkerSize',5)
hold on
plot(zeros(1,length(PCContribution(:,1))),'k--','LineWidth',3)
j=1;
for i=nodors:nodors:(length(PCContribution)-1)
    plot((i+0.5)*ones(1,5), linspace(2*min(PCContribution),2*max(PCContribution),5),'--','Color',[0.7 0.7 0.7],'LineWidth',2)
    %text(i+floor(nodors/3),min(PCContribution),num2str(glomsFound(j)),'FontSize',15)
    j=j+1;
end

%set(gca,'xtick',(1:nodors:length(PCContribution))+floor(nodors/2),'xticklabel',string(gNames),'FontSize',10)
xtickangle(30)
ylabel('PC 2 loadings')
axis([0 66 -.34 .34])
set(gca,'FontSize',15)
set(gca,'xtick','')
set(gca,'ytick','')
currpc=COEFF(:,pcstouse);

% FIG 2g loadings
figure; %9
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
ylabel('PC 2 loadings')
axis([0 66 -.34 .34])
set(gca,'FontSize',15)

octmchselect=ones(1,65);
octmchselect([14:39 53:65])=0;
octmchselect(40:52)=-1;
% FIG 2h
figure; %10
bar(octmchselect,'FaceColor',pcolor)
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
ylabel('PC 2 loadings')
axis([0 66 -1.1 1.1])
set(gca,'FontSize',15)


%nyidalur=interp1([1 128 129 256],[1 0 1; .3091 .2303 .3091; .2306 .3082 .2306; 0 1 0],1:256);
skalafell=interp1([1 128 129 256],[1 0 1; 1 .9921 1; .9921 1 .9921; 0 1 0],1:256);
figure %11
scatter(flypc(:,1),flypc(:,2),50,flyTruePref,'filled')
colormap(skalafell)
box on
xlabel(['PC 1 (' num2str(EXPLAINED(1),'%2.1f') '%)'])
ylabel(['PC 2 (' num2str(EXPLAINED(2),'%2.1f') '%)'])
colorbar
axis square
set(gca,'FontSize',15)

%% use fitlm to fit left and right lobe separately

pcstouse=[2];

behaviorprediction=(SCORE(:,pcstouse));
flyTruePrefL=zeros(1,flyNum);
flyTruePrefR=zeros(1,flyNum);

ally=behaviorOcc';

nactivityL=zeros(flyNum,length(pcstouse));
nactivityR=zeros(flyNum,length(pcstouse));
for i=1:flyNum
    flyTruePrefL(i)=mean(ally(flyindicesL{i}));
    nactivityL(i,:)=mean(behaviorprediction(flyindicesL{i},:));
    flyTruePrefR(i)=mean(ally(flyindicesR{i}));
    nactivityR(i,:)=mean(behaviorprediction(flyindicesR{i},:));
end
linmodelL=fitlm(nactivityL,flyTruePrefL);
linmodelR=fitlm(nactivityR,flyTruePrefR);
mypredictionL=predict(linmodelL,nactivityL);
mypredictionR=predict(linmodelR,nactivityR);
figure %12
plot(mypredictionL,flyTruePrefL,'*','LineWidth',3,'Markersize',10)
hold on
plot(mypredictionR,flyTruePrefR,'p','LineWidth',3,'Markersize',10)
plot(myprediction,flyTruePref,'d','LineWidth',3,'Markersize',10)
legend('Left Lobe Only','Right Lobe Only','Both Lobes')
legend boxoff
% for i=1:flyNum
%    hold on
%    text(myprediction(i)+0.01,flyTruePref(i),num2str(i),'FontSize',15)
% end
xlabel('Predicted Preference')
ylabel('Measured Preference')
set(gca,'FontSize',15)
box off
linmodelL
linmodelR

%% use  left and right lobe as separate data points

pcstouse=[2];

behaviorprediction=(SCORE(:,pcstouse));
flyTruePrefL=zeros(1,flyNum);
flyTruePrefR=zeros(1,flyNum);

ally=behaviorOcc';

nactivityL=zeros(flyNum,length(pcstouse));
nactivityR=zeros(flyNum,length(pcstouse));
for i=1:flyNum
    flyTruePrefL(i)=mean(ally(flyindicesL{i}));
    nactivityL(i,:)=mean(behaviorprediction(flyindicesL{i},:));
    flyTruePrefR(i)=mean(ally(flyindicesR{i}));
    nactivityR(i,:)=mean(behaviorprediction(flyindicesR{i},:));
end
linmodel=fitlm([nactivityL; nactivityR],[flyTruePrefL flyTruePrefR]);

myprediction=predict(linmodel,[nactivityL; nactivityR]);

figure %13
plot(myprediction,[flyTruePrefL flyTruePrefR],'*','LineWidth',3,'Markersize',10)

% for i=1:flyNum
%    hold on
%    text(myprediction(i)+0.01,flyTruePref(i),num2str(i),'FontSize',15)
% end
xlabel('Predicted Preference')
ylabel('Measured Preference')
set(gca,'FontSize',15)
box off
linmodel


%% use average DM2 - DC2 activity and relative oct/mch activation

behaviorprediction=100*(-mean(responsesNoResponseRemoved(1:13,:),1)+mean(responsesNoResponseRemoved(40:52,:),1))./(mean([mean(responsesNoResponseRemoved(1:13,:),1); mean(responsesNoResponseRemoved(40:52,:),1)]));

behaviorprediction=[behaviorprediction'];

flyTruePref=zeros(1,flyNum);
flyPredictedPref=zeros(1,flyNum);
ally=behaviorOcc';

linmodel=fitlm(behaviorprediction,ally);
myprediction=predict(linmodel,behaviorprediction);
figure %14
plot(myprediction,ally,'o','LineWidth',3)
xlabel('Predicted Preference')
ylabel('Measured Preference')
set(gca,'FontSize',15)
box off
linmodel

nactivity=zeros(flyNum,size(behaviorprediction,2));
for i=1:flyNum
   flyTruePref(i)=mean(ally(flyindices{i}));
   flyPredictedPref(i)=mean(mean(behaviorprediction(flyindices{i},:)));
   nactivity(i,:)=mean(behaviorprediction(flyindices{i},:));
end

% plot histogram of DM2 - DC2
% SUP FIG 13e
figure %15
histogram(nactivity,10)
ylabel('# flies')
xlabel('DM2 - DC2 (% df/f difference)')
axis square

% plot raw values
% SUP FIG 13f
figure %16
plot(nactivity,flyTruePref,'.','Color',pcolor, 'LineWidth',3,'MarkerSize',20)
xlabel('DM2 - DC2 (% df/f difference)')
ylabel('measured preference')
axis square

linmodel=fitlm(nactivity,flyTruePref);
myprediction=predict(linmodel,nactivity);
figure %17
plot(myprediction,flyTruePref,'ko','LineWidth',3)
for i=1:flyNum
   hold on
   %text(myprediction(i)+0.01,flyTruePref(i),num2str(i),'FontSize',15)
end
xlabel('Predicted Preference')
ylabel('Measured Preference')
set(gca,'FontSize',15)
box off
linmodel


% FIG 2i
figure %18
hold on;
xVals = (myprediction-mean(myprediction))/(std(myprediction));
yVals = (flyTruePref-mean(flyTruePref))/(std(flyTruePref));
linreg = linearRegressionCI2(xVals, yVals.', 1, 0, -3, 3);

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
axis([-3 3 -3 3])
xlabel('predicted preference (z-scored)')
ylabel('measured preference (z-scored)')
set(gca,'xtick','')
set(gca,'ytick','')
axis square

