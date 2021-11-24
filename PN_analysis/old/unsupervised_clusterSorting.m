%% use unsupervised algorithm to sort glomeruli across individual trials
clear all 
close all


manualLabelHome='/Users/mattchurgin/Dropbox/flyimaging/analysis/PN_analysis/manuallyLabelled_allThrough181210';


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
        behaviorpreOcc=[behaviorpreOcc preocc];

        gs=median(grnResponse(:,:,odortimes),3); % use median
        responseTemp=NaN*zeros(length(publishedOR.gh146glomerulusNames),nodors+3);
        responseTempTimeCourse=NaN*zeros(length(publishedOR.gh146glomerulusNames)*length(odortimesfortimecourse),size(grnResponse,1));
        responsesTensorTemp=NaN*zeros(length(publishedOR.gh146glomerulusNames),size(grnResponse,1),length(odortimesfortimecourse));
        centroids=[];
        for j=1:size(gs,2)
            coords=clusterInfoU{j}.Centroid;
            centroids=[centroids; coords];
        end
        % normalize centroids
        mc=mean(centroids);
        sc=std(centroids);
        centroids=(centroids-mc)./sc;
        
        % fill in response matrix
        for j=1:size(gs,2)
            
            k=j;    
            coords=centroids(j,:);
            
            responseTemp(k,:)=[(gs(:,j)-mean(gs,2))./transpose(std(gs')); coords']; % add x y z location
            glomsannotated(i)=glomsannotated(i)+1;
            glomfound(k)=glomfound(k)+1;
            for oo=1:nodors
                responseTempTimeCourse((((k-1)*length(odortimesfortimecourse)+1):k*length(odortimesfortimecourse)),oo)=grnResponse(oo,j,odortimesfortimecourse);
                responsesTensorTemp(k,oo,:)=grnResponse(oo,j,odortimesfortimecourse);
            end
        end
        
        responseTempT=responseTemp';
        
        responsesGlomByOdor=[responsesGlomByOdor responseTempT(:)];
       
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


%% reshape responseglombyodor
allgloms=responsesGlomByOdor(:);
allgloms=reshape(allgloms,nodors+3,length(allgloms)/(nodors+3));

% remove air and ethanol
allgloms([1 8],:)=0;


% run kmeans on allgloms
km=39;
[IDX, C, SUMD, D] = kmeans(allgloms', km);

% fill in new matrix with sorted gloms
unsupervisedSortedGlomResponse=NaN*zeros(nodors*39,size(allgloms,2)/39);
for i=1:size(allgloms,2)/39
    for j=1:39
        try
        currglom=IDX((i-1)*39+j);
        for k=1:nodors
            if ~isnan(unsupervisedSortedGlomResponse((currglom-1)*nodors+k,i)) % if multiple clusters map to same glomerulus, average reuslts
                unsupervisedSortedGlomResponse((currglom-1)*nodors+k,i)=(unsupervisedSortedGlomResponse((currglom-1)*nodors+k,i)+allgloms(k,(i-1)*39+j))/2;
            else
            unsupervisedSortedGlomResponse((currglom-1)*nodors+k,i)=allgloms(k,(i-1)*39+j);
            end
        end
        catch
        end
    end
end

%% perform pca on responses

clear responsesNoResponseRemoved

fracIn=0.05; % best results when fracIn is high, ~0.5, only using high confidence glomeruli


medianResponseOrTimeCourse=1; % 1 for median response only, 0 for full time course

if medianResponseOrTimeCourse
    responsesNoResponseRemoved=unsupervisedSortedGlomResponse;
else
    responsesNoResponseRemoved=responsesTimeCourseGlomByOdor;
end

% blank odors (when solenoids failed for example)
% if medianResponseOrTimeCourse
%     odorstoremove=[2 6 9 10 12];
%     blankThresh=[0.3 100 100 0.3 0.3]; % blank odor if 75% prctile response is below this (for odors that sporadically worked)
%     
%     % hard coded all flies before fly 23.  1:88 !!! Hard coded!!
%     for i=1:length(odorstoremove)
%         flyodorstoblank=find(prctile(responsesNoResponseRemoved(odorstoremove(i):nodors:end,:),75)<=blankThresh(i));
%         for j=1:(size(responsesNoResponseRemoved,1)/nodors)
%             temp=find(isfinite(responsesNoResponseRemoved(((j-1)*nodors+odorstoremove(i)),[1:88])));
%             tofill=intersect(flyodorstoblank,temp);  
%             if isfinite(nanmean(responsesNoResponseRemoved(((j-1)*nodors+odorstoremove(i)),89:end)))
%                 responsesNoResponseRemoved(((j-1)*nodors+odorstoremove(i)),tofill)=nanmean(responsesNoResponseRemoved(((j-1)*nodors+odorstoremove(i)),89:end));
%             end
%         end
%     end
%     
%     % remove MCH for flies on 181108
%     for j=1:(size(responsesNoResponseRemoved,1)/nodors)
%         tofill=find(isfinite(responsesNoResponseRemoved(((j-1)*nodors+11),89:100)))+88;
%         if isfinite(nanmean(responsesNoResponseRemoved(((j-1)*nodors+11),[1:89 101:end])))
%         responsesNoResponseRemoved(((j-1)*nodors+11),tofill)=nanmean(responsesNoResponseRemoved(((j-1)*nodors+11),[1:89 101:end]));
%         end
%     end
%     
% %  odorstoremove=[1 8]; % delete air and ethanol
% %  for i=1:length(odorstoremove)
% %      for j=1:nodors
% %      responsesNoResponseRemoved(((j-1)*nodors+odorstoremove(i)),:)=[];
% %      end
% %  end
% else
%     
% end

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

% remove air and ethanol
responsesNoResponseRemoved(1:nodors:end,:)=0;
responsesNoResponseRemoved(8:nodors:end,:)=0;

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
varianceToKeep=20; % percent of variance to keep


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


%% find predictive power of each pc
warning('off')

iters=10;

trainsize=20;

RON=zeros(size(co,2),iters);
PON=zeros(size(co,2),iters);
ROFF=zeros(size(co,2),iters);
POFF=zeros(size(co,2),iters);

totalCorr=zeros(1,size(co,2));
totalCorrS=zeros(1,size(co,2));
totalP=zeros(1,size(co,2));
totalPs=zeros(1,size(co,2));

for pcstouse=1:size(co,2)
    predictor=(co(:,pcstouse));
    predictor=predictor';
    
    fitCorr=zeros(1,iters);
    fitCorrs=zeros(1,iters);
    
    fitP=zeros(1,iters);
    fitPs=zeros(1,iters);
    
    for ii=1:iters
        temp=randperm(flyNum);
        trainflies=temp(1:floor(length(temp)*trainsize/100));
        validationflies=temp(ceil(length(temp)*trainsize/100):end);
        
        flyTruePref=zeros(1,length(trainflies));
        flyPredictedPref=zeros(1,length(trainflies));
        flyTruePrefShuffled=zeros(1,length(trainflies));
        for i=1:length(trainflies)
            flyTruePref(i)=mean(behaviorOcc(flyindices{trainflies(i)}));
            flyPredictedPref(i)=mean(predictor(1,flyindices{trainflies(i)}));
            
            flyTruePrefShuffled(i)=mean(behaviorOcc(flyindicesShuffled{trainflies(i)}));
        end
        
        validationflyTruePref=zeros(1,length(validationflies));
        validationflyPredictedPref=zeros(1,length(validationflies));
        validationflyTruePrefShuffled=zeros(1,length(validationflies));
        for i=1:length(validationflies)
            validationflyTruePref(i)=mean(behaviorOcc(flyindices{validationflies(i)}));
            validationflyPredictedPref(i)=mean(predictor(1,flyindices{validationflies(i)}));
            
            validationflyTruePrefShuffled(i)=mean(behaviorOcc(flyindicesShuffled{validationflies(i)}));
        end
        
        % generate linear model
        linModel=fitlm(flyPredictedPref,flyTruePref);
        % evaluate predictions for training data
        neuralPrediction=predict(linModel,flyPredictedPref');
        % evaluate predictions for test data
        neuralPredictionValidation=predict(linModel,validationflyPredictedPref');
        
        [r p]=corrcoef(neuralPredictionValidation,validationflyTruePref);
        %[r p]=corrcoef(neuralPrediction,flyTruePref);
        fitCorr(ii)=r(1,2);
        fitP(ii)=p(1,2);
        
        % generate linear model with shuffled data
        linModelShuffled=fitlm(flyPredictedPref,flyTruePrefShuffled);
        % evaluate predictions for shuffled training data
        neuralPredictionShuffled=predict(linModelShuffled,flyPredictedPref');
        % evaluate predictions for shuffled test data
        neuralPredictionValidationShuffled=predict(linModelShuffled,validationflyPredictedPref');
        
        [rs ps]=corrcoef(neuralPredictionValidationShuffled,validationflyTruePrefShuffled);
        %[rs ps]=corrcoef(neuralPredictionShuffled,flyTruePrefShuffled);
        fitCorrs(ii)=rs(1,2);
        fitPs(ii)=ps(1,2);
    end
    
    totalCorr(pcstouse)=mean(fitCorr);
    totalCorrS(pcstouse)=mean(fitCorrs);
    totalP(pcstouse)=mean(fitP);
    totalPs(pcstouse)=mean(fitPs);
    
    if mod(pcstouse,1)==0
        disp(['calculated fit using ' num2str(pcstouse) ' PCs of ' num2str(size(co,2))])
    end
end

left_color=[0 0 0];
right_color=[0.5 0 0.9];

fig=figure;
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(totalCorr,'o-','LineWidth',2)
ylabel('correlation')
hold on
yyaxis right
ylabel('p-value')
plot(log10(totalP),'-','LineWidth',2)
xlabel('PC #')
set(gca,'FontSize',15)

fig=figure;
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(totalCorrS,'o-','LineWidth',2)
ylabel('correlation')
hold on
yyaxis right
ylabel('p-value')
plot(log10(totalPs),'-','LineWidth',2)
xlabel('# PCs')
set(gca,'FontSize',15)

figure
subplot(2,1,1)
hist(totalCorr,20)
subplot(2,1,2)
hist(totalCorrS,20)
%% use fitlm

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
plot(myprediction,flyTruePref,'o','LineWidth',3)
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
ylabel('PC 2 loadings')
box off
set(gca,'FontSize',15)

currpc=COEFF(:,pcstouse);


