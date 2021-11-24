% Compare PN and ORN data
% use PN "test set" flies acquired during similar time period as ORNs

% load PNs
clear all
close all

pncolor=[0 0.8 0];
orncolor=[0.9 0 0.7];
rng('default')

manualLabelHome='/Users/mattchurgin/Dropbox/flyimaging/analysis/ORNvsPN_analysis/pnsadded181218plus2019';

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

% perform pca on PN responses

clear responsesNoResponseRemoved

fracIn=0.4;% best results when fracIn is high, ~0.5, only using high confidence glomeruli

responsesNoResponseRemoved=responsesGlomByOdor;




% remove fully empty rows
%totallyempty=~any(isfinite(responsesNoResponseRemoved),2);
%responsesNoResponseRemoved(totallyempty,:)=[];

gNames=publishedOR.gh146glomerulusNames;
glomsFound=glomfound;
numFinite=sum(isfinite(responsesNoResponseRemoved),2);
toRemove=find(numFinite/size(responsesNoResponseRemoved,2)<=fracIn);
responsesNoResponseRemoved(toRemove,:)=[];

temp=nodors;
fp=toRemove(find(mod(toRemove,temp)==1));
glomsremoved=((fp-1)/temp)+1;
gNames(glomsremoved)=[];
glomsFound(glomsremoved)=[];


% manually remove D
gNames(1)=[];
responsesNoResponseRemoved(1:nodors,:)=[];

gh146rawdata=responsesNoResponseRemoved;

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
% load ORNs

manualLabelHome='/Users/mattchurgin/Dropbox/flyimaging/analysis/ORNvsPN_analysis/ornflies';

publishedOdorPath='/Users/mattchurgin/Desktop/dblab/mattFunctions/odorpanelResponsesDoOR/odorPanel_12/odorPanel_12_DoORData.mat';
load(publishedOdorPath);
manualLabelledFolders=dir(manualLabelHome);
manualLabelledFolders=manualLabelledFolders(3:end);

labels=cell(1,length(manualLabelledFolders));
glomsannotatedO=zeros(1,length(manualLabelledFolders));
totalPutativeGloms=zeros(1,length(manualLabelledFolders));

nodors=13;
odortimes=[6:9]; % hard-coded specific time interval for summarizing odor response
odortimesfortimecourse=[1:18];
responsesGlomByOdorO=[];
responsesTimeCourseGlomByOdor=[];
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
        
        %behaviorOcc=[behaviorOcc occ_zeroed];
        %behaviorpreOcc=[behaviorpreOcc preocc_zeroed];
        behaviorOccO=[behaviorOccO occ-preocc];
        %behaviorOcc=[behaviorOcc occ];
        behaviorpreOccO=[behaviorpreOccO preocc];
        
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
                    glomsannotatedO(i)=glomsannotatedO(i)+1;
                    glomfoundO(k)=glomfoundO(k)+1;
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
        
        responsesGlomByOdorO=[responsesGlomByOdorO responseTempT(:)];
        responsesTimeCourseGlomByOdor=[responsesTimeCourseGlomByOdor responseTempTimeCourseT(:)];
        glombyodorfliesO=[glombyodorfliesO (flyNumO+(cL-1)*leftRightSplitVal)];
        
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


% remove fully empty rows
%totallyempty=~any(isfinite(responsesNoResponseRemoved),2);
%responsesNoResponseRemoved(totallyempty,:)=[];

gNamesO=publishedOR.gh146glomerulusNames;
glomsFoundO=glomfoundO;
numFinite=sum(isfinite(responsesNoResponseRemovedO),2);
toRemove=find(numFinite/size(responsesNoResponseRemovedO,2)<=fracIn);
responsesNoResponseRemovedO(toRemove,:)=[];

temp=nodors;
fp=toRemove(find(mod(toRemove,temp)==1));
glomsremoved=((fp-1)/temp)+1;
gNamesO(glomsremoved)=[];
glomsFoundO(glomsremoved)=[];

orcoRawData=responsesNoResponseRemovedO;

% % fill nans with mean
for i=1:size(responsesNoResponseRemovedO,1)
    for j=1:size(responsesNoResponseRemovedO,2)
        if isnan(responsesNoResponseRemovedO(i,j))
            responsesNoResponseRemovedO(i,j)=nanmean(responsesNoResponseRemovedO(i,:));
        end
    end
end

% remove air and ethanol
%responsesNoResponseRemoved(1:nodors:end,:)=0;
%responsesNoResponseRemoved(8:nodors:end,:)=0;

% fill missing values with linear interpolation
% data=responsesNoResponseRemoved';
% dz=data;
% dataFilled=fillWithRegressedValues(dz);
% responsesNoResponseRemoved=dataFilled';

yesZscore=0;
if yesZscore
    responsesNoResponseRemovedO=responsesNoResponseRemovedO';
    responsesNoResponseRemovedO=(responsesNoResponseRemovedO-mean(responsesNoResponseRemovedO))./std(responsesNoResponseRemovedO);
    responsesNoResponseRemovedO=responsesNoResponseRemovedO';
end

glombyodorfliesTemp=glombyodorflies;
glombyodorfliesOTemp=glombyodorfliesO;

%% keep only one glomerulus
figure(1)
figure(2)
for glomToKeep = 1:5
    glombyodorflies=glombyodorfliesTemp;
    glombyodorfliesO=glombyodorfliesOTemp;
    
    gh146SingleGlom=gh146rawdata(((glomToKeep-1)*nodors)+(1:(nodors)),:);
    orcoSingleGlom=orcoRawData(((glomToKeep-1)*nodors)+(1:(nodors)),:);
    
    deleteTrialsWithMissingData=1; % 1: delete trials with missing data
    
    if deleteTrialsWithMissingData
        ghToDelete=isnan(mean(gh146SingleGlom));
        orcoToDelete=isnan(mean(orcoSingleGlom));
        gh146SingleGlom(:,ghToDelete)=[];
        orcoSingleGlom(:,orcoToDelete)=[];
        glombyodorflies(ghToDelete)=[];
        glombyodorfliesO(orcoToDelete)=[];
    end
    
    responsesNoResponseRemoved=gh146SingleGlom;
    responsesNoResponseRemovedO=orcoSingleGlom;
    
    % calculate PCs
    
    % PNs
    
    opt = statset('pca');
    opt.Display='iter';
    [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(responsesNoResponseRemoved','Options',opt);
    
    
    % ORNs
    opt = statset('pca');
    opt.Display='iter';
    [COEFFO, SCOREO, LATENTO, TSQUAREDO, EXPLAINEDO] = pca(responsesNoResponseRemovedO','Options',opt);
    
    
    % figure;
    % plot(EXPLAINED,'-','Color',pncolor,'LineWidth',3)
    % hold on
    % plot(EXPLAINEDO,'o-','Color',orncolor,'LineWidth',3)
    % ylabel('variance explained (%)')
    % xlabel('PC #')
    % legend('PNs','ORNs')
    % legend boxoff
    % box off
    % set(gca,'FontSize',20)
    %
    % figure;
    % semilogy(EXPLAINED,'-','Color',pncolor,'LineWidth',3)
    % hold on
    % semilogy(EXPLAINEDO,'o-','Color',orncolor,'LineWidth',3)
    % ylabel('variance explained (%)')
    % xlabel('PC #')
    % legend('PNs','ORNs')
    % legend boxoff
    % box off
    % set(gca,'FontSize',20)
    
    % plot pc loadings
    % figure;
    % for i=1:8
    %     subplot(2,4,i)
    %     plot(COEFF(:,i),'*','Color',pncolor,'LineWidth',2,'MarkerSize',8)
    %     hold on
    %     plot(zeros(1,length(COEFF(:,i))),'k--','LineWidth',3)
    %     j=1;
    %     for ii=1:nodors:length(COEFF)
    %         plot((ii-0.5)*ones(1,5), linspace(min(COEFF(:,i)),max(COEFF(:,i)),5),'k--','LineWidth',2)
    %         %text(i+floor(nodors/3),min(PCContribution),num2str(glomsFound(j)),'FontSize',15)
    %         j=j+1;
    %     end
    %
    %     set(gca,'xtick',(1:nodors:length(COEFF(:,i)))+floor(nodors/2),'xticklabel',string(gNames),'FontSize',10)
    %     xtickangle(30)
    %     ylabel(['PC ' num2str(i) ' loadings'])
    %     box off
    %     set(gca,'FontSize',15)
    % end
    %
    % figure;
    % for i=1:8
    %     subplot(2,4,i)
    %     plot(COEFFO(:,i),'*','Color',orncolor,'LineWidth',2,'MarkerSize',8)
    %     hold on
    %     plot(zeros(1,length(COEFFO(:,i))),'k--','LineWidth',3)
    %     j=1;
    %     for ii=1:nodors:length(COEFFO)
    %         plot((ii-0.5)*ones(1,5), linspace(min(COEFFO(:,i)),max(COEFFO(:,i)),5),'k--','LineWidth',2)
    %         %text(i+floor(nodors/3),min(PCContribution),num2str(glomsFound(j)),'FontSize',15)
    %         j=j+1;
    %     end
    %
    %     set(gca,'xtick',(1:nodors:length(COEFFO(:,i)))+floor(nodors/2),'xticklabel',string(gNamesO),'FontSize',10)
    %     xtickangle(30)
    %     ylabel(['PC ' num2str(i) ' loadings'])
    %     box off
    %     set(gca,'FontSize',15)
    % end
    
    % calculate PN distances
    
    varianceToKeep=100; % percent of variance to keep
    
    
    co = SCORE;
    totalVarianceExplained=cumsum(EXPLAINED);
    pcsWithinVariance=find(totalVarianceExplained<varianceToKeep);
    pcstouse=pcsWithinVariance(end);
    %disp(['Using ' num2str(pcstouse) ' PCs'])
    
    
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
    % figure
    % hold on
    %
    % % plot(0,0,'Marker',m{1},'Color','k')
    % % plot(0,0,'Marker',m{2},'Color','k')
    % %legend('Left Lobe','Right Lobe')
    % for j=1:(flyNum)
    %
    %     ltemps=find(glombyodorflies==(j));
    %     lcurr=ltemps(1:end);
    %     rtemps=find(glombyodorflies==(j+1000));
    %     rcurr=rtemps(1:end);
    %
    %     if length(lcurr)>0
    %         plot(co(lcurr(1:size(lcurr,2)),1),co(lcurr(1:size(lcurr,2)),2),'Color',mycmap(j,:),'Marker','*','LineWidth',lsize,'MarkerSize',msize);
    %
    %     end
    %     if length(rcurr)>0
    %         plot(co(rcurr(1:size(rcurr,2)),1),co(rcurr(1:size(rcurr,2)),2),'Color',mycmap(j,:),'Marker','o','LineWidth',lsize,'MarkerSize',msize);
    %     end
    %
    %     if length(lcurr)>0  && length(rcurr)>0
    %         plot([co(lcurr(1:size(lcurr,2)),1)' co(rcurr(1:size(rcurr,2)),1)' co(lcurr(1),1)],[co(lcurr(1:size(lcurr,2)),2)' co(rcurr(1:size(rcurr,2)),2)' co(lcurr(1),2)],'Color',mycmap(j,:),'LineStyle','--','LineWidth',lsize,'MarkerSize',msize);
    %     end
    %
    % end
    % xlabel('PC 1 Score')
    % ylabel('PC 2 Score')
    % set(gca,'FontSize',15)
    %
    % figure
    % hold on
    % for j=1:(flyNum)
    %
    %     ltemps=find(glombyodorflies==(j));
    %     lcurr=ltemps(1:end);
    %     rtemps=find(glombyodorflies==(j+1000));
    %     rcurr=rtemps(1:end);
    %
    %
    %     if length(lcurr)>0
    %         h3=plot(mean(co(lcurr(1:size(lcurr,2)),1)),mean(co(lcurr(1:size(lcurr,2)),2)),'Color',mycmap(j,:),'Marker','*','LineWidth',lsize,'MarkerSize',msize);
    %         %text(mean(co(lcurr(1:size(lcurr,2)),1)),mean(co(lcurr(1:size(lcurr,2)),2)),[num2str(j)],'FontSize',15)
    %     end
    %     if  length(rcurr)>0
    %         h4=plot(mean(co(rcurr(1:size(rcurr,2)),1)),mean(co(rcurr(1:size(rcurr,2)),2)),'Color',mycmap(j,:),'Marker','o','LineWidth',lsize,'MarkerSize',msize);
    %     end
    %     if  length(lcurr)>0 && length(rcurr)>0
    %         h5=plot([mean(co(lcurr(1:size(lcurr,2)),1))' mean(co(rcurr(1:size(rcurr,2)),1))'],[mean(co(lcurr(1:size(lcurr,2)),2))' mean(co(rcurr(1:size(rcurr,2)),2))'],'Color',mycmap(j,:),'LineStyle','--','LineWidth',lsize,'MarkerSize',1);
    %     end
    %
    % end
    % h1=plot(-100,-100,'Marker','o','MarkerSize',msize,'LineStyle','none','LineWidth',lsize,'Color','k');
    % h2=plot(-100,-100,'Marker','*','MarkerSize',msize,'LineStyle','none','LineWidth',lsize,'Color','k');
    % legend([h1 h2],{'Left Lobe','Right Lobe'})
    % legend boxoff
    % box off
    % axis([min(co(:,1)) max(co(:,1)) min(co(:,2)) max(co(:,2))])
    % xlabel('PC 1 Score')
    % ylabel('PC 2 Score')
    % set(gca,'FontSize',15)
    
    
    
    %Boxplots of distances
    % summarize across flies
    
    withinleft=withinLLobe;
    withinright=withinRLobe;
    withinacross=withinDifferentLobe;
    acrossleft=acrossLLobe;
    acrossright=acrossRLobe;
    acrossall=acrossAllLobe;
    
    % figure
    % boxplot([withinleft(:) withinright(:) withinacross(:) acrossleft(:) acrossright(:) acrossall(:)])
    % ylabel('Distance in Coding Space')
    % xlabels{1}='Within Fly (Left Lobe)';
    % xlabels{2}='Within Fly (Right Lobe)';
    % xlabels{3}='Within Fly (Opposite Lobes)';
    %
    % xlabels{4}='Across Fly (Left Lobe)';
    % xlabels{5}='Across Fly (Right Lobe)';
    % xlabels{6}='Across Fly (Both Lobes)';
    % set(gca,'xtick',1:6,'xticklabel',xlabels,'FontSize',10)
    % xtickangle(30)
    % box off
    % set(gca,'FontSize',15)
    
    % calculate ORN distances
    
    co = SCOREO;
    totalVarianceExplained=cumsum(EXPLAINEDO);
    pcsWithinVariance=find(totalVarianceExplained<varianceToKeep);
    pcstouse=pcsWithinVariance(end);
    %disp(['Using ' num2str(pcstouse) ' PCs'])
    
    
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
    %
    % % plot average across all odors for flies with both lobe data
    % figure
    % hold on
    %
    % % plot(0,0,'Marker',m{1},'Color','k')
    % % plot(0,0,'Marker',m{2},'Color','k')
    % %legend('Left Lobe','Right Lobe')
    % for j=1:(flyNumO)
    %
    %     ltemps=find(glombyodorfliesO==(j));
    %     lcurr=ltemps(1:end);
    %     rtemps=find(glombyodorfliesO==(j+1000));
    %     rcurr=rtemps(1:end);
    %
    %     if length(lcurr)>0
    %         plot(co(lcurr(1:size(lcurr,2)),1),co(lcurr(1:size(lcurr,2)),2),'Color',mycmapO(j,:),'Marker','*','LineWidth',lsize,'MarkerSize',msize);
    %
    %     end
    %     if length(rcurr)>0
    %         plot(co(rcurr(1:size(rcurr,2)),1),co(rcurr(1:size(rcurr,2)),2),'Color',mycmapO(j,:),'Marker','o','LineWidth',lsize,'MarkerSize',msize);
    %     end
    %
    %     if length(lcurr)>0  && length(rcurr)>0
    %         plot([co(lcurr(1:size(lcurr,2)),1)' co(rcurr(1:size(rcurr,2)),1)' co(lcurr(1),1)],[co(lcurr(1:size(lcurr,2)),2)' co(rcurr(1:size(rcurr,2)),2)' co(lcurr(1),2)],'Color',mycmapO(j,:),'LineStyle','--','LineWidth',lsize,'MarkerSize',msize);
    %     end
    %
    % end
    % xlabel('PC 1 Score')
    % ylabel('PC 2 Score')
    % set(gca,'FontSize',15)
    %
    % figure
    % hold on
    % for j=1:(flyNumO)
    %
    %     ltemps=find(glombyodorfliesO==(j));
    %     lcurr=ltemps(1:end);
    %     rtemps=find(glombyodorfliesO==(j+1000));
    %     rcurr=rtemps(1:end);
    %
    %
    %     if length(lcurr)>0
    %         h3=plot(mean(co(lcurr(1:size(lcurr,2)),1)),mean(co(lcurr(1:size(lcurr,2)),2)),'Color',mycmapO(j,:),'Marker','*','LineWidth',lsize,'MarkerSize',msize);
    %         %text(mean(co(lcurr(1:size(lcurr,2)),1)),mean(co(lcurr(1:size(lcurr,2)),2)),[num2str(j)],'FontSize',15)
    %     end
    %     if  length(rcurr)>0
    %         h4=plot(mean(co(rcurr(1:size(rcurr,2)),1)),mean(co(rcurr(1:size(rcurr,2)),2)),'Color',mycmapO(j,:),'Marker','o','LineWidth',lsize,'MarkerSize',msize);
    %     end
    %     if  length(lcurr)>0 && length(rcurr)>0
    %         h5=plot([mean(co(lcurr(1:size(lcurr,2)),1))' mean(co(rcurr(1:size(rcurr,2)),1))'],[mean(co(lcurr(1:size(lcurr,2)),2))' mean(co(rcurr(1:size(rcurr,2)),2))'],'Color',mycmapO(j,:),'LineStyle','--','LineWidth',lsize,'MarkerSize',1);
    %     end
    %
    % end
    % h1=plot(-100,-100,'Marker','o','MarkerSize',msize,'LineStyle','none','LineWidth',lsize,'Color','k');
    % h2=plot(-100,-100,'Marker','*','MarkerSize',msize,'LineStyle','none','LineWidth',lsize,'Color','k');
    % legend([h1 h2],{'Left Lobe','Right Lobe'})
    % legend boxoff
    % box off
    % axis([min(co(:,1)) max(co(:,1)) min(co(:,2)) max(co(:,2))])
    % xlabel('PC 1 Score')
    % ylabel('PC 2 Score')
    % set(gca,'FontSize',15)
    
    
    %Boxplots of distances
    % summarize across flies
    %
    withinleftO=withinLLobe;
    withinrightO=withinRLobe;
    withinacrossO=withinDifferentLobe;
    acrossleftO=acrossLLobe;
    acrossrightO=acrossRLobe;
    acrossallO=acrossAllLobe;
    
    % figure
    % boxplot([withinleftO(:) withinrightO(:) withinacrossO(:) acrossleftO(:) acrossrightO(:) acrossallO(:)])
    % ylabel('Distance in Coding Space')
    % xlabels{1}='Within Fly (Left Lobe)';
    % xlabels{2}='Within Fly (Right Lobe)';
    % xlabels{3}='Within Fly (Opposite Lobes)';
    %
    % xlabels{4}='Across Fly (Left Lobe)';
    % xlabels{5}='Across Fly (Right Lobe)';
    % xlabels{6}='Across Fly (Both Lobes)';
    % set(gca,'xtick',1:6,'xticklabel',xlabels,'FontSize',10)
    % xtickangle(30)
    % box off
    % set(gca,'FontSize',15)
    
    %
    
    
    % compare within versus across fly
    PNtoboxplot=[withinleft(:) withinright(:) withinacross(:) acrossleft(:) acrossright(:) acrossall(:)];
    ORNtoboxplot=[ withinleftO(:) withinrightO(:) withinacrossO(:) acrossleftO(:) acrossrightO(:) acrossallO(:)];
    
    combined=[PNtoboxplot(:); ORNtoboxplot(:)];
    groupings=[];
    groupings=[groupings 1*ones(1,3*flyNum)];
    groupings=[groupings 3*ones(1,3*flyNum)];
    
    
    groupings=[groupings 2*ones(1,3*flyNumO)];
    groupings=[groupings 4*ones(1,3*flyNumO)];
    
    [aa PNp cc]=ttest2(combined(groupings==1),combined(groupings==3));
    [aa ORNp cc]=ttest2(combined(groupings==2),combined(groupings==4));
    
    colors=[0 0.8 0; 0.9 0 0.7; 0 0.8 0; 0.9 0 0.7];
    figure(1)
    subplot(1,5,glomToKeep)
    boxplot(combined,groupings,'PlotStyle','compact','Colors',colors)
    %boxplot(combined,groupings,'BoxStyle','filled','Colors',colors)
    if glomToKeep ==1
        ylabel('distance')
    else
        set(gca,'ytick','')
    end
%     text(1,2,['PN p = ' num2str(PNp)],'FontSize',10)
%     text(1,1,['ORN p = ' num2str(ORNp)],'FontSize',10)
    xlabels{1}='within';
    xlabels{2}='across';
    set(gca,'xtick',[1.5 3.5],'xticklabel',xlabels,'FontSize',10)
    xtickangle(30)
    axis([0 4.5 0 3])
    box off
    title([gNames{glomToKeep}])
    set(gca,'FontSize',15)
    
    PNwithin=[withinleft(:) withinright(:) withinacross(:)];
    PNacross=[acrossleft(:) acrossright(:) acrossall(:)];
    PNind=(nanmean(PNacross,2)-(nanmean(PNwithin,2)))./(nanmean(PNacross,2)+nanmean(PNwithin,2));
    
    ORNwithin=[withinleftO(:) withinrightO(:) withinacrossO(:)];
    ORNacross=[acrossleftO(:) acrossrightO(:) acrossallO(:)];
    ORNind=(nanmean(ORNacross,2)-(nanmean(ORNwithin,2)))./(nanmean(ORNacross,2)+nanmean(ORNwithin,2));
    
    figure(2)
    subplot(1,5,glomToKeep)
     boxplot([PNind(:); ORNind(:)],[1*ones(1,length(PNind)) 2*ones(1,length(ORNind))],'PlotStyle','compact','Colors',colors)
    %boxplot(combined,groupings,'BoxStyle','filled','Colors',colors)
    if glomToKeep ==1
        ylabel('relative distance (across-within)')
    else
        set(gca,'ytick','')
    end
    [aa diffP cc]=ttest2(PNind(:),ORNind(:));
     text(1,0,['p = ' num2str(diffP)],'FontSize',10)
    set(gca,'xtick',[1.5 3.5],'FontSize',10)
    xtickangle(30)
    axis([0 2.5 -0.5 1])
    box off
    title([gNames{glomToKeep}])
    set(gca,'FontSize',15)
end

