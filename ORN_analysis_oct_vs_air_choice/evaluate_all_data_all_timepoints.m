% look at manually labelled gh146 flies in odor space
clear all
%close all
load ORN_PN_colors
load analysis_dir_path
manualLabelHome=fullfile(analysis_dir_path, 'ORN_analysis_oct_vs_air_choice/ornflies_oct_vs_air');

publishedOdorPath=fullfile(analysis_dir_path, 'odorPanel_12_DoORData.mat');
load(publishedOdorPath);

manualLabelledFolders=dir(manualLabelHome);
manualLabelledFolders=manualLabelledFolders(3:end);

labels=cell(1,length(manualLabelledFolders));
glomsannotated=zeros(1,length(manualLabelledFolders));
totalPutativeGloms=zeros(1,length(manualLabelledFolders));

nodors=13;
odortimes=[6:9]; % hard-coded specific time interval for summarizing odor response
odortimesfortimecourse=[1:20];
tc = length(odortimesfortimecourse);
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
        behaviorpreOcc=[behaviorpreOcc preocc];
        
        manualClusterLabels=clusterLabels;
        totalPutativeGloms(i)=length(manualClusterLabels);
        
        
        gs=median(grnResponse(:,:,odortimes),3); % use median
        %gs=prctile(grnResponse(:,:,odortimes),75); % use percentile
        responseTemp=NaN*zeros(length(publishedOR.gh146glomerulusNames),nodors);
        responseTempTimeCourse=NaN*zeros(length(publishedOR.gh146glomerulusNames)*length(odortimesfortimecourse),size(grnResponse,1));
        responseTempTimeCourse=NaN*zeros(length(publishedOR.gh146glomerulusNames)*nodors,length(odortimesfortimecourse));
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
                        %responseTempTimeCourse((((k-1)*length(odortimesfortimecourse)+1):k*length(odortimesfortimecourse)),oo)=grnResponse(oo,j,odortimesfortimecourse);
                        responsesTensorTemp(k,oo,:)=grnResponse(oo,j,odortimesfortimecourse);
                    end
                    for oo=1:length(odortimesfortimecourse)
                        responseTempTimeCourse((((k-1)*nodors+1):k*nodors),oo)=grnResponse(:,j,odortimesfortimecourse(oo))-nanmedian(grnResponse(:,j,1:5),3);
                        %responsesTensorTemp(k,oo,:)=grnResponse(oo,j,odortimesfortimecourse);
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


medianResponseOrTimeCourse=0; % 1 for median response only, 0 for full time course

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



% % fill nans with mean
for i=1:size(responsesNoResponseRemoved,1)
    for j=1:size(responsesNoResponseRemoved,2)
        if isnan(responsesNoResponseRemoved(i,j))
            responsesNoResponseRemoved(i,j)=nanmean(responsesNoResponseRemoved(i,:));
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

%% build linear model using individual pc

pcstouse=[1];

behaviorprediction=(SCORE(:,pcstouse));
flyTruePref=zeros(1,flyNum);
flyPredictedPref=zeros(1,flyNum);
ally=behaviorOcc';

linmodel=fitlm(behaviorprediction,ally);
myprediction=predict(linmodel,behaviorprediction);
figure %2
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
figure %3
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

PCContribution=COEFF(:,pcstouse);
figure; %4
plot(PCContribution,'*','LineWidth',2,'MarkerSize',8)
hold on
plot(zeros(1,length(PCContribution(:,1))),'k--','LineWidth',3)
j=1;
for i=1:(nodors*tc):length(PCContribution)
    plot((i-0.5)*ones(1,5), linspace(min(PCContribution),max(PCContribution),5),'k--','LineWidth',3)
    %text(i+floor(nodors/3),min(PCContribution),num2str(glomsFound(j)),'FontSize',15)
    j=j+1;
end
for i=1:(tc):length(PCContribution)
    plot((i-0.5)*ones(1,5), linspace(min(PCContribution),max(PCContribution),5),'m--','LineWidth',1)
    %text(i+floor(nodors/3),min(PCContribution),num2str(glomsFound(j)),'FontSize',15)
    j=j+1;
end

set(gca,'xtick',(1:(nodors*tc):length(PCContribution))+floor((nodors*tc)/2),'xticklabel',string(gNames),'FontSize',10)
xtickangle(30)
ylabel('PC loadings')
box off
set(gca,'FontSize',15)

currpc=COEFF(:,pcstouse);


mycolors=parula(tc);
figure; %5
for i =1:tc
plot(PCContribution(i:tc:end),'*','Color',mycolors(i,:),'LineWidth',2,'MarkerSize',8)
hold on
end
j=1;
for i=1:(nodors):length(PCContribution(1:tc:end))
    plot((i-0.5)*ones(1,5), linspace(min(PCContribution),max(PCContribution),5),'k--','LineWidth',3)
    j=j+1;
end
set(gca,'xtick',(1:(nodors):length(PCContribution(1:tc:end)))+floor((nodors)/2),'xticklabel',string(gNames),'FontSize',10)
xtickangle(30)
ylabel('PC loadings')
box off
set(gca,'FontSize',15)

%% decode behavior using PCs corresponding to each time point individually
pcstouse=[1];
for j = 1:tc
    
    tempCoeff=zeros(size(COEFF));
    tempCoeff(j:tc:end,:)=COEFF(j:tc:end,:);
    
    behaviorprediction=(responsesNoResponseRemoved')*tempCoeff(:,pcstouse);
    flyTruePref=zeros(1,flyNum);
    flyPredictedPref=zeros(1,flyNum);
    ally=behaviorOcc';
    
    
    nactivity=zeros(flyNum,length(pcstouse));
    for i=1:flyNum
        flyTruePref(i)=mean(ally(flyindices{i}));
        flyPredictedPref(i)=mean(mean(behaviorprediction(flyindices{i},:)));
        nactivity(i,:)=mean(behaviorprediction(flyindices{i},:));
    end
    linmodel=fitlm(nactivity,flyTruePref);
    r2(j) = linmodel.Rsquared.Ordinary;
    
end

figure; %6
plot(1.2:1.2:25,r2,'-','Color',ocolor,'LineWidth',3,'MarkerSize',15)
hold on
plot(6*ones(5,1),linspace(0,0.18,5),'k--')
text(6.2,0.18,'odor on','FontSize',15)
plot(10*ones(5,1),linspace(0,0.18,5),'k--')
text(10.2,0.18,'odor off','FontSize',15)
xlabel('time (s)')
ylabel('behavior decoding (R^2)')
axis([0 25 0 0.3])
set(gca,'FontSize',15)

%% build cross temporal decoder

pcstouse=[1];
for j = 1:tc
    
    tempCoeff=zeros(size(COEFF));
    tempCoeff(j:tc:end,:)=COEFF(j:tc:end,:);
    
    
    behaviorprediction=(responsesNoResponseRemoved')*tempCoeff(:,pcstouse);
    flyTruePref=zeros(1,flyNum);
    flyPredictedPref=zeros(1,flyNum);
    ally=behaviorOcc';
    
    
    nactivity=zeros(flyNum,length(pcstouse));
    for i=1:flyNum
        flyTruePref(i)=mean(ally(flyindices{i}));
        nactivity(i,:)=mean(behaviorprediction(flyindices{i},:));
    end
    linmodel=fitlm(nactivity,flyTruePref);

    for k = 1:tc
        tempCoeff2=zeros(size(COEFF));
        tempCoeff2(k:tc:end,:)=COEFF(j:tc:end,:);
        behaviorprediction2=(responsesNoResponseRemoved')*tempCoeff2(:,pcstouse);
        nactivity=zeros(flyNum,length(pcstouse));
        for i=1:flyNum
            flyTruePref(i)=mean(ally(flyindices{i}));
            nactivity(i,:)=mean(behaviorprediction2(flyindices{i},:));
        end
        preds = predict(linmodel,nactivity);
        ctd(j,k)=(corr(preds,flyTruePref')^2)*sign(corr(preds,flyTruePref'));
    end
end

figure %7
imagesc(1.2:1.2:25,1.2:1.2:25,ctd,[-0.025 0.2])
ylabel('train time (s)')
xlabel('test time (s)')
hcb=colorbar;
title(hcb,'R^2')
set(gca,'FontSize',15)



%% perform pca on each time point alone
clear all

load ORN_PN_colors
load analysis_dir_path
manualLabelHome=fullfile(analysis_dir_path, 'ORN_analysis_oct_vs_air_choice/ornflies_oct_vs_air');

publishedOdorPath=fullfile(analysis_dir_path, 'odorPanel_12_DoORData.mat');
load(publishedOdorPath);

manualLabelledFolders=dir(manualLabelHome);
manualLabelledFolders=manualLabelledFolders(3:end);

labels=cell(1,length(manualLabelledFolders));
glomsannotated=zeros(1,length(manualLabelledFolders));
totalPutativeGloms=zeros(1,length(manualLabelledFolders));

for odortimes = 1:20
    disp(['loading time point ' num2str(odortimes)])
    nodors=13;
    
    odortimesfortimecourse=[1:20];
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
            behaviorpreOcc=[behaviorpreOcc preocc];
            
            manualClusterLabels=clusterLabels;
            totalPutativeGloms(i)=length(manualClusterLabels);
            
            gs=zeros(13,size(grnResponse,2));
            for j=1:size(grnResponse,2)
                temp = squeeze(grnResponse(:,j,odortimes)-nanmedian(grnResponse(:,j,1:5),3));
                gs(:,j)= transpose((temp'));
            end
            
            responseTemp=NaN*zeros(length(publishedOR.gh146glomerulusNames),nodors);
            responseTempTimeCourse=NaN*zeros(length(publishedOR.gh146glomerulusNames)*length(odortimesfortimecourse),size(grnResponse,1));
            responseTempTimeCourse=NaN*zeros(length(publishedOR.gh146glomerulusNames)*nodors,length(odortimesfortimecourse));
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
                            %responseTempTimeCourse((((k-1)*length(odortimesfortimecourse)+1):k*length(odortimesfortimecourse)),oo)=grnResponse(oo,j,odortimesfortimecourse);
                            responsesTensorTemp(k,oo,:)=grnResponse(oo,j,odortimesfortimecourse);
                        end
                        for oo=1:length(odortimesfortimecourse)
                            responseTempTimeCourse((((k-1)*nodors+1):k*nodors),oo)=grnResponse(:,j,odortimesfortimecourse(oo))-nanmedian(grnResponse(:,j,1:5),3);
                            %responsesTensorTemp(k,oo,:)=grnResponse(oo,j,odortimesfortimecourse);
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
    
    % perform pca on responses
    
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
   
    
    
    % % fill nans with mean
    for i=1:size(responsesNoResponseRemoved,1)
        for j=1:size(responsesNoResponseRemoved,2)
            if isnan(responsesNoResponseRemoved(i,j))
                responsesNoResponseRemoved(i,j)=nanmean(responsesNoResponseRemoved(i,:));
            end
        end
    end
    
    responsesNoResponseRemoved=responsesNoResponseRemoved';
    responsesNoResponseRemoved=(responsesNoResponseRemoved-mean(responsesNoResponseRemoved));
    responsesNoResponseRemoved=responsesNoResponseRemoved';
    
    opt = statset('pca');
    opt.Display='iter';
    [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(responsesNoResponseRemoved','Options',opt);
    
    
    pcstouse=[1];
    
    behaviorprediction=(SCORE(:,pcstouse));
    flyTruePref=zeros(1,flyNum);
    flyPredictedPref=zeros(1,flyNum);
    ally=behaviorOcc';
    
    nactivity=zeros(flyNum,length(pcstouse));
    for i=1:flyNum
        flyTruePref(i)=mean(ally(flyindices{i}));
        flyPredictedPref(i)=mean(mean(behaviorprediction(flyindices{i},:)));
        nactivity(i,:)=mean(behaviorprediction(flyindices{i},:));
    end
    linmodel=fitlm(nactivity,flyTruePref);
    myprediction=predict(linmodel,nactivity);
    
    myr2(odortimes)=linmodel.Rsquared.Ordinary;
end


% SUP FIG 10b
figure %8
plot(1:20,myr2,'Color',ocolor,'LineWidth',3)
xlabel('time')
ylabel('R^2')
%axis([0 20 0 0.25])
