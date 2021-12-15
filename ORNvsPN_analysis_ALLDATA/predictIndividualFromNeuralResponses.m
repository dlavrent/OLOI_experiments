% Predict individual identity from  neural activity

% load PNs
clear all
close all

load ORN_PN_COLORS
allColorMaps

rng('default')

load analysis_dir_path

manualLabelHome=fullfile(analysis_dir_path, 'ORNvsPN_analysis_ALLDATA/pn_alldata');

publishedOdorPath=fullfile(analysis_dir_path, 'utilities/odorPanel_12_DoORData.mat');
load(publishedOdorPath);

manualLabelledFolders=dir(manualLabelHome);
manualLabelledFolders=manualLabelledFolders(3:end);

labels=cell(1,length(manualLabelledFolders));
glomsannotated=zeros(1,length(manualLabelledFolders));
totalPutativeGloms=zeros(1,length(manualLabelledFolders));
conmat=zeros(13,13,20);
pihatcg=cell(1,20);
trsplit=0.5;
fv=30; % fraction of variance to retain
fvs=[10 30 50 60 80];

nodors=13;
odortimes=[6:9];
odortimesfortimecourse=[1:18];
%responsesTensorFullyOrthogonal=NaN*zeros(22,2,39,nodors,length(odortimesfortimecourse)); % flynum HARD CODED
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
            gs(:,j)= transpose(max(temp'));
        end
        responseTemp=NaN*zeros(length(publishedOR.gh146glomerulusNames),nodors);
        %             responseTempTimeCourse=NaN*zeros(length(publishedOR.gh146glomerulusNames)*length(odortimesfortimecourse),size(grnResponse,1));
        %             responsesTensorTemp=NaN*zeros(length(publishedOR.gh146glomerulusNames),size(grnResponse,1),length(odortimesfortimecourse));
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
                    %                         for oo=1:nodors
                    %                             responseTempTimeCourse((((k-1)*length(odortimesfortimecourse)+1):k*length(odortimesfortimecourse)),oo)=grnResponse(oo,j,odortimesfortimecourse);
                    %                             responsesTensorTemp(k,oo,:)=grnResponse(oo,j,odortimesfortimecourse)-nanmedian(grnResponse(oo,j,1:5),3);
                    %                         end
                    break
                end
            end
        end
        
        responseTempT=responseTemp';
        %responseTempTimeCourseT=responseTempTimeCourse';
        
        %responsesTensorFullyOrthogonal(flyNum,cL,:,:,:)=responsesTensorTemp;
        
        responsesGlomByOdor=[responsesGlomByOdor responseTempT(:)];
        %responsesTimeCourseGlomByOdor=[responsesTimeCourseGlomByOdor responseTempTimeCourseT(:)];
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
glombyodorflies(glombyodorflies>1000)=glombyodorflies(glombyodorflies>1000)-1000;

mycmap=distinguishable_colors(flyNum);

% perform pca on PN responses

clear responsesNoResponseRemoved

fracIn=0.4;% best results when fracIn is high, ~0.5, only using high confidence glomeruli

responsesNoResponseRemoved=responsesGlomByOdor;
%responsesNoResponseRemovedT=responsesTimeCourseGlomByOdor;
%responsesTensor=responsesTensorFullyOrthogonal;

% blank odors (when solenoids failed for example)

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

gNames=publishedOR.gh146glomerulusNames;
glomsFound=glomfound;
numFinite=sum(isfinite(responsesNoResponseRemoved),2);
toRemove=find(numFinite/size(responsesNoResponseRemoved,2)<=fracIn);
responsesNoResponseRemoved(toRemove,:)=[];
%numFiniteT=sum(isfinite(responsesNoResponseRemovedT),2);
%toRemoveT=find(numFiniteT/size(responsesNoResponseRemovedT,2)<=fracIn);
%responsesNoResponseRemovedT(toRemoveT,:)=[];

temp=nodors;
fp=toRemove(find(mod(toRemove,temp)==1));
glomsremoved=((fp-1)/temp)+1;
gNames(glomsremoved)=[];
glomsFound(glomsremoved)=[];
%responsesTensor(:,:,glomsremoved,:,:)=[];


% manually remove D
gNames(1)=[];
responsesNoResponseRemoved(1:nodors,:)=[];
%responsesNoResponseRemovedT(1:((nodors)*length(odortimesfortimecourse)),:)=[];
%responsesTensor(:,:,1,:,:)=[];

gh146rawdata=responsesNoResponseRemoved;

% % fill nans with mean
for i=1:size(responsesNoResponseRemoved,1)
    for j=1:size(responsesNoResponseRemoved,2)
        if isnan(responsesNoResponseRemoved(i,j))
            responsesNoResponseRemoved(i,j)=nanmean(responsesNoResponseRemoved(i,:));
        end
    end
end

g=responsesNoResponseRemoved;
g=g./mean(g);

temp = randperm(length(glombyodorflies));
sp = round(trsplit*length(temp));
tr = temp(1:sp); te = temp((sp+1):end);

tr = temp;
te = temp;
disp(['fitting logistic regression. gh146'])


[coe sco as asd exp]=pca(g');

jj=0;
for fv=fvs
    tic
    disp(['gh146 fv = ' num2str(fv)])
    jj=jj+1;
    inds = find(cumsum(exp)>fv);
    toppc=inds(1);
    actualfvg(jj)=sum(exp(1:toppc));
     
    [B]=mnrfit(sco(tr,1:toppc),categorical((transpose(glombyodorflies(tr)))));
    pihat = mnrval(B,sco(te,1:toppc));
    [val in] = max(pihat,[],2);
    classact = sum(transpose(in) == glombyodorflies(te))/length(glombyodorflies(te));
    gh146acc(jj)=classact;
    
    [B]=mnrfit(sco(tr,1:toppc),categorical((transpose(glombyodorflies(tr(randperm(length(tr))))))));
    pihat = mnrval(B,sco(te,1:toppc));
    [val in] = max(pihat,[],2);
    classact = sum(transpose(in) == glombyodorflies(te(randperm(length(tr)))))/length(glombyodorflies(te));
    gh146acc_shuffled(jj)=classact;
    disp(['gh146 fv = ' num2str(fv) ' time elapsed = ' num2str(toc/60,'%2.1f') ' minutes'])
end




% now for orn
rng('default')


manualLabelHome=fullfile(analysis_dir_path, 'ORNvsPN_analysis_ALLDATA/orn_alldata');

publishedOdorPath=fullfile(analysis_dir_path, 'utilities/odorPanel_12_DoORData.mat');
load(publishedOdorPath);

manualLabelledFolders=dir(manualLabelHome);
manualLabelledFolders=manualLabelledFolders(3:end);

labels=cell(1,length(manualLabelledFolders));
glomsannotated=zeros(1,length(manualLabelledFolders));
totalPutativeGloms=zeros(1,length(manualLabelledFolders));
conmato=zeros(13,13,20);
pihatco=cell(1,20);

odortimes=[6:9];
nodors=13;
odortimesfortimecourse=[1:18];
%responsesTensorFullyOrthogonal=NaN*zeros(22,2,39,nodors,length(odortimesfortimecourse)); % flynum HARD CODED
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
            gs(:,j)= transpose(max(temp'));
        end
        
        responseTemp=NaN*zeros(length(publishedOR.gh146glomerulusNames),nodors);
        %             responseTempTimeCourse=NaN*zeros(length(publishedOR.gh146glomerulusNames)*length(odortimesfortimecourse),size(grnResponse,1));
        %             responsesTensorTemp=NaN*zeros(length(publishedOR.gh146glomerulusNames),size(grnResponse,1),length(odortimesfortimecourse));
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
                    %                         for oo=1:nodors
                    %                             responseTempTimeCourse((((k-1)*length(odortimesfortimecourse)+1):k*length(odortimesfortimecourse)),oo)=grnResponse(oo,j,odortimesfortimecourse);
                    %                             responsesTensorTemp(k,oo,:)=grnResponse(oo,j,odortimesfortimecourse)-nanmedian(grnResponse(oo,j,1:5),3);
                    %                         end
                    break
                end
            end
        end
        
        responseTempT=responseTemp';
        %responseTempTimeCourseT=responseTempTimeCourse';
        
        %responsesTensorFullyOrthogonal(flyNum,cL,:,:,:)=responsesTensorTemp;
        
        responsesGlomByOdor=[responsesGlomByOdor responseTempT(:)];
        %responsesTimeCourseGlomByOdor=[responsesTimeCourseGlomByOdor responseTempTimeCourseT(:)];
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
glombyodorflies(glombyodorflies>1000)=glombyodorflies(glombyodorflies>1000)-1000;
mycmap=distinguishable_colors(flyNum);

clear responsesNoResponseRemoved

fracIn=0.25;

responsesNoResponseRemoved=responsesGlomByOdor;
%responsesNoResponseRemovedT=responsesTimeCourseGlomByOdor;
%responsesTensor=responsesTensorFullyOrthogonal;

gNames=publishedOR.gh146glomerulusNames;
glomsFound=glomfound;
numFinite=sum(isfinite(responsesNoResponseRemoved),2);
toRemove=find(numFinite/size(responsesNoResponseRemoved,2)<=fracIn);
responsesNoResponseRemoved(toRemove,:)=[];
%numFiniteT=sum(isfinite(responsesNoResponseRemovedT),2);
%toRemoveT=find(numFiniteT/size(responsesNoResponseRemovedT,2)<=fracIn);
%responsesNoResponseRemovedT(toRemoveT,:)=[];

temp=nodors;
fp=toRemove(find(mod(toRemove,temp)==1));
glomsremoved=((fp-1)/temp)+1;
gNames(glomsremoved)=[];
glomsFound(glomsremoved)=[];
%responsesTensor(:,:,glomsremoved,:,:)=[];


% % fill nans with mean
for i=1:size(responsesNoResponseRemoved,1)
    for j=1:size(responsesNoResponseRemoved,2)
        if isnan(responsesNoResponseRemoved(i,j))
            responsesNoResponseRemoved(i,j)=nanmean(responsesNoResponseRemoved(i,:));
        end
    end
end

disp(['fitting logistic regression. orco'])
g=responsesNoResponseRemoved;
g=g./mean(g);

temp = randperm(length(glombyodorflies));
sp = round(trsplit*length(temp));
tr = temp(1:sp); te = temp((sp+1):end);

tr = temp;
te = temp;

[coe sco as asd exp]=pca(g');

jj=0;
for fv=fvs
    tic
    jj=jj+1;
    disp(['orco fv = ' num2str(fv)])
    inds = find(cumsum(exp)>fv);
    toppc=inds(1);
    actualfvo(jj)=sum(exp(1:toppc));
    
    [B]=mnrfit(sco(tr,1:toppc),categorical((transpose(glombyodorflies(tr)))));
    pihat = mnrval(B,sco(te,1:toppc));
    [val in] = max(pihat,[],2);
    classact = sum(transpose(in) == glombyodorflies(te))/length(glombyodorflies(te));
    orcoacc(jj)=classact;
    
    [B]=mnrfit(sco(tr,1:toppc),categorical((transpose(glombyodorflies(tr(randperm(length(tr))))))));
    pihat = mnrval(B,sco(te,1:toppc));
    [val in] = max(pihat,[],2);
    classact = sum(transpose(in) == glombyodorflies(te(randperm(length(tr)))))/length(glombyodorflies(te));
    orcoacc_shuffled(jj)=classact;
    
    disp(['orco fv = ' num2str(fv) ' time elapsed = ' num2str(toc/60,'%2.1f') ' minutes'])
end

%%
figure %1
plot(actualfvo,orcoacc,'x-','Color',ocolor,'LineWidth',3,'MarkerSize',10)
hold on
plot(actualfvo,orcoacc_shuffled,'x--','Color',ocolor,'LineWidth',3,'MarkerSize',10)
plot(actualfvg,gh146acc,'o-','Color',pcolor,'LineWidth',3,'MarkerSize',10)
plot(actualfvg,gh146acc_shuffled,'o--','Color',pcolor,'LineWidth',3,'MarkerSize',10)
legend('ORN','ORN shuffled','PN','PN shuffled')
legend boxoff
xlabel('% variance retained')
ylabel('individual identity decoding accuracy')
axis([10 90 0 1])
axis square
set(gca,'FontSize',15)

figure %2
plot(actualfvo,orcoacc,'x','Color',ocolor,'LineWidth',3,'MarkerSize',10)
hold on
plot(actualfvo,orcoacc_shuffled,'x--','Color',ocolor,'LineWidth',3,'MarkerSize',10)
plot(actualfvg,gh146acc,'o','Color',pcolor,'LineWidth',3,'MarkerSize',10)
plot(actualfvg,gh146acc_shuffled,'o--','Color',pcolor,'LineWidth',3,'MarkerSize',10)
legend('ORN','ORN shuffled','PN','PN shuffled')
legend boxoff
xlabel('% variance retained')
ylabel('individual identity decoding accuracy')
axis([10 90 0 1])
set(gca,'FontSize',15)
%%
% plot
x=[1:4];
y=[orcoacc orcoacc_shuffled gh146acc gh146acc_shuffled];
figure %3
bar(x,y,'FaceColor','k')