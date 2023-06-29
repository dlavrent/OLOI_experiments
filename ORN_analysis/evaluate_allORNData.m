% look at manually labelled gh146 flies in odor space
clear all
%close all

load ORN_PN_colors

load analysis_dir_path

manualLabelHome=fullfile(analysis_dir_path, 'ORN_analysis/ornflies');
trainedModel = load(fullfile(analysis_dir_path, 'ORN_analysis/trainDataModel.mat'));
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
        %behaviorOcc=[behaviorOcc occ-preocc];
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
        %gs=prctile(grnResponse(:,:,odortimes),90,3); % use percentile
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

cell_strs = cell(1, flyNum);
for i = 1:flyNum
    ls = flyindicesL{i};
    rs = flyindicesR{i};
    curstr = sprintf('L%dR%d', length(ls), length(rs));
    cell_strs{i} = curstr;
end 
CHAR_CELL_STRs = categorical(cell_strs);

curstr_categories = categories(CHAR_CELL_STRs);
num_in_categories = countcats(CHAR_CELL_STRs);
ncats = length(curstr_categories);
FLYNUM_SUMMARY = cell(2, ncats);
FLYNUM_SUMMARY(1, :) = curstr_categories;
for csi = 1:ncats
    FLYNUM_SUMMARY{2, csi} = num_in_categories(csi);
end 
FLYNUM_SUMMARY = FLYNUM_SUMMARY'

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

fracIn=0.25; % best results when fracIn is high, ~0.5, only using high confidence glomeruli


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

responsesNoResponseRemoved=responsesNoResponseRemoved';
responsesNoResponseRemoved=(responsesNoResponseRemoved-mean(responsesNoResponseRemoved));
responsesNoResponseRemoved=responsesNoResponseRemoved';

opt = statset('pca');
opt.Display='iter';
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(responsesNoResponseRemoved','Options',opt);


figure; %1
plot(log10(EXPLAINED),'o-','LineWidth',3)
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

num_L2R2 = num_in_categories(ncats);
disp(sprintf('num L2R2: %d', num_L2R2));
nactivity=zeros(flyNum,length(pcstouse));
lrtri1s = zeros(1, num_L2R2);
lrtri2s = zeros(1, num_L2R2);
lrcntr = 1;
for i=1:flyNum
    flyTruePref(i)=mean(ally(flyindices{i}));
    flyPredictedPref(i)=mean(mean(behaviorprediction(flyindices{i},:)));
    nactivity(i,:)=mean(behaviorprediction(flyindices{i},:));
    
    mycurls = behaviorprediction(flyindicesL{i});
    mycurrs = behaviorprediction(flyindicesR{i});
    if (length(mycurls) == 2 ) && (length(mycurrs) == 2)       
        lrtri1s(lrcntr) = (mycurls(1) + mycurrs(1))/2;
        lrtri2s(lrcntr) = (mycurls(2) + mycurrs(2))/2;
        lrcntr = lrcntr + 1;
    end
end
[rlrtri, plrtri] = corrcoef(lrtri1s, lrtri2s);
disp(sprintf('r for %d mean L+R tria1 1 to 2 scores: %f', num_L2R2, rlrtri(1,2)));
linmodel=fitlm(nactivity,flyTruePref);
myprediction=predict(linmodel,nactivity);
figure %7
plot(myprediction,flyTruePref,'*','Color',ocolor,'LineWidth',3)
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
currpc=COEFF(:,pcstouse);

PCContribution=COEFF(:,pcstouse);
figure; %8
plot(PCContribution,'*','Color',ocolor,'LineWidth',2,'MarkerSize',8)
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

% FIG measured pref vs predicted OCT-MCH ORN PC 1
figure %9
hold on;
xVals = (myprediction-mean(myprediction))/(std(myprediction));
yVals = (flyTruePref-mean(flyTruePref))/(std(flyTruePref));
linreg = linearRegressionCI2(xVals, yVals.', 1, 0, -3, 2);

areaBar(linreg.xVals,polyval(linreg.pOverall,linreg.xVals),2*std(linreg.fits),[0 0 0],[0.9 0.9 0.9])
plot(xVals,yVals,'.','Color',ocolor, 'LineWidth',3)
for i=1:flyNum
   hold on
   %text(myprediction(i)+0.01,flyTruePref(i),num2str(i),'FontSize',15)
end
[r p]=corrcoef(xVals,yVals);
text(-.25,-.25,['r = ' num2str(r(1,2),'%2.2f')],'FontSize',15)
text(-.3,-.25,['p = ' num2str(p(1,2),'%2.2f')],'FontSize',15)
set(gca,'FontSize',15)
box on
axis([-3 2 -3 2])
axis square
xlabel('predicted preference (z-scored)')
ylabel('measured preference (z-scored)')
set(gca,'xtick','')
set(gca,'ytick','')

linmodel
% 
% % get odor valences
% for i=1:nodors
%     odorValence(i)=mean(PCContribution(i:nodors:end));
% end
% %  only use when using one PC for prediction
% figure
% hold on
% for i=1:length(gNames)
%     plot(currpc(2+13*(i-1)),currpc(11+13*(i-1)),'o','LineWidth',2)
%     text(1.05*currpc(2+13*(i-1)),1.05*currpc(11+13*(i-1)),gNames{i},'FontSize',15)
% end
% xlabel('PC 5 OCT loading')
% ylabel('PC 5 MCH loading')
% box off
% set(gca,'FontSize',15)

%%

matchedResponses=responsesNoResponseRemoved;
COEFFmatched=COEFF;
trainedPC=trainedModel.mypc;
trainedGNames=trainedModel.gNames;
% match available gloms in test data with gloms from trainedPC
testGlomInTrainGlom=zeros(1,length(gNames));
for i=1:length(gNames)
    for j=1:length(trainedGNames)
        if strcmp(gNames{i},trainedGNames{j})
            testGlomInTrainGlom(i)=1;
            break
        end
    end
end
% remove data for test gloms not in train set
todeleteTest=[];
for i=1:length(gNames)
   if ~testGlomInTrainGlom(i)
       todeleteTest=[todeleteTest ((nodors*(i-1)+1):(nodors*i))];
   end
end
matchedResponses(todeleteTest,:)=[];


trainGlomInTestGlom=zeros(1,length(trainedGNames));
for i=1:length(trainedGNames)
    for j=1:length(gNames)
        if strcmp(gNames{j},trainedGNames{i})
            trainGlomInTestGlom(i)=1;
            break
        end
    end
end
% remove data for train gloms PC not in test set
todeleteTrain=[];
for i=1:length(trainedGNames)
   if ~trainGlomInTestGlom(i)
       todeleteTrain=[todeleteTrain ((nodors*(i-1)+1):(nodors*i))];
   end
end
trainedPC(todeleteTrain)=[];

COEFFmatched(todeleteTest,:)=[];
% measure correlation between the trained pc and test data PC
pccorr=zeros(1,10);
for i=1:10
    tempr=corr(trainedPC,COEFFmatched(:,i));
    pccorr(i)=tempr;
end

behaviorprediction=matchedResponses'*trainedPC;

flyTruePref=zeros(1,flyNum);
flyPredictedPref=zeros(1,flyNum);
ally=behaviorOcc';

nactivity=zeros(flyNum,1);
for i=1:flyNum
    flyTruePref(i)=mean(ally(flyindices{i}));
    flyPredictedPref(i)=mean(mean(behaviorprediction(flyindices{i},:)));
    nactivity(i,:)=mean(behaviorprediction(flyindices{i},:));
end
linmodel=trainedModel.linmodelPrecorrected;
myprediction=predict(linmodel,nactivity);
figure %10
plot(myprediction,flyTruePref,'.','Color',ocolor,'LineWidth',3)
hold on
xlabel('Predicted Preference')
ylabel('Measured Preference')
box on
axis([-.6 .3 -.6 .3])
axis square
set(gca,'FontSize',15)


corrcoef(myprediction,flyTruePref)

figure %11
plot((myprediction-mean(myprediction))/(std(myprediction)),(flyTruePref-mean(flyTruePref))/(std(flyTruePref)),'.','Color',ocolor, 'LineWidth',3,'MarkerSize',15)
for i=1:flyNum
   hold on
   %text(myprediction(i)+0.01,flyTruePref(i),num2str(i),'FontSize',15)
end
[r p]=corrcoef((myprediction-mean(myprediction))/(std(myprediction)),(flyTruePref-mean(flyTruePref))/(std(flyTruePref)));
text(-.25,-.25,['r = ' num2str(r(1,2),'%2.2f')],'FontSize',15)
text(-.3,-.25,['p = ' num2str(p(1,2),'%2.2f')],'FontSize',15)
set(gca,'FontSize',15)
box on
axis([-2.8 2 -2.8 2])
axis square
xlabel('predicted preference (z-scored)')
ylabel('measured preference (z-scored)')
set(gca,'xtick','')
set(gca,'ytick','')


%% use average DC2 - DM2 activity and relative oct/mch activation
orncolor=[0.9 0 0.7];

behaviorprediction=(mean(responsesNoResponseRemoved(1:13,:),1)-mean(responsesNoResponseRemoved(40:52,:),1));
%behaviorprediction=(mean(responsesNoResponseRemoved(1:13,:),1)-mean(responsesNoResponseRemoved(40:52,:),1))./(mean(responsesNoResponseRemoved(1:13,:),1)+mean(responsesNoResponseRemoved(40:52,:),1));

behaviorprediction=[behaviorprediction'];

%uncomment to get ORN projection onto PN PC2
%pnmodel = load('PN_analysis\trainDataModel.mat');
%pnpc2 = pnmodel.mypc;
%ndat = length(ally);
%ornprojtopnpc2 = zeros(ndat, 1);
%for i=1:ndat
%    ornprojtopnpc2(i) = sum(pnpc2 .* responsesNoResponseRemoved(:, i));
%end
%behaviorprediction=ornprojtopnpc2;

flyTruePref=zeros(1,flyNum);
flyPredictedPref=zeros(1,flyNum);
ally=behaviorOcc';

linmodel=fitlm(behaviorprediction,ally);
myprediction=predict(linmodel,behaviorprediction);
figure %12
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
figure %13
plot(myprediction,flyTruePref,'o','Color',orncolor,'LineWidth',3)
for i=1:flyNum
   hold on
   text(myprediction(i)+0.01,flyTruePref(i),num2str(i),'FontSize',15)
end

xlabel('Predicted Preference')
ylabel('Measured Preference')
text(0,0,['R^2 = ' num2str(linmodel.Rsquared.Adjusted)],'FontSize',15)
set(gca,'FontSize',15)
box off
linmodel



