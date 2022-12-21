%% load data from preprocessed file
clear all
close all
load orcobrpshort_all_data
%% load data from raw file 
clear all
% close all

load analysis_dir_path

apple=0;
if apple==1
    startdir=fullfile(analysis_dir_path, 'IHC/OrcoBrpshort/autoSegmentation');
    manualSegDir=fullfile(analysis_dir_path, 'IHC/OrcoBrpshort/orcoOX5brpshortOXZ');
else
    startdir=fullfile(analysis_dir_path, 'IHC\OrcoBrpshort\autoSegmentation');
    manualSegDir=fullfile(analysis_dir_path, 'IHC\OrcoBrpshort\orcoOX5brpshortOXZ');
end
cd(startdir)

folders{1}='190227_orcoOX5_brpshort';
folders{2}='190404_orcobrpshort_behaviorAndImaging';
folders{3}='190516_orcoBrpshortOXZ';
folders{4}='190608_orcoOX5brpshort_behaviorandihc';


testfolders{1}='190710_orcobrpshort_ihc_age8days';
testfolders{2}='190717_orcobrpshort_age13days';
testfolders{3}='190724_orcobrpshort_ihc_age9days';
testfolders{4}='191119_orco_brpshort_behaviorAndIHC_age_day8';
sizeConversion=((.09*11)^2)*0.49;
densityOrTotal=1;

flyNum=0;
for i=1:length(folders)
    cd(folders{i})
    
    currdate=folders{i}(1:6);
    
    currfiles=dir(pwd);
    for j=1:length(currfiles)
        if strfind(currfiles(j).name,'labelled')
            load(currfiles(j).name)
            load(currfiles(j).name(1:(end-13)))
            us=strfind(currfiles(j).name,'_');
            currFlyName=currfiles(j).name(1:(us(1)-1));
            
            
            
            % load behavior
            load([manualSegDir '/' currdate '_' currFlyName '/' currFlyName '_behavior.mat'])
            
            
            
            flyNum=flyNum+1;
            
            odorb(flyNum)=occ-preocc;
            odorpre(flyNum)=preocc;
            
            
            for leftright=[1 2]
                for k=1:length(glomeruliToSegment)
                    glomSize(2*(flyNum-1)+leftright,k)=sizeConversion*sum(glomeruliSegs{leftright,k}(:));
                    glomTotalF(2*(flyNum-1)+leftright,k)=sum(b(:).*glomeruliSegs{leftright,k}(:));
                    glomRelativeF(2*(flyNum-1)+leftright,k)=glomTotalF(2*(flyNum-1)+leftright,k)/sum(glomeruliSegs{leftright,k}(:));
                end
            end
        
            disp(['loaded training fly ' num2str(flyNum)])
            
        end
    end
    
    cd ..
end


% load test flies

flyNumTest=0;
for i=1:length(testfolders)
    cd(testfolders{i})
    
    currdate=testfolders{i}(1:6);
    
    currfiles=dir(pwd);
    for j=1:length(currfiles)
        if strfind(currfiles(j).name,'labelled')
            load(currfiles(j).name)
            load(currfiles(j).name(1:(end-13)))
            us=strfind(currfiles(j).name,'_');
            currFlyName=currfiles(j).name(1:(us(1)-1));
            
            % load behavior
            load([manualSegDir '/' currdate '_' currFlyName '/' currFlyName '_behavior.mat'])
            
            flyNumTest=flyNumTest+1;
            odorbTest(flyNumTest)=occ-preocc;
            odorpreTest(flyNumTest)=preocc;
            
            
            for leftright=[1 2]
                for k=1:length(glomeruliToSegment)
                    glomSizeTest(2*(flyNumTest-1)+leftright,k)=sizeConversion*sum(glomeruliSegs{leftright,k}(:));
                    glomTotalFTest(2*(flyNumTest-1)+leftright,k)=sum(b(:).*glomeruliSegs{leftright,k}(:));
                    glomRelativeFTest(2*(flyNumTest-1)+leftright,k)=glomTotalFTest(2*(flyNumTest-1)+leftright,k)/sum(glomeruliSegs{leftright,k}(:));
                end
            end

            disp(['loaded test fly ' num2str(flyNumTest)])
            
        end
    end
    
    cd ..
end
disp('done loading')


glomSize(:,4)=[];
glomTotalF(:,4)=[];
glomRelativeF(:,4)=[];
glomSizeTest(:,4)=[];
glomTotalFTest(:,4)=[];
glomRelativeFTest(:,4)=[];
glomeruliToSegment(4)=[];

allGlomSize = [glomSize; glomSizeTest];
allGlomTotalF = [glomTotalF; glomTotalFTest];
allGlomRelativeF = [glomRelativeF; glomRelativeFTest];
allodorb = [odorb odorbTest];
allodorpre = [odorpre odorpreTest];

glomSizeOrig = glomSize;
glomTotalFOrig = glomTotalF;
glomRelativeFOrig = glomRelativeF;
odorbOrig = odorb;
odorpreOrig = odorpre;
%% set to plot only training data or all data
load ORN_PN_colors
trainingonly =0; % 0 to include train + test data, 1 to include only train data for model building
if trainingonly
    glomSize = glomSizeOrig;
    glomTotalF = glomTotalFOrig;
    glomRelativeF = glomRelativeFOrig;
    odorb = odorbOrig;
    odorpre = odorpreOrig;
else
    glomSize = allGlomSize;
    glomTotalF = allGlomTotalF;
    glomRelativeF = allGlomRelativeF;
    odorb = allodorb;
    odorpre = allodorpre;
end

%%  plot glomerulus properties

figure %1
subplot(1,3,1)
% SUP FIG 12a
plot(glomSize((1:flyNum)*2-1,:),glomSize((1:flyNum)*2,:),'.','LineWidth',2,'MarkerSize',15)
title('volume (um^3)')
xlabel('left lobe')
ylabel('right lobe')
set(gca,'FontSize',15)
axis square
axis([0 6200 0 6200])
subplot(1,3,2)
% SUP FIG 12b
plot(glomTotalF((1:flyNum)*2-1,:),glomTotalF((1:flyNum)*2,:),'.','LineWidth',2,'MarkerSize',15)
title('total fluorescence')
xlabel('left lobe')
ylabel('right lobe')
set(gca,'FontSize',15)
axis square
axis([0 5e5 0 5e5])
subplot(1,3,3)
% SUP FIG 12c
plot(glomRelativeF((1:flyNum)*2-1,:),glomRelativeF((1:flyNum)*2,:),'.','LineWidth',2,'MarkerSize',15)
legend(glomeruliToSegment)
legend boxoff
title('normalized fluorescence')
xlabel('left lobe')
ylabel('right lobe')
set(gca,'FontSize',15)
axis square
axis([10 80 10 80])

% figure 2
% SUP FIG 12d
violinPlot(glomSize,[ocolor])
set(gca,'XTickLabel',glomeruliToSegment)
xtickangle(30)
ylabel('volume (um^3)')
axis([0 5 0 1.2*max(max(glomSize))])
box on
set(gca,'FontSize',15)

% figure 3
% SUP FIG 12e
violinPlot(glomTotalF,[ocolor])
set(gca,'XTickLabel',glomeruliToSegment)
xtickangle(30)
ylabel('total fluorescence')
axis([0 5 0 1.2*max(max(glomTotalF))])
box on
set(gca,'FontSize',15)

% figure 4
% SUP FIG 12f
violinPlot(glomRelativeF,[ocolor])
set(gca,'XTickLabel',glomeruliToSegment)
xtickangle(30)
ylabel('fluorescence/volume')
axis([0 5 0 1.2*max(max(glomRelativeF))])
box on
set(gca,'FontSize',15)

%% perform pca on auto segmented data

normalize = 0; % 1 for normalize, 0 for don't
if ~normalize
    % standardize
    glomS=glomSize;
    glomF=glomTotalF;
    glomR=glomRelativeF;
else
    
    % normalize
    normGlomS =zeros(size(glomRelativeF,1),size(glomRelativeF,2));
    normGlomF =zeros(size(glomRelativeF,1),size(glomRelativeF,2));
    normGlomR =zeros(size(glomRelativeF,1),size(glomRelativeF,2));
    for i = 1:size(glomRelativeF,1)
        normGlomS(i,:)= glomSize(i,:)/sum(glomSize(i,:));
        normGlomF(i,:)= glomTotalF(i,:)/sum(glomTotalF(i,:));
        normGlomR(i,:)= glomRelativeF(i,:)/sum(glomRelativeF(i,:));
    end
    glomS=normGlomS;
    glomF=normGlomF;
    glomR=normGlomR;
end

[coeffS scoreS latentS tsqS explainedS] = pca(glomS); % volume

[coeffF scoreF latentF tsqF explainedF] = pca(glomF); % total F

[coeffR scoreR latentR tsqR explainedR] = pca(glomR); % relative F (F/volume)

figure %5
for i=1:size(coeffF,2)
    
    subplot(3,size(coeffF,2),i)
    plot(coeffS(:,i),'mo','LineWidth',3,'MarkerSize',10)
    hold on
    plot(0:(size(coeffS(:,i),1)+1),zeros(1,size(coeffS(:,i),1)+2),'k--','LineWidth',2)
    axis([0 size(coeffS(:,i),1)+1 -1 1])
    ylabel('loading')
    set(gca,'XTick',[1:size(coeffF,2)])
    set(gca,'XTickLabel',glomeruliToSegment)
    title(['PC ' num2str(i)])
    xtickangle(30)
    box off
    set(gca,'FontSize',15)
    
    
    subplot(3,size(coeffF,2),i+size(coeffF,2))
    plot(coeffF(:,i),'ro','LineWidth',3,'MarkerSize',10)
    hold on
    plot(0:(size(coeffF(:,i),1)+1),zeros(1,size(coeffF(:,i),1)+2),'k--','LineWidth',2)
    axis([0 size(coeffF(:,i),1)+1 -1 1])
    ylabel('loading')
    set(gca,'XTick',[1:size(coeffF,2)])
    set(gca,'XTickLabel',glomeruliToSegment)
    title(['PC ' num2str(i)])
    xtickangle(30)
    box off
    set(gca,'FontSize',15)
    
    
    subplot(3,size(coeffF,2),i+2*size(coeffF,2))
    plot(coeffR(:,i),'bo','LineWidth',3,'MarkerSize',10)
    hold on
    plot(0:(size(coeffR(:,i),1)+1),zeros(1,size(coeffR(:,i),1)+2),'k--','LineWidth',2)
    axis([0 size(coeffR(:,i),1)+1 -1 1])
    ylabel('loading')
    set(gca,'XTick',[1:size(coeffF,2)])
    set(gca,'XTickLabel',glomeruliToSegment)
    title(['PC ' num2str(i)])
    xtickangle(30)
    box off
    set(gca,'FontSize',15)
    
end

% SUP FIG 12g/h (toggle trainingonly = 0/1)
figure %6
for i=1:size(coeffF,2)
    
    subplot(1,size(coeffF,2),i)
    plot(coeffR(:,i),'.','Color',ocolor,'LineWidth',2,'MarkerSize',15)
    hold on
    plot(0:(size(coeffR(:,i),1)+1),zeros(1,size(coeffR(:,i),1)+2),'k--','LineWidth',2)
    axis([0 size(coeffR(:,i),1)+1 -1 1])
    ylabel('loading')
    set(gca,'XTick',[1:size(coeffF,2)])
    set(gca,'XTickLabel',glomeruliToSegment)
    title(['PC ' num2str(i)])
    xtickangle(30)
    box on
    set(gca,'FontSize',15)
    
end

% FIG 3e if trainingonly=1
figure %7
plot(coeffR(:,2),'.','Color',ocolor,'LineWidth',3,'MarkerSize',20)
hold on
plot(0:(size(coeffR(:,i),1)+1),zeros(1,size(coeffR(:,i),1)+2),'k--','LineWidth',2)
ylabel('PC 2 loadings')
set(gca,'XTick',[1:size(coeffF,2)])
set(gca,'XTickLabel',glomeruliToSegment)
xtickangle(30)
axis([0 5 -.6 .85])
set(gca,'FontSize',15)
%% plot principal component 2 vs behavior
% close all

componentR=2;
measureToUse=scoreR;

go=find(odorb);

if trainingonly
    pcmodel = fitlm(mean([measureToUse(2*(go-1)+1,componentR)'; measureToUse(2*go,componentR)']),odorb);
    save('IHC_training_model.mat', 'pcmodel');
else
    load('IHC_training_model.mat');
end

ypredtrain = predict(pcmodel,transpose(mean([measureToUse(2*go-1,componentR)'; measureToUse(2*go,componentR)'])));

myprediction=ypredtrain;
flyTruePref=odorb(go);
% FIG 3f if toggletrainingonly=1
figure %8
hold on;
xVals = (myprediction-mean(myprediction))/(std(myprediction));
yVals = (flyTruePref-mean(flyTruePref))/(std(flyTruePref));
linreg = linearRegressionCI2(xVals, yVals.', 1, 0, -2.2, 2.2);

areaBar(linreg.xVals,polyval(linreg.pOverall,linreg.xVals),2*std(linreg.fits),[0 0 0],[0.9 0.9 0.9])
plot(xVals,yVals,'.','Color',ocolor, 'LineWidth',3,'MarkerSize',15)

[r p]=corrcoef(xVals,yVals);
text(-.25,-.25,['r = ' num2str(r(1,2),'%2.2f')],'FontSize',15)
text(-.3,-.25,['p = ' num2str(p(1,2),'%2.2f')],'FontSize',15)
text(-.3,-.75,['n = ' num2str(size(xVals),'%2.2f')],'FontSize',15)
xlabel('predicted preference (z-scored)')
ylabel('measured preference (z-scored)')
set(gca,'FontSize',15)
box on
axis([-2.2 2.2 -2.2 2.2])
set(gca,'xtick','')
set(gca,'ytick','')
axis square


%% plot DM2 - DC2 versus behavior

dif = mean([glomR(2*(go-1)+1,1)'-glomR(2*(go-1)+1,3)'; glomR(2*(go),1)'-glomR(2*(go),3)']);
percentdif = 100*mean([(glomR(2*(go-1)+1,1)'-glomR(2*(go-1)+1,3)')./(mean([glomR(2*(go-1)+1,1)'; glomR(2*(go-1)+1,3)'])); (glomR(2*(go),1)'-glomR(2*(go),3)')./(mean([glomR(2*(go),1)'; glomR(2*(go),3)']))]);
% SUP FIG 13g
figure %9
histogram(percentdif,10)
xlabel('DM2 - DC2 (% Brp-Short density)')
ylabel('# flies')
axis square

% SUP FIG 13h
figure %10
plot(percentdif,odorb,'.','Color',ocolor, 'LineWidth',3,'MarkerSize',20)
xlabel('DM2 - DC2 (% Brp-Short density)')
ylabel('measured preference')
axis square

figure %11
hold on
for i=go
    plot(mean([(glomR(2*(i-1)+1,1)'-glomR(2*(i-1)+1,3)'); (glomR(2*(i),1)'-glomR(2*(i),3)')]),odorb(i),'mo','LineWidth',3)
end
[rd pd]=corrcoef(mean([glomR(2*(go-1)+1,1)'-glomR(2*(go-1)+1,3)'; glomR(2*(go),1)'-glomR(2*(go),3)']),odorb(go));
text(0.0,0.1,['r = ' num2str(rd(1,2))],'FontSize',15)
text(0.0,0.1,['p = ' num2str(pd(1,2))],'FontSize',15)
xlabel('DM2-DC2 Brpshort density')
ylabel('odor preference')
box off
set(gca,'FontSize',15)

calciummodel = fitlm(mean([glomR(2*(go-1)+1,1)'-glomR(2*(go-1)+1,3)'; glomR(2*(go),1)'-glomR(2*(go),3)']),odorb);

calciumpredtrain = predict(calciummodel,transpose(mean([glomR(2*(go-1)+1,1)'-glomR(2*(go-1)+1,3)'; glomR(2*(go),1)'-glomR(2*(go),3)'])));

myprediction=calciumpredtrain;
flyTruePref=odorb(go);
% FIG 3g when trainingonly=0
figure %12
hold on;
xVals = (myprediction-mean(myprediction))/(std(myprediction));
yVals = (flyTruePref-mean(flyTruePref))/(std(flyTruePref));
linreg = linearRegressionCI2(xVals, yVals.', 1, 0, -2.5, 2.5);

areaBar(linreg.xVals,polyval(linreg.pOverall,linreg.xVals),2*std(linreg.fits),[0 0 0],[0.9 0.9 0.9])
plot(xVals,yVals,'.','Color',ocolor, 'LineWidth',3,'MarkerSize',15)

[r p]=corrcoef(xVals,yVals);
text(-.25,-.25,['r = ' num2str(r(1,2),'%2.2f')],'FontSize',15)
text(-.3,-.25,['p = ' num2str(p(1,2),'%2.2f')],'FontSize',15)
text(-.3,-.75,['n = ' num2str(size(xVals),'%2.2f')],'FontSize',15)
xlabel('predicted preference (z-scored)')
ylabel('measured preference (z-scored)')
set(gca,'FontSize',15)
box on
axis([-2.5 2.5 -2.5 2.5])
set(gca,'xtick','')
set(gca,'ytick','')
axis square


%% plot single glomerulus brp-short density alone versus behavior
% close all

glomtouse=3; % 1 for DM2, 3 for DC2
predictormatrix=glomR; 

go=find(odorb);

% use relative fluorescence
figure %13
hold on
for i=go
    plot(mean([predictormatrix(2*(i-1)+1,glomtouse)'; predictormatrix(2*i,glomtouse)']),odorb(i),'o','Color',ocolor,'LineWidth',3)
end
[rr pr]=corrcoef(mean([predictormatrix(2*go-1,glomtouse)'; predictormatrix(2*go,glomtouse)']),odorb(go));
xlabel(['dc2'])
ylabel('odor preference')
title('fluorescence/volume')
set(gca,'FontSize',15)

singleGlomModel = fitlm(mean([predictormatrix(2*(go-1)+1,glomtouse)'; predictormatrix(2*go,glomtouse)']),odorb);

ypredtrain = predict(singleGlomModel,transpose(mean([predictormatrix(2*go-1,glomtouse)'; predictormatrix(2*go,glomtouse)'])));

figure %14
plot(ypredtrain,odorb(go),'o','Color',[0.9 0.2 0.9],'LineWidth',3)
text(0.0,0.1,['r = ' num2str(rr(1,2),'%02.2f')],'FontSize',15)
text(0.0,0.1,['p = ' num2str(pr(1,2),'%02.2f')],'FontSize',15)
xlabel('predicted preference')
ylabel('measured preference')
box off
set(gca,'FontSize',15)


myprediction=ypredtrain;
flyTruePref=odorb(go);
figure %15
plot((myprediction-mean(myprediction))/(std(myprediction)),(flyTruePref-mean(flyTruePref))/(std(flyTruePref)),'.','Color',ocolor, 'LineWidth',3,'MarkerSize',15)
[r p]=corrcoef((myprediction-mean(myprediction))/(std(myprediction)),(flyTruePref-mean(flyTruePref))/(std(flyTruePref)));
text(-.25,-.25,['r = ' num2str(r(1,2),'%2.2f')],'FontSize',15)
text(-.3,-.25,['p = ' num2str(p(1,2),'%2.2f')],'FontSize',15)
text(-.3,-.75,['n = ' num2str(size(xVals),'%2.2f')],'FontSize',15)
xlabel('predicted preference (z-scored)')
ylabel('measured preference (z-scored)')
set(gca,'FontSize',15)
box on
axis([-2.5 2.5 -2.5 2.5])
set(gca,'xtick','')
set(gca,'ytick','')
axis square




%% bootstrap for model selection (for use with training data)


iters=100;
bootstrapsamples=flyNum;

highestPCtouse=4;

measureToUse=glomR;

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
rawPC=cell(1,highestPCtouse);
for jjj=1:iters
    
    % select bootstrap sample flies
    trainflies=zeros(1,bootstrapsamples);
    for i=1:bootstrapsamples
        trainflies(i)=1+round((flyNum-1)*rand(1));
    end
    traindata=zeros(bootstrapsamples,size(glomR,2));
    for i =1:2:(2*bootstrapsamples)
        traindata(i,:)=measureToUse((trainflies((i+1)/2)-1)*2+1,:);
        traindata(i+1,:)=measureToUse(trainflies((i+1)/2)*2,:);
    end

    [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(traindata);
    
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
        behaviorprediction=SCORE(:,pcstouse);
        flyTruePref=zeros(1,length(trainflies));
        flyTruePrefShuffled=zeros(1,length(trainflies));
        ally=odorb;
        
        nactivity=zeros(length(trainflies),length(pcstouse));
        shuffledindtrain=randperm(length(trainflies)); % for comparing to shuffled null model

       for i=1:length(trainflies)
            flyTruePref(i)=ally(trainflies(i));
            flyTruePrefShuffled(i)=ally(trainflies(shuffledindtrain(i)));
            nactivity(i,:)=mean(behaviorprediction([(i-1)*2+1 (i*2)]));
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
            averagePredictor{pcstouse}=currpc';
        else
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
   
    if mod(jjj,10)==0
        disp(['iteration ' num2str(jjj)])
    end
end


for i=1:highestPCtouse
    averagePredictor{i}=averagePredictor{i}/iters;
end

testR2t=testR2';
testR2shuffledt=testR2shuffled';

labs=[];
for j=1:highestPCtouse
    labs=[labs j*ones(1,iters)];
end

% FIG 3d if trainingonly=1
figure %16
boxplot(testR2t(:),labs,'plotstyle','compact','BoxStyle','filled','Colors',ocolor,'medianstyle','target','symbol','','outliersize',1)
xlabel('PC')
ylabel('Unshuffled R^2')
set(gca,'xtick','')
set(gca,'ytick','')
axis([0 highestPCtouse+1 0 0.7])
set(gca,'FontSize',15)
figure %17
boxplot(testR2shuffledt(:),labs,'plotstyle','compact','BoxStyle','filled','Colors',ocolor,'medianstyle','target','symbol','','outliersize',1)
xlabel('PC')
ylabel('Shuffled R^2')
axis([0 highestPCtouse+1 0 0.7])
set(gca,'FontSize',15)
set(gca,'xtick','')
set(gca,'ytick','')
%% apply predictor to test data (assuming model trained with train data only)
pctouse=componentR;

if ~normalize
    glomSt=glomSizeTest;
    glomFt=glomTotalFTest;
    glomRt=glomRelativeFTest;
else
    normGlomRTest =zeros(size(glomRelativeFTest,1),size(glomRelativeFTest,2));
    for i = 1:size(glomRelativeFTest,1)
        normGlomRTest(i,:)= glomRelativeFTest(i,:)/sum(glomRelativeFTest(i,:));
    end
    glomRt=normGlomRTest;    
end
glomSt=glomSt-mean(glomSt);
glomFt=glomFt-mean(glomFt);
glomRt=glomRt-mean(glomRt);
scoreTe=glomRt*coeffR(:,pctouse);


% use relative fluorescence
figure %18
hold on
go2=find(odorbTest);

for i=1:size(go2,2)
    plot(mean([scoreTe(2*go2(i)-1)'; scoreTe(2*go2(i))']),odorbTest(go2(i)),'ko','LineWidth',3)
end
[rr pr]=corrcoef(mean([scoreTe(2*(go2-1)+1)'; scoreTe(2*go2)']),odorbTest(go2));
xlabel(['pc ' num2str(pctouse)])
ylabel('odor preference')
title('fluorescence/volume')
set(gca,'FontSize',15)

ypred = predict(pcmodel,transpose(mean([scoreTe(2*go2-1)'; scoreTe(2*go2)'])));
figure %19
plot(ypred,odorbTest(go2),'o','Color',[0.95 0.2 0.95],'LineWidth',3)
text(0.0,0.1,['r = ' num2str(rr(1,2))],'FontSize',15)
text(0.0,0.1,['p = ' num2str(pr(1,2))],'FontSize',15)
xlabel('predicted preference')
ylabel('measured preference')
box off
set(gca,'FontSize',15)


myprediction=ypred;
flyTruePref=odorbTest(go2);
% FIG 3f when trainingonly=0
figure %20
hold on;
xVals = (myprediction-mean(myprediction))/(std(myprediction));
yVals = (flyTruePref-mean(flyTruePref))/(std(flyTruePref));
linreg = linearRegressionCI2(xVals, yVals.', 1, 0, -2.2, 2.2);

areaBar(linreg.xVals,polyval(linreg.pOverall,linreg.xVals),2*std(linreg.fits),[0 0 0],[0.9 0.9 0.9])
plot(xVals,yVals,'.','Color',ocolor, 'LineWidth',3,'MarkerSize',15)

[r p]=corrcoef(xVals,yVals);
text(-.25,-.25,['r = ' num2str(r(1,2),'%2.2f')],'FontSize',15)
text(-.3,-.25,['p = ' num2str(p(1,2),'%2.2f')],'FontSize',15)
text(-.3,-.75,['n = ' num2str(size(xVals),'%2.2f')],'FontSize',15)
xlabel('predicted preference (z-scored)')
ylabel('measured preference (z-scored)')
set(gca,'FontSize',15)
box on
axis([-2.2 2.2 -2.2 2.2])
set(gca,'xtick','')
set(gca,'ytick','')
axis square


figure %21
hold on
for i=go2
    plot(mean([(glomRt(2*(i-1)+1,1)'-glomRt(2*(i-1)+1,3)'); (glomRt(2*(i),1)'-glomRt(2*(i),3)')]),odorbTest(i),'ko','LineWidth',3)
end
[rd pd]=corrcoef(mean([glomRt(2*(go2-1)+1,1)'-glomRt(2*(go2-1)+1,3)'; glomRt(2*(go2),1)'-glomRt(2*(go2),3)']),odorbTest(go2));
text(0.0,0.1,['r = ' num2str(rd(1,2))],'FontSize',15)
text(0.0,0.1,['p = ' num2str(pd(1,2))],'FontSize',15)
xlabel('DM2-DC2 fluorescence density')
ylabel('odor preference')
box off
set(gca,'FontSize',15)



calciumpred = predict(calciummodel,transpose(mean([(glomRt(2*(go2-1)+1,1)'-glomRt(2*(go2-1)+1,3)'); (glomRt(2*(go2),1)'-glomRt(2*(go2),3)')])));
figure %22
hold on
plot(calciumpred,odorbTest(go2),'o','Color',[0.95 0.2 0.95],'LineWidth',3)
text(0.0,0.1,['r = ' num2str(-rd(1,2),'%02.2f')],'FontSize',15)
text(0.0,0.1,['p = ' num2str(pd(1,2),'%02.2f')],'FontSize',15)
xlabel('predicted preference')
ylabel('measured preference')
box off
set(gca,'FontSize',15)

%% bootstrap resample test data
nsamples = length(go2);
niters = 10000;

for i = 1:niters
    inds = round(rand(nsamples,1)*(nsamples-1))+1;
    [rr pr]=corrcoef(mean([scoreTe(2*(inds-1)+1)'; scoreTe(2*inds)']),odorbTest(inds));
    [rd pd]=corrcoef(mean([glomRt(2*(inds-1)+1,1)'-glomRt(2*(inds-1)+1,3)'; glomRt(2*(inds),1)'-glomRt(2*(inds),3)']),odorbTest(inds));
    
    rpc(i) = rr(1,2);
    ppcp(i) = pr(1,2);
    rdm2dc2(i) = rd(1,2);
    
end
ppc = sum(rpc<0)/length(rpc)
pdm2dc2 = sum(rdm2dc2>0)/length(rdm2dc2)

figure %23
histogram(rpc,'Normalization','probability')
xlabel('r')
ylabel('probability')
box off
hold on
histogram(-rdm2dc2,'Normalization','probability')
xlabel('r')
ylabel('probability')
box off
legend('top PC-based prediction','calcium prior-based prediction')
legend boxoff
text(0,0.01,['p1 = ' num2str(ppc)],'FontSize',15)
text(0,0.01,['p2 = ' num2str(pdm2dc2)],'FontSize',15)
set(gca,'FontSize',15)


figure %24
histogram(rpc,'Normalization','probability')
title('bootstrap resampling')
xlabel('r')
ylabel('probability')
box off
text(0,0.01,['p = ' num2str(ppc)],'FontSize',15)
set(gca,'FontSize',15)

% FIG 3c
% figure 25
example_fly_mask = fullfile(analysis_dir_path, '\IHC\OrcoBrpshort\autoSegmentation\190404_orcobrpshort_behaviorAndImaging\fly12__autoseg_xystep11_labelled.mat');
showIHCsegmented(example_fly_mask)
