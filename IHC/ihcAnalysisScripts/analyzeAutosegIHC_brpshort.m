%% load data from preprocessed file
clear all
close all
load orcobrpshort_all_data
%% load data from raw file 
clear all
close all
apple=1;
if apple==1
    startdir='/Users/mattchurgin/Dropbox/flyimaging/analysis/IHC/OrcoBrpshort/autoSegmentation';
    manualSegDir='/Users/mattchurgin/Dropbox/flyimaging/analysis/IHC/OrcoBrpshort/orcoOX5brpshortOXZ';
else
    startdir='C:\Users\mac0456\Dropbox\flyimaging\analysis\IHC\OrcoBrpshort\autoSegmentation';
    manualSegDir='C:\Users\mac0456\Dropbox\flyimaging\analysis\IHC\OrcoBrpshort\orcoOX5brpshortOXZ';
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
trainingonly =1; % 0 to include train + test data, 1 to include only train data for model building
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






%% plot manual segmentation results vs. automated

figure;
plot(dm2lv,glomSize(1:2:end,1),'o','Color',[0 0 0],'LineWidth',3,'MarkerSize',10)
hold on
plot(dm2rv,glomSize(2:2:end,1),'x','Color',[0 0 0],'LineWidth',3,'MarkerSize',10)
plot(dc2lv,glomSize(1:2:end,3),'o','Color',[0.7 0.7 0.7],'LineWidth',3,'MarkerSize',10)
plot(dc2rv,glomSize(2:2:end,3),'x','Color',[0.7 0.7 0.7],'LineWidth',3,'MarkerSize',10)
xlabel('manual')
ylabel('automated')
legend('DM2 left','DM2 right','DC2 left','DC2 right')
legend boxoff
box off
title('volume (voxels)')
set(gca,'FontSize',15)

figure;
plot(dm2lp,glomTotalF(1:2:end,1),'o','Color',[0 0 0],'LineWidth',3,'MarkerSize',10)
hold on
plot(dm2rp,glomTotalF(2:2:end,1),'x','Color',[0 0 0],'LineWidth',3,'MarkerSize',10)
plot(dc2lp,glomTotalF(1:2:end,3),'o','Color',[0.7 0.7 0.7],'LineWidth',3,'MarkerSize',10)
plot(dc2rp,glomTotalF(2:2:end,3),'x','Color',[0.7 0.7 0.7],'LineWidth',3,'MarkerSize',10)
xlabel('manual')
ylabel('automated')
legend('DM2 left','DM2 right','DC2 left','DC2 right')
legend boxoff
box off
title('total fluorescence (a.u.)')
set(gca,'FontSize',15)

%%  plot glomerulus properties

figure
subplot(1,3,1)
plot(glomSize((1:flyNum)*2-1,:),glomSize((1:flyNum)*2,:),'.','LineWidth',2,'MarkerSize',15)
title('volume (um^3)')
xlabel('left lobe')
ylabel('right lobe')
set(gca,'FontSize',15)
axis square
axis([0 6200 0 6200])
subplot(1,3,2)
plot(glomTotalF((1:flyNum)*2-1,:),glomTotalF((1:flyNum)*2,:),'.','LineWidth',2,'MarkerSize',15)
title('total fluorescence')
xlabel('left lobe')
ylabel('right lobe')
set(gca,'FontSize',15)
axis square
axis([0 5e5 0 5e5])
subplot(1,3,3)
plot(glomRelativeF((1:flyNum)*2-1,:),glomRelativeF((1:flyNum)*2,:),'.','LineWidth',2,'MarkerSize',15)
legend(glomeruliToSegment)
legend boxoff
title('normalized fluorescence')
xlabel('left lobe')
ylabel('right lobe')
set(gca,'FontSize',15)
axis square
axis([10 80 10 80])

violinPlot(glomSize,[ocolor])
set(gca,'XTickLabel',glomeruliToSegment)
xtickangle(30)
ylabel('volume (um^3)')
axis([0 5 0 1.2*max(max(glomSize))])
box on
set(gca,'FontSize',15)

violinPlot(glomTotalF,[ocolor])
set(gca,'XTickLabel',glomeruliToSegment)
xtickangle(30)
ylabel('total fluorescence')
axis([0 5 0 1.2*max(max(glomTotalF))])
box on
set(gca,'FontSize',15)

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

figure
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

figure
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

figure
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
close all

componentR=2;

go=find(odorb);

% use relative fluorescence
figure
hold on
for i=go
    plot(mean([scoreR(2*(i-1)+1,componentR)'; scoreR(2*i,componentR)']),odorb(i),'o','Color',[0.9 0.2 0.9],'LineWidth',3)
    %text(mean([scoreR(2*(i-1)+1,componentR)'; scoreR(2*i,componentR)'])+0.001,odorb(i),num2str(i),'FontSize',15)
end
[rr pr]=corrcoef(mean([scoreR(2*go-1,componentR)'; scoreR(2*go,componentR)']),odorb(go));
xlabel(['pc ' num2str(componentR)])
ylabel('odor preference')
title('fluorescence/volume')
set(gca,'FontSize',15)

pcmodel = fitlm(mean([scoreR(2*(go-1)+1,componentR)'; scoreR(2*go,componentR)']),odorb);

ypredtrain = predict(pcmodel,transpose(mean([scoreR(2*go-1,componentR)'; scoreR(2*go,componentR)'])));

myprediction=ypredtrain;
flyTruePref=odorb(go);
figure
plot((myprediction-mean(myprediction))/(std(myprediction)),(flyTruePref-mean(flyTruePref))/(std(flyTruePref)),'.','Color',ocolor, 'LineWidth',3,'MarkerSize',15)
[r p]=corrcoef((myprediction-mean(myprediction))/(std(myprediction)),(flyTruePref-mean(flyTruePref))/(std(flyTruePref)));
text(-.25,-.25,['r = ' num2str(r(1,2),'%2.2f')],'FontSize',15)
text(-.3,-.25,['p = ' num2str(p(1,2),'%2.2f')],'FontSize',15)
xlabel('predicted preference (z-scored)')
ylabel('measured preference (z-scored)')
set(gca,'FontSize',15)
box on
axis([-2.2 2.2 -2.2 2.2])
set(gca,'xtick','')
set(gca,'ytick','')
axis square



%% plot DM2 - DC2 versus behavior

figure
hold on
for i=go
    % plot(mean([glomS(2*(i-1)+1,1)'-glomS(2*(i-1)+1,3)'; glomS(2*(i),1)'-glomS(2*(i),3)']),odorb(i),'mo','LineWidth',3)
    % plot(mean([glomF(2*(i-1)+1,1)'-glomF(2*(i-1)+1,3)'; glomF(2*(i),1)'-glomF(2*(i),3)']),odorb(i),'ro','LineWidth',3)
    plot(mean([(glomR(2*(i-1)+1,1)'-glomR(2*(i-1)+1,3)'); (glomR(2*(i),1)'-glomR(2*(i),3)')]),odorb(i),'mo','LineWidth',3)
    %text(mean([(glomR(2*(i-1)+1,1)'-glomR(2*(i-1)+1,3)'); (glomR(2*(i),1)'-glomR(2*(i),3)')]),odorb(i),num2str(i),'FontSize',15)
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
figure
plot((myprediction-mean(myprediction))/(std(myprediction)),(flyTruePref-mean(flyTruePref))/(std(flyTruePref)),'.','Color',ocolor, 'LineWidth',3,'MarkerSize',15)
[r p]=corrcoef((myprediction-mean(myprediction))/(std(myprediction)),(flyTruePref-mean(flyTruePref))/(std(flyTruePref)));
text(-.25,-.25,['r = ' num2str(r(1,2),'%2.2f')],'FontSize',15)
text(-.3,-.25,['p = ' num2str(p(1,2),'%2.2f')],'FontSize',15)
xlabel('predicted preference (z-scored)')
ylabel('measured preference (z-scored)')
set(gca,'FontSize',15)
box on
axis([-2.5 2.5 -2.5 2.5])
set(gca,'xtick','')
set(gca,'ytick','')
axis square


%% plot single glomerulus brp-short density alone versus behavior
close all

glomtouse=3; % 1 for DM2, 3 for DC2
predictormatrix=glomR; 

go=find(odorb);

% use relative fluorescence
figure
hold on
for i=go
    plot(mean([predictormatrix(2*(i-1)+1,glomtouse)'; predictormatrix(2*i,glomtouse)']),odorb(i),'o','Color',ocolor,'LineWidth',3)
    %text(mean([scoreR(2*(i-1)+1,componentR)'; scoreR(2*i,componentR)'])+0.001,odorb(i),num2str(i),'FontSize',15)
end
[rr pr]=corrcoef(mean([predictormatrix(2*go-1,glomtouse)'; predictormatrix(2*go,glomtouse)']),odorb(go));
xlabel(['dc2'])
ylabel('odor preference')
title('fluorescence/volume')
set(gca,'FontSize',15)

singleGlomModel = fitlm(mean([predictormatrix(2*(go-1)+1,glomtouse)'; predictormatrix(2*go,glomtouse)']),odorb);

ypredtrain = predict(singleGlomModel,transpose(mean([predictormatrix(2*go-1,glomtouse)'; predictormatrix(2*go,glomtouse)'])));
figure
plot(ypredtrain,odorb(go),'o','Color',[0.9 0.2 0.9],'LineWidth',3)
text(0.0,0.1,['r = ' num2str(rr(1,2),'%02.2f')],'FontSize',15)
text(0.0,0.1,['p = ' num2str(pr(1,2),'%02.2f')],'FontSize',15)
xlabel('predicted preference')
ylabel('measured preference')
box off
set(gca,'FontSize',15)


myprediction=ypredtrain;
flyTruePref=odorb(go);
figure
plot((myprediction-mean(myprediction))/(std(myprediction)),(flyTruePref-mean(flyTruePref))/(std(flyTruePref)),'.','Color',ocolor, 'LineWidth',3,'MarkerSize',15)
[r p]=corrcoef((myprediction-mean(myprediction))/(std(myprediction)),(flyTruePref-mean(flyTruePref))/(std(flyTruePref)));
text(-.25,-.25,['r = ' num2str(r(1,2),'%2.2f')],'FontSize',15)
text(-.3,-.25,['p = ' num2str(p(1,2),'%2.2f')],'FontSize',15)
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
        %traindata(i,:)=glomR([(trainflies(i)-1)*2+1 trainflies(i)*2],:);
        traindata(i,:)=glomR((trainflies((i+1)/2)-1)*2+1,:);
        traindata(i+1,:)=glomR(trainflies((i+1)/2)*2,:);
    end
%     
%     flyindicestrain=cell(1,length(trainflies));
%     datapts=cell(1,length(trainflies));
%     lastdatapt=0;
%     for i=1:length(trainflies)
%         [temp temp2]=find(glombyodorflies==trainflies(i));
%         [tem3 temp4]=find(glombyodorflies==(trainflies(i)+leftRightSplitVal));
%         flyindicestrain{i}=[temp2 temp4];
%         datapts{i}=[(lastdatapt+1):(lastdatapt+length(flyindicestrain{i}))];
%         lastdatapt=(lastdatapt+length(flyindicestrain{i}));
%     end
%     
%     traindata=[];
%     for i=1:length(trainflies)
%         for j=1:length(flyindicestrain{i})
%             traindata=[traindata responsesNoResponseRemoved(:,flyindicestrain{i}(j))];
%         end
%     end
%     
    
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
        %behaviorprediction=(traindata*currpc);
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
%labs=[ones(1,iters) 2*ones(1,iters) 3*ones(1,iters)];
labs=[labs j*ones(1,iters)];
end

figure
boxplot(testR2t(:),labs,'plotstyle','compact','BoxStyle','filled','Colors',ocolor,'medianstyle','target','symbol','','outliersize',1)
xlabel('PC')
ylabel('Unshuffled R^2')
set(gca,'xtick','')
set(gca,'ytick','')
axis([0 highestPCtouse+1 0 0.7])
set(gca,'FontSize',15)
figure
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
glomRt=glomRt-mean(glomRt);
scoreTe=glomRt*coeffR(:,pctouse);


% use relative fluorescence
figure
hold on
go2=find(odorbTest);

for i=1:size(go2,2)
    plot(mean([scoreTe(2*go2(i)-1)'; scoreTe(2*go2(i))']),odorbTest(go2(i)),'ko','LineWidth',3)
    %plot(scoreTe(2*(go2(i))-1)',odorbTest(go2(i)),'ko','LineWidth',3)
    %text(mean([scoreTe(2*go2(i)-1)'; scoreTe(2*go2(i))'])+0.01,odorbTest(go2(i)),num2str(i),'FontSize',15)
end
[rr pr]=corrcoef(mean([scoreTe(2*(go2-1)+1)'; scoreTe(2*go2)']),odorbTest(go2));
xlabel(['pc ' num2str(pctouse)])
ylabel('odor preference')
title('fluorescence/volume')
set(gca,'FontSize',15)

ypred = predict(pcmodel,transpose(mean([scoreTe(2*go2-1)'; scoreTe(2*go2)'])));
figure
plot(ypred,odorbTest(go2),'o','Color',[0.95 0.2 0.95],'LineWidth',3)
text(0.0,0.1,['r = ' num2str(rr(1,2))],'FontSize',15)
text(0.0,0.1,['p = ' num2str(pr(1,2))],'FontSize',15)
xlabel('predicted preference')
ylabel('measured preference')
box off
set(gca,'FontSize',15)


myprediction=ypred;
flyTruePref=odorbTest(go2);
figure
plot((myprediction-mean(myprediction))/(std(myprediction)),(flyTruePref-mean(flyTruePref))/(std(flyTruePref)),'.','Color',ocolor, 'LineWidth',3,'MarkerSize',15)
[r p]=corrcoef((myprediction-mean(myprediction))/(std(myprediction)),(flyTruePref-mean(flyTruePref))/(std(flyTruePref)));
text(-.25,-.25,['r = ' num2str(r(1,2),'%2.2f')],'FontSize',15)
text(-.3,-.25,['p = ' num2str(p(1,2),'%2.2f')],'FontSize',15)
xlabel('predicted preference (z-scored)')
ylabel('measured preference (z-scored)')
set(gca,'FontSize',15)
box on
axis([-2.2 2.2 -2.2 2.2])
set(gca,'xtick','')
set(gca,'ytick','')
axis square

figure
hold on
for i=go2
    plot(mean([(glomRt(2*(i-1)+1,1)'-glomRt(2*(i-1)+1,3)'); (glomRt(2*(i),1)'-glomRt(2*(i),3)')]),odorbTest(i),'ko','LineWidth',3)
    %text(mean([(glomRt(2*(i-1)+1,1)'-glomRt(2*(i-1)+1,3)'); (glomRt(2*(i),1)'-glomRt(2*(i),3)')]),odorbTest(i),num2str(i),'FontSize',15)
end
[rd pd]=corrcoef(mean([glomRt(2*(go2-1)+1,1)'-glomRt(2*(go2-1)+1,3)'; glomRt(2*(go2),1)'-glomRt(2*(go2),3)']),odorbTest(go2));
text(0.0,0.1,['r = ' num2str(rd(1,2))],'FontSize',15)
text(0.0,0.1,['p = ' num2str(pd(1,2))],'FontSize',15)
xlabel('DM2-DC2 fluorescence density')
ylabel('odor preference')
box off
set(gca,'FontSize',15)



calciumpred = predict(calciummodel,transpose(mean([(glomRt(2*(go2-1)+1,1)'-glomRt(2*(go2-1)+1,3)'); (glomRt(2*(go2),1)'-glomRt(2*(go2),3)')])));
figure
%plot(calciumpredtrain,odorb(go),'k*','LineWidth',3)
hold on
plot(calciumpred,odorbTest(go2),'o','Color',[0.95 0.2 0.95],'LineWidth',3)
text(0.0,0.1,['r = ' num2str(-rd(1,2),'%02.2f')],'FontSize',15)
text(0.0,0.1,['p = ' num2str(pd(1,2),'%02.2f')],'FontSize',15)
xlabel('predicted preference')
ylabel('measured preference')
%legend('train','test')
%legend boxoff
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

figure
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


figure
histogram(rpc,'Normalization','probability')
title('bootstrap resampling')
xlabel('r')
ylabel('probability')
box off
text(0,0.01,['p = ' num2str(ppc)],'FontSize',15)
set(gca,'FontSize',15)
