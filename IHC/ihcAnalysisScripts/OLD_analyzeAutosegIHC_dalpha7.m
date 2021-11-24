clear all
close all
apple=0;
if apple==1
    startdir='/Users/mattchurgin/Dropbox/flyimaging/analysis/IHC/Gh146Dalpha7/autoSegmentation';
    manualSegDir='/Users/mattchurgin/Dropbox/flyimaging/analysis/IHC/Gh146Dalpha7/gh146OX10Dalpha7OXZ';
else
    startdir='C:\Users\mac0456\Dropbox\flyimaging\analysis\IHC\Gh146Dalpha7\autoSegmentation';
    manualSegDir='C:\Users\mac0456\Dropbox\flyimaging\analysis\IHC\Gh146Dalpha7\gh146OX10Dalpha7OXZ';
end
cd(startdir)

folders{1}='190304_gh146_dalpha7OXZ';
folders{2}='190301_gh146_dalpha7OXZ';
folders{3}='190214_gh146_Dalpha7gfp_myr_behaviorANDIHC';
folders{4}='190507_gh146dalpha7_behaviorPlusIHC';
folders{5}='190610_gh146Dalpha7_behaviorandihc';

densityOrTotal=1;
sizeConversion=((.09*13)^2)*0.49;

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
            flyNum=flyNum+1;
            for leftright=[1 2]
                for k=1:length(glomeruliToSegment)
                    glomSize(2*(flyNum-1)+leftright,k)=sizeConversion*sum(glomeruliSegs{leftright,k}(:));
                    glomTotalF(2*(flyNum-1)+leftright,k)=sum(b(:).*glomeruliSegs{leftright,k}(:));
                    glomRelativeF(2*(flyNum-1)+leftright,k)=sum(b(:).*glomeruliSegs{leftright,k}(:))/sum(glomeruliSegs{leftright,k}(:));
                end
            end
            
            % load behavior
            load([manualSegDir '/' currdate '_' currFlyName '/' currFlyName '_behavior.mat'])
            
            odorb(flyNum)=occ-preocc;
            odorpre(flyNum)=preocc;
            
            % load manual segmentation data
            manualfiles=dir([manualSegDir '/' currdate '_' currFlyName]);
            for kkk=1:length(manualfiles)
                if strfind(manualfiles(kkk).name,'manuallysegmented')
                    load([manualSegDir '/' currdate '_' currFlyName '/' manualfiles(kkk).name])
                    break
                end
            end
            
            if densityOrTotal==0
                punctadiffL(flyNum)=100*(dm2l_volnormed-dc2l_volnormed)/(dm2l_volnormed+dc2l_volnormed);
                punctadiffR(flyNum)=100*(dm2r_volnormed-dc2r_volnormed)/(dm2r_volnormed+dc2r_volnormed);
                punctadiffM(flyNum)=nanmean([punctadiffL(i) punctadiffR(i)]);
                % save puncta density for each glomerulus
                
                dm2lp(flyNum)=dm2l_volnormed;
                dm2rp(flyNum)=dm2r_volnormed;
                dc2lp(flyNum)=dc2l_volnormed;
                dc2rp(flyNum)=dc2r_volnormed;
                dm2lv(flyNum)=sum(sum(sum(dm2Lmask)));
                dm2rv(flyNum)=sum(sum(sum(dm2Rmask)));
                dc2lv(flyNum)=sum(sum(sum(dc2Lmask)));
                dc2rv(flyNum)=sum(sum(sum(dc2Rmask)));
            else
                % calculate total fluorescence rather than density
                dm2lp(flyNum)=dm2l_volnormed*sum(sum(sum(dm2Lmask)));
                dm2rp(flyNum)=dm2r_volnormed*sum(sum(sum(dm2Rmask)));
                dc2lp(flyNum)=dc2l_volnormed*sum(sum(sum(dc2Lmask)));
                dc2rp(flyNum)=dc2r_volnormed*sum(sum(sum(dc2Rmask)));
                punctadiffL(flyNum)=100*(dm2lp(i)-dc2lp(i))/(dm2lp(i)+dc2lp(i));
                punctadiffR(flyNum)=100*(dm2rp(i)- dc2rp(i))/(dm2rp(i)+ dc2rp(i));
                punctadiffM(flyNum)=nanmean([punctadiffL(i) punctadiffR(i)]);
                
                dm2lv(flyNum)=sum(sum(sum(dm2Lmask)));
                dm2rv(flyNum)=sum(sum(sum(dm2Rmask)));
                dc2lv(flyNum)=sum(sum(sum(dc2Lmask)));
                dc2rv(flyNum)=sum(sum(sum(dc2Rmask)));
            end
            
            disp(['loaded fly ' num2str(flyNum)])
        end
    end
    
    cd ..
end


disp('done loading')

%% correct for odor preference bimodality
bimodalSplit=-0.1;
odorOrig=odorb;
odorb(odorOrig<bimodalSplit)=odorb(odorOrig<bimodalSplit)-mean(odorb(odorOrig<bimodalSplit));
odorb(odorOrig>=bimodalSplit)=odorb(odorOrig>=bimodalSplit)-mean(odorb(odorOrig>=bimodalSplit));
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

%% plot glomerulus properties

nbins=5;
figure

subplot(1,3,3)
h1 = histogram(glomRelativeF(:,1),'Normalization','pdf');
h1.NumBins=nbins;
hold on
h2 = histogram(glomRelativeF(:,2),'Normalization','pdf');
h2.NumBins=nbins;
h3 = histogram(glomRelativeF(:,3),'Normalization','pdf');
h3.NumBins=nbins;
h4 = histogram(glomRelativeF(:,4),'Normalization','pdf');
h4.NumBins=nbins;
h5 = histogram(glomRelativeF(:,5),'Normalization','pdf');
h5.NumBins=nbins;
legend('DM2','DM3','DC2','DL5','DM1')
legend boxoff
xlabel('fluorescence/volume')
ylabel('probability')
set(gca,'FontSize',15)
box off


subplot(1,3,1)

h1 = histogram(glomSize(:,1),'Normalization','pdf');
h1.NumBins=nbins;
hold on
h2 = histogram(glomSize(:,2),'Normalization','pdf');
h2.NumBins=nbins;
h3 = histogram(glomSize(:,3),'Normalization','pdf');
h3.NumBins=nbins;
h4 = histogram(glomSize(:,4),'Normalization','pdf');
h4.NumBins=nbins;
h5 = histogram(glomSize(:,5),'Normalization','pdf');
h5.NumBins=nbins;
legend('DM2','DM3','DC2','DL5','DM1')
legend boxoff
xlabel('volume (voxels)')
ylabel('probability')
set(gca,'FontSize',15)
box off



subplot(1,3,2)

h1 = histogram(glomTotalF(:,1),'Normalization','pdf');
h1.NumBins=nbins;
hold on
h2 = histogram(glomTotalF(:,2),'Normalization','pdf');
h2.NumBins=nbins;
h3 = histogram(glomTotalF(:,3),'Normalization','pdf');
h3.NumBins=nbins;
h4 = histogram(glomTotalF(:,4),'Normalization','pdf');
h4.NumBins=nbins;
h5 = histogram(glomTotalF(:,5),'Normalization','pdf');
h5.NumBins=nbins;
legend('DM2','DM3','DC2','DL5','DM1')
legend boxoff
xlabel('total fluorescence')
ylabel('probability')
set(gca,'FontSize',15)
box off



figure
subplot(1,3,1)
distributionPlot(glomSize,'histOpt',1,'colormap',1-gray(64),'showMM',0)
set(gca,'XTickLabel',glomeruliToSegment)
xtickangle(30)
ylabel('volume (um^3)')
axis([0 6 0 1.2*max(max(glomSize))])
box off
set(gca,'FontSize',15)

subplot(1,3,2)
distributionPlot(glomTotalF,'histOpt',1,'colormap',1-gray(64),'showMM',0)
set(gca,'XTickLabel',glomeruliToSegment)
xtickangle(30)
ylabel('total fluorescence')
axis([0 6 0 1.2*max(max(glomTotalF))])
box off
set(gca,'FontSize',15)

subplot(1,3,3)
distributionPlot(glomRelativeF,'histOpt',1,'colormap',1-gray(64),'showMM',0)
set(gca,'XTickLabel',glomeruliToSegment)
xtickangle(30)
ylabel('fluorescence/volume')
axis([0 6 0 1.2*max(max(glomRelativeF))])
box off
set(gca,'FontSize',15)



figure
subplot(1,2,1)
distributionPlot(glomSize,'histOpt',1,'colormap',1-gray(64),'showMM',0)
set(gca,'XTickLabel',glomeruliToSegment)
xtickangle(30)
ylabel('volume (um^3)')
axis([0 6 0 1.2*max(max(glomSize))])
box off
set(gca,'FontSize',15)

subplot(1,2,2)
distributionPlot(glomRelativeF,'histOpt',1,'colormap',1-gray(64),'showMM',0)
set(gca,'XTickLabel',glomeruliToSegment)
xtickangle(30)
ylabel('fluorescence/volume')
axis([0 6 0 1.2*max(max(glomRelativeF))])
box off
set(gca,'FontSize',15)
%% perform pca on auto segmented data


% standardize
glomS=(glomSize-mean(glomSize))./(std(glomSize));
glomF=(glomTotalF-mean(glomTotalF))./(std(glomTotalF));
glomR=(glomRelativeF-mean(glomRelativeF))./(std(glomRelativeF));

% glomS=glomSize;
% slomF=glomTotalF;
%glomR=glomRelativeF;
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
%% plot principal component vs behavior. use pc with largest dm2-dc2 difference
componentS=5;
componentF=5;
componentR=1;
go=find(odorOrig);%>-0.45);

%go=find(odorb<0.5);

% use glomerulus volume
figure
hold on
for i=go
    plot(mean([scoreS(2*(i-1)+1,componentS)'; scoreS(2*i,componentS)']),odorb(i),'ko','LineWidth',3)
    text(mean([scoreS(2*(i-1)+1,componentS)'; scoreS(2*i,componentS)'])+0.01,odorb(i),num2str(i),'FontSize',15)
end
[rs ps]=corrcoef(mean([scoreS(2*go-1,componentS)'; scoreS(2*go,componentS)']),odorb(go));
xlabel(['pc ' num2str(componentS)])
ylabel('odor preference')
title('volume')
set(gca,'FontSize',15)

% use total fluorescence
figure
hold on
for i=go
    plot(mean([scoreF(2*(i-1)+1,componentF)'; scoreF(2*i,componentF)']),odorb(i),'ko','LineWidth',3)
    text(mean([scoreF(2*(i-1)+1,componentF)'; scoreF(2*i,componentF)'])+0.01,odorb(i),num2str(i),'FontSize',15)
end
[rf pf]=corrcoef(mean([scoreF(2*go-1,componentF)'; scoreF(2*go,componentF)']),odorb(go));
xlabel(['pc ' num2str(componentF)])
ylabel('odor preference')
title('total fluoresence')
set(gca,'FontSize',15)

% use relative fluorescence
figure
hold on
for i=go
    plot(mean([scoreR(2*(i-1)+1,componentR)'; scoreR(2*i,componentR)']),odorb(i),'ko','LineWidth',3)
    %text(mean([scoreR(2*(i-1)+1,componentR)'; scoreR(2*i,componentR)'])+0.01,odorb(i),num2str(i),'FontSize',15)
end
[rr pr]=corrcoef(mean([scoreR(2*go-1,componentR)'; scoreR(2*go,componentR)']),odorb(go));
xlabel(['pc ' num2str(componentR)])
ylabel('odor preference')
title('fluorescence/volume')
set(gca,'FontSize',15)


% use dm2-dc2 relative fluorescence
figure
hold on
for i=go
   % plot(mean([glomS(2*(i-1)+1,1)'-glomS(2*(i-1)+1,3)'; glomS(2*(i),1)'-glomS(2*(i),3)']),odorb(i),'mo','LineWidth',3)
   % plot(mean([glomF(2*(i-1)+1,1)'-glomF(2*(i-1)+1,3)'; glomF(2*(i),1)'-glomF(2*(i),3)']),odorb(i),'ro','LineWidth',3)
    plot(mean([(glomR(2*(i-1)+1,1)'-glomR(2*(i-1)+1,3)')./(glomR(2*(i-1)+1,1)'+glomR(2*(i-1)+1,3)'); (glomR(2*(i),1)'-glomR(2*(i),3)')./(glomR(2*(i),1)'-glomR(2*(i),3)')]),odorb(i),'ko','LineWidth',3)
    text(mean([(glomR(2*(i-1)+1,1)'-glomR(2*(i-1)+1,3)')./(glomR(2*(i-1)+1,1)'+glomR(2*(i-1)+1,3)'); (glomR(2*(i),1)'-glomR(2*(i),3)')./(glomR(2*(i),1)'-glomR(2*(i),3)')]),odorb(i),num2str(i),'FontSize',15)
end
[rdn pdn]=corrcoef(mean([(glomR(2*(go-1)+1,1)'-glomR(2*(go-1)+1,3)')./(glomR(2*(go-1)+1,1)'+glomR(2*(go-1)+1,3)'); (glomR(2*(go),1)'-glomR(2*(go),3)')./(glomR(2*(go),1)'-glomR(2*(go),3)')]),odorb(go));
xlabel('(DM2-DC2)/(DM2+DC2)')
ylabel('odor preference')
box off
set(gca,'FontSize',15)

figure
hold on
for i=go
   % plot(mean([glomS(2*(i-1)+1,1)'-glomS(2*(i-1)+1,3)'; glomS(2*(i),1)'-glomS(2*(i),3)']),odorb(i),'mo','LineWidth',3)
   % plot(mean([glomF(2*(i-1)+1,1)'-glomF(2*(i-1)+1,3)'; glomF(2*(i),1)'-glomF(2*(i),3)']),odorb(i),'ro','LineWidth',3)
    plot(mean([(glomR(2*(i-1)+1,1)'-glomR(2*(i-1)+1,3)'); (glomR(2*(i),1)'-glomR(2*(i),3)')]),odorb(i),'ko','LineWidth',3)
    text(mean([(glomR(2*(i-1)+1,1)'-glomR(2*(i-1)+1,3)'); (glomR(2*(i),1)'-glomR(2*(i),3)')]),odorb(i),num2str(i),'FontSize',15)
end
[rd pd]=corrcoef(mean([glomR(2*(go-1)+1,1)'-glomR(2*(go-1)+1,3)'; glomR(2*(go),1)'-glomR(2*(go),3)']),odorb(go));
xlabel('(DM2-DC2)')
ylabel('odor preference')
box off
set(gca,'FontSize',15)

%% run train/test many times and save replicate data
warning('off')


go=find(odorOrig>-0.45);

newFlyNum=length(go);
newOdorb=odorb(go);

iters=100;
highestPCtouse=5;

datam=glomR(sort([(go-1)*2+1 go*2]),:);

testR=zeros(highestPCtouse,iters);
testRshuffled=zeros(highestPCtouse,iters);
testNRSS=zeros(highestPCtouse,iters); % normalized residual sum of squares
averagePredictor=cell(highestPCtouse,1);
bestPredictorLoading=cell(iters,1); % loadings for the best predictor
bestPredictor=zeros(iters,1); % PC which is the best predictor
bestPredictorR=zeros(iters,1);
pcPredictorRank=zeros(iters,highestPCtouse);
vExplained=cell(1,highestPCtouse);
rawPC=cell(1,highestPCtouse);

for jjj=1:iters
    
    % Split data into randomly assigned train/test sets
    testsize=40;
    trainsize=100-testsize;
    randomizedOrder=randperm(newFlyNum);
    holdoutflies=randomizedOrder(1:round((length(randomizedOrder)*testsize/100)));
    trainflies=setxor(1:newFlyNum,holdoutflies);
    traintoremove=[];
      
       
    traindata=datam(sort([(trainflies-1)*2+1 trainflies*2]),:);
   testdata=datam(sort([(holdoutflies-1)*2+1 holdoutflies*2]),:);

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
        behaviorprediction=(datam*currpc);
        flyTruePref=zeros(1,length(trainflies));
        flyTruePrefShuffled=zeros(1,length(trainflies));
        ally=newOdorb';
               
        nactivity=zeros(length(trainflies),length(pcstouse));
        shuffledindtrain=randperm(length(trainflies)); % for comparing to shuffled null model
        for i=1:length(trainflies)
            flyTruePref(i)=ally(trainflies(i));
            flyTruePrefShuffled(i)=ally(trainflies(shuffledindtrain(i)));
            nactivity(i,:)=mean(behaviorprediction(((trainflies(i)-1)*2+1):(trainflies(i)*2)));
        end
        linmodel=fitlm(nactivity,flyTruePref);
        
        beta=linmodel.Coefficients.Estimate;
        
        PCContribution=currpc*beta(2:end);
        
        % apply predictor to testdata and evaluate how held out points fit in
        % model
        
        %behaviorprediction=testdata'*PCContribution;
        
        testflyTruePref=zeros(1,length(holdoutflies));
        testflyTruePrefShuffled=zeros(1,length(holdoutflies));
        flyPredictedPref=zeros(1,length(holdoutflies));
        
        testnactivity=zeros(length(holdoutflies),1);
        shuffledindtest=randperm(length(holdoutflies)); % for comparing to shuffled null model
        for i=1:length(holdoutflies)
            testflyTruePref(i)=ally(holdoutflies(i));
            testflyTruePrefShuffled(i)=ally(holdoutflies(shuffledindtest(i)));
            testnactivity(i,:)=mean(behaviorprediction(((holdoutflies(i)-1)*2+1):(holdoutflies(i)*2)));
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

figure
subplot(1,2,1)
distributionPlot(myr2,'histOpt',1,'colormap',1-gray(64),'showMM',0)
xlabel('PC used for linear model')
ylabel('Unshuffled R^2')
axis([0 highestPCtouse+1 -1 1])
box off
set(gca,'FontSize',15)
subplot(1,2,2)
distributionPlot(myr2shuffled,'histOpt',1,'colormap',1-gray(64),'showMM',0)
xlabel('PC used for linear model')
ylabel('Shuffled R^2')
axis([0 highestPCtouse+1 -1 1])
box off
set(gca,'FontSize',15)


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
% 
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