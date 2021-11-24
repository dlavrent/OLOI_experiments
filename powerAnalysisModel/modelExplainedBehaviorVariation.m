% model behavioral variance explained as a function of behavioral persistence
% and true correlation between latent variable of interest with behavior
% Matt Churgin, de Bivort lab 2019
clear all
close all

numflies=50; % number of flies
iters=100; % number of iterations
nsteps=40; % num steps for x-y axis. 40 is a good chioce

behaviorcorr=[0:1/nsteps:1];
trueLatentCorr=[0:1/nsteps:1];

% behavior is assumed to come from a normal distribution
meanBehavior=0.4; % population mean behavioral measurement
standarddevBehavior=0.1; % population behavioral standard deviation 

varianceExplained=zeros(iters,length(behaviorcorr),length(trueLatentCorr));
corrc=zeros(iters,length(behaviorcorr),length(trueLatentCorr));

for i=1:iters
    flyData{i}=zeros(numflies,3,length(behaviorcorr),length(trueLatentCorr));
    for k=1:length(behaviorcorr)
        for kk=1:length(trueLatentCorr)
            
            firstrandomnumber=standarddevBehavior*randn(numflies,1);
            secondrandomnumber=standarddevBehavior*randn(numflies,1);
            thirdrandomnumber=standarddevBehavior*randn(numflies,1);
            
            % generate each fly's behavior measurement
            measuredbehavior=meanBehavior+firstrandomnumber*standarddevBehavior;
            
            % generate expected fly behavior at time of latent predictor measurement (dependent on
            % behavioral persistence correlation value, behaviorcorr(k))
            expectedBehavior=meanBehavior+sqrt(behaviorcorr(k))*firstrandomnumber+sqrt(1-sqrt(behaviorcorr(k))^2)*secondrandomnumber;
            
            % generate ideal behavior predicted from imaging (dependent on
            % true correlation between latent predictor and behavior, trueLatentCorr(kk))
            behaviorPredictionIdeal=meanBehavior+(expectedBehavior-meanBehavior)*sqrt(trueLatentCorr(kk))+thirdrandomnumber*(sqrt(1-sqrt(trueLatentCorr(kk))^2));
            
            flyData{i}(:,:,k,kk)=[measuredbehavior expectedBehavior behaviorPredictionIdeal];
            
            % predict measuredbehavior from behaviorPredictionIdeal
            linmodel=fitlm(flyData{i}(:,3,k,kk),flyData{i}(:,1,k,kk));
            
            % save R^2
            varianceExplained(i,k,kk)=linmodel.Rsquared.Ordinary;
            
            % save corr coefficient
            [r p]=corrcoef(flyData{i}(:,3,k,kk),flyData{i}(:,1,k,kk));
            corrc(i,k,kk)=r(1,2);
        end
    end
    disp(['iteration ' num2str(i) ' of ' num2str(iters)])
end
finalVarianceExplained=squeeze(mean(varianceExplained,1));
finalCorrc=squeeze(mean(corrc,1));
disp('done :)')

save neuralactivitybehavior_model
%%
clear all 
close all
load neuralactivitybehavior_model

figure
imagesc(finalVarianceExplained)
hold on
[C,h] = contour(finalVarianceExplained,[0.1 0.2 0.3 0.4 0.5]);
clabel(C,h,'FontSize',15)
h.LineWidth = 3;
h.LineColor='k';
h.LineStyle=':';
ylabel('behavioral persistence r')
xlabel('neural activity-behavior r')
hcb=colorbar;
title(hcb,'R^2')
set(gca,'XTick',[1:round(size(finalVarianceExplained,2)/10):size(finalVarianceExplained,2)])
set(gca,'YTick',[1:round(size(finalVarianceExplained,1)/10):size(finalVarianceExplained,1)])
set(gca,'XTickLabel',[trueLatentCorr(1:round(length(trueLatentCorr)/10):end)])
set(gca,'YTickLabel',[behaviorcorr(1:round(length(behaviorcorr)/10):end)])
set(gca,'FontSize',15)

%% plot latent variable - behavior correlation as heat map
% y variable is persistence, x variable is measured r-squared
rsquaredbins=20;
latentVariableCorr=NaN*zeros(length(behaviorcorr),rsquaredbins);
for i=1:length(behaviorcorr)
    for j=1:rsquaredbins
        temp=find(finalVarianceExplained(i,:)<(j/rsquaredbins));
        temp2=find(finalVarianceExplained(i,:)>((j-1)/rsquaredbins));
        temp3=intersect(temp,temp2);
        if ~isempty(trueLatentCorr(temp3))
            latentVariableCorr(i,j)=min(trueLatentCorr(temp3));
        end
    end
end


figure
imagesc(latentVariableCorr)
hold on
[C,h] = contour(latentVariableCorr,[0.5 0.75 0.9]);
clabel(C,h,'FontSize',15)
h.LineWidth = 3;
h.LineColor='k';
h.LineStyle=':';
set(gca,'XTick',[1:round(size(latentVariableCorr,2)/10):size(latentVariableCorr,2)])
set(gca,'YTick',[1:round(size(finalVarianceExplained,1)/10):size(finalVarianceExplained,1)])
set(gca,'XTickLabel',[0:round(size(latentVariableCorr,2)/10):(size(latentVariableCorr,2)-1)]/rsquaredbins)
set(gca,'YTickLabel',[behaviorcorr(1:round(length(behaviorcorr)/10):end)])
ylabel('behavior persistence R^2')
xlabel('measured neural activity-behavior R^2')
%hcb=colorbar;
%title(hcb,'true variance of behavior explained by neural activity')
set(gca,'FontSize',15)


g2=latentVariableCorr;
g2(isnan(g2))=-1;
minc = -0.01;
g2=g2-minc;
maxc = 1;
g2=255*g2/maxc;
g2=g2+1;
vectorPixels(g2,parula(256),[0 0 0])

figure
imagesc(zeros(size(latentVariableCorr)))
hold on
[C,h] = contour(latentVariableCorr,[0.5 0.75 0.9]);
clabel(C,h,'FontSize',15)
h.LineWidth = 3;
h.LineColor='k';
h.LineStyle=':';
