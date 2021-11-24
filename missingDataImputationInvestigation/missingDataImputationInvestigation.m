% generate toy data with r normal clusters to investigate the effects of different
% missing data imputation methods
% Matt Churgin
clear all
close all

dims=50; % # dimensions
n=2000; % # observations

sig=[5/65]; % sigma multiplier 
holdout=[0.4]; % fraction of data to hold out when assessing affects of missing data imputation
clusters=[1 10]; % # of clusters in ground truth data

figure

for r = 1:length(clusters)
    data=[];
    for i=1:clusters(r)
        mu=randn(1,dims);
        cov=randn(dims,dims);
        cov=sig*cov'*cov;
        
        data=[data; mvnrnd(mu,cov,n/clusters(r))];
    end
    
    % pca on original data
    [COEFF SCORE LATENT TSQUARED EXPLAINED MU] = pca(data);
    
    % hold out data at random
    temp=rand(1,(dims*n));
    holdoutinds=temp<=holdout;
    
    dataDecimated=data;
    dataDecimated(holdoutinds)=NaN;
    
    % mean fill missing data
    dataDecimatedMeanfilled=dataDecimated;
    for i=1:dims
        temp=isnan(dataDecimated(:,i));
        dataDecimatedMeanfilled(temp,i)=nanmean(dataDecimated(:,i));
    end
    
    % pca on mean-filled data
    [COEFFMF SCOREMF LATENTMF TSQUAREDMF EXPLAINEDMF MUMF] = pca(dataDecimatedMeanfilled);
    
    % pca on decimated data using ALS
    [COEFFALS SCOREALS LATENTALS TSQUAREDALS EXPLAINEDALS MUALS] = pca(dataDecimated,'algorithm','als');
    
    % pca on original data using ALS
    [COEFFALSO SCOREALSO LATENTALSO TSQUAREDALSO EXPLAINEDALSO MUALSO] = pca(data,'algorithm','als');
    
    
    % reconstruct each matrix 
    
%     origReconstructed=SCORE(:,1:clusters(r))*COEFF(:,1:clusters(r))' + MU;
%     mfReconstructed=SCOREMF(:,1:clusters(r))*COEFFMF(:,1:clusters(r))' + MUMF;
%     alsReconstructed=SCOREALS(:,1:clusters(r))*COEFFALS(:,1:clusters(r))' + MUALS;
%     alsoReconstructed=SCOREALSO(:,1:clusters(r))*COEFFALSO(:,1:clusters(r))' + MUALSO;
%     
    origReconstructed=SCORE*COEFF' + MU;
    mfReconstructed=SCOREMF*COEFFMF' + MUMF;
    alsReconstructed=SCOREALS*COEFFALS' + MUALS;
    alsoReconstructed=SCOREALSO*COEFFALSO' + MUALSO;
    
    origError(r)=nansum(nansum((data(holdoutinds)-origReconstructed(holdoutinds)).^2));
    mfError(r)=nansum(nansum((data(holdoutinds)-mfReconstructed(holdoutinds)).^2));
    alsError(r)=nansum(nansum((data(holdoutinds)-alsReconstructed(holdoutinds)).^2));
    alsoError(r)=nansum(nansum((data(holdoutinds)-alsoReconstructed(holdoutinds)).^2));
    
    ocolor=[0 0 0];
    mfcolor=[0.8 0.4 0.8];
    alscolor=[0.4 0.4 0.8];
    
    subplot(1,length(clusters),r)
    plot(EXPLAINED/sum(EXPLAINED),'Color',ocolor,'LineWidth',3)
    hold on
    plot(EXPLAINEDMF/sum(EXPLAINEDMF),'o','Color',mfcolor,'LineWidth',3)
    plot(EXPLAINEDALS/sum(EXPLAINEDALS),'--','Color',alscolor,'LineWidth',3)
    plot(EXPLAINEDALSO/sum(EXPLAINEDALSO),'x-','Color',ocolor,'LineWidth',3)
    legend('Original','Mean-filled','ALS','Original + ALS')
    title([num2str(clusters(r)) ' clusters'])
    xlabel('PC #')
    ylabel('var. expl.')
    box off
    legend boxoff
    set(gca,'FontSize',15)
    drawnow
    
end

% plot reconstruction error
figure;
plot(clusters,origError,'o','LineWidth',2,'MarkerSize',10)
hold on
plot(clusters,mfError,'x','LineWidth',2,'MarkerSize',10)
plot(clusters,alsError,'s','LineWidth',2,'MarkerSize',10)
plot(clusters,alsoError,'*','LineWidth',2,'MarkerSize',10)
legend('Original','Mean-filled','ALS','Original + ALS')
legend boxoff
box off
ylabel('mean squared reconstruction error')
xlabel('ground truth clusters')
set(gca,'FontSize',15)