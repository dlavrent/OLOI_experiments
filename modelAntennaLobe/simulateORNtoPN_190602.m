%% use hallem and carlson data to model antenna lobe
% analysis based off Luo et al PNAS 2010
clear all
close all

fdir='C:\Users\mac0456\Dropbox\flyimaging\analysis\modelAntennaLobe';
fname='hallemAndCarlson2004_tableS1data.xlsx';
[NUM,TXT,RAW]=xlsread([fdir '\' fname]);

orstoremove=[8 13 16 23];
ORNr=NUM(1:110,:)'+NUM(111,:)';  % add background firing rate for each OR

ORNr(orstoremove,:)=[]; 
TXT(orstoremove)=[];

%ORNr(ORNr<0)=0;


meanF=mean(ORNr,2);  % mean firing rate per OR


nodor=size(ORNr,2); % number of odors
norn=size(ORNr,1); % number of ORNs



% generate PN responses
rmax=165;
sigmap=12;
m=0.05;
for j=1:nodor
    PNr(:,j)=rmax*ORNr(:,j).^(1.5)./(sigmap^(1.5)+ORNr(:,j).^(1.5)+m*sum(ORNr(:,j))^(1.5));
    PNnol(:,j)=rmax*ORNr(:,j).^(1.5)./(sigmap^(1.5)+ORNr(:,j).^(1.5)); % pn firing rate with no lateral inhibition/divisive normalization
end
PNr=real(PNr);
PNnol=real(PNnol);

meanFPN=mean(PNr,2);  % mean firing rate per OR

orntopntransition=nanmean(PNr./ORNr,2)./sum(isfinite(PNr./ORNr),2);
figure;
plot(meanF,orntopntransition,'o')
xlabel('mean firing rate')
ylabel('magnitude of ORN to PN transition')


% plot PN vs ORN firing
figure;
plot(PNr(:),ORNr(:),'o')

figure
subplot(1,2,1)
imagesc(ORNr,[0 300])
set(gca,'YTick',[1:length(meanF)])
set(gca,'YTickLabel',TXT)
ylabel('OR')
xlabel('odorant')
title('ORN firing rate')
h=colorbar;
title(h,'Hz')
set(gca,'FontSize',15)
subplot(1,2,2)
imagesc(PNr,[0 300])
title('PN firing rate')
set(gca,'YTick',[1:length(meanF)])
set(gca,'YTickLabel',TXT)
xlabel('odorant')
h2=colorbar;
title(h2,'Hz')
set(gca,'FontSize',15)

nbins=25;
figure
h1 = histogram(ORNr(:),'Normalization','pdf');
h1.NumBins=nbins;
hold on
h2 = histogram(PNr(:),'Normalization','pdf');
h2.NumBins=nbins;
legend('ORN','PN')
legend boxoff
xlabel('firing rate (hz)')
ylabel('probability')
box off
set(gca,'FontSize',15)

[coeff, score, latent, tsquared, explained] = pca(ORNr'); % ORN pca
[coeffp, scorep, latentp, tsquaredp, explainedp] = pca(PNr'); % PN with lateral inhibition and divisive normalization
[coeffn, scoren, latentn, tsquaredn, explainedn] = pca(PNnol'); % PN without divisive normalization

figure
subplot(1,3,1)
plot(100*explained/sum(explained),'ko','LineWidth',3)
ylabel('variance explained (%)')
axis([0 size(coeff,1) 0 50])
title('ORN')
box off
set(gca,'FontSize',15)
subplot(1,3,2)
plot(100*explainedn/sum(explainedn),'ko','LineWidth',3)
xlabel('pc')
axis([0 size(coeff,1) 0 50])
title('PN (without lateral inhibition)')
box off
set(gca,'FontSize',15)
subplot(1,3,3)
plot(100*explainedp/sum(explainedp),'ko','LineWidth',3)
axis([0 size(coeff,1) 0 50])
title('PN (with lateral inhibition)')
box off
set(gca,'FontSize',15)






%% add noise to each ORN and PN (uniform across all odors)


nflies=100;
zscore=0;

% noise parameters

% ORN odor-specific noise parameters
eta=0;
alpha=0.025;

% ORN odor-independent noise parameters
eta2=10;
alpha2=0.025;

% PN odor-specific noise
eta3=0;
alpha3=0.025;

% PN odor-independent noise parameters
eta4=0;
alpha4 = 0.025;



ornio2=zeros(nflies,norn*nodor); % individual ORN samples
pnio2=zeros(nflies,norn*nodor);  % individual pn samples
pnnonoise=zeros(nflies,norn*nodor);  % individual PN samples without noise
pntemp=zeros(norn,nodor); 
pntempnonoise=zeros(norn,nodor); 




for i=1:nflies
    
    % add ORN odor-specific noise (represents trial variation plus lateral
    % connections)
    temp=ORNr+eta*tanh(alpha*ORNr).*randn(size(ORNr,1),size(ORNr,2));
    tempnonoise=ORNr;
    
    
    % add ORN odor-independent noise (represents ORN->PN synapse variation)
    for j=1:norn
        temp(j,:)=(temp(j,:)+eta2*tanh(alpha2*meanF(j))*randn(1));%*mean(meanF)/eta2;
    end
   
    
    %temp(temp<0)=0;
    %tempnonoise(tempnonoise<0)=0;
    for j=1:nodor
        pntemp(:,j)=rmax*temp(:,j).^(1.5)./(sigmap^(1.5)+temp(:,j).^(1.5)+m*sum(temp(:,j))^(1.5));
        pntempnonoise(:,j)=rmax*tempnonoise(:,j).^(1.5)./(sigmap^(1.5)+tempnonoise(:,j).^(1.5)+m*sum(tempnonoise(:,j))^(1.5));
    end
    
    
    pntemp=real(pntemp);
    pntempnonoise=real(pntempnonoise);
    
    % add PN odor-specific noise  (represents trial variation plus lateral
    % connections)
    pntemp=pntemp+eta3*tanh(alpha3*pntemp).*randn(size(pntemp,1),size(pntemp,2));
    
    for j=1:norn
        % add PN odor-independent noise (represents ORN->PN synapse variation)
         pnio2(i,((j-1)*nodor+1):j*nodor)=(pntemp(j,:)'+eta4*tanh(alpha4*meanFPN(j))*randn(1));%*mean(meanFPN)/eta2;
         pnnonoise(i,((j-1)*nodor+1):j*nodor)=pntempnonoise(j,:)';
         ornio2(i,((j-1)*nodor+1):j*nodor)=temp(j,:)';
    end
    
    %pnio2(pnio2<0)=0;
    %pnnonoise(pnnonoise<0)=0;
end
pnio2=real(pnio2);
pnnonoise=real(pnnonoise);

% z score
if zscore==1
    pnio2=(pnio2-mean(pnio2))./std(pnio2);
    pnio2(isnan(pnio2))=0;
    
    ornio2=(ornio2-mean(ornio2))./std(ornio2);
    ornio2(isnan(ornio2))=0;
end
[coeffo2, scoreo2, latento2, tsquaredo2, explainedo2] = pca(pnio2);

[coefforn2, scoreorn2, latentorn2, tsquaredorn2, explainedorn2] = pca(ornio2);

% look at average absolute loadings for each glomerulus
abscoeffo2=zeros(norn,size(coeffo2,2));
abscoefforn2=zeros(norn,size(coefforn2,2));
for i=1:size(coeffo2,2)
    for j=1:norn
        abscoeffo2(j,i)=mean(abs(coeffo2(((j-1)*nodor+1):j*nodor,i)));
         abscoefforn2(j,i)=mean(abs(coefforn2(((j-1)*nodor+1):j*nodor,i)));
    end
end

varianceRetained=80;
temp=find(cumsum(explainedo2)<varianceRetained);
highestPCo2=temp(end);
temp=find(cumsum(explainedorn2)<varianceRetained);
highestPCorn2=temp(end);

figure;
subplot(1,2,1)
imagesc(pnnonoise,[0 150])
ylabel('individual')
xlabel('dimension')
title('PN firing rate (no noise)')
set(gca,'FontSize',15)
h1=colorbar;
title(h1,'Hz')
subplot(1,2,2)
imagesc(pnio2,[0 150])
xlabel('dimension')
title('PN firing rate (with noise)')
h2=colorbar;
title(h2,'Hz')
set(gca,'FontSize',15)

figure;
plot(100*explainedo2/sum(explainedo2),'ko','LineWidth',2)
xlabel('pc')
ylabel('variance explained (%)')
box off
set(gca,'FontSize',15)

figure;
subplot(1,2,1)
plot(meanF,mean(abscoeffo2(:,1:highestPCo2),2),'ko','LineWidth',2,'MarkerSize',10)
box off
xlabel('mean ORN firing rate')
ylabel('abs(loading weight)')
set(gca,'FontSize',15)
subplot(1,2,2)
plot(meanF,mean(abscoeffo2(:,1:highestPCo2),2)./meanF,'ko','LineWidth',2,'MarkerSize',10)
box off
xlabel('mean ORN firing rate')
ylabel('normalized abs(loading weight)')
set(gca,'FontSize',15)

%% look at differences in total AL output for all odor pairs across individuals
odorpairstd=zeros(nodor,nodor);
for i=1:nodor
    for j=1:nodor
        temp1=mean(pnio2(:,i:nodor:end),2);
        temp2=mean(pnio2(:,j:nodor:end),2);
        tempdiff=temp1-temp2;
        odorpairstd(i,j)=std(tempdiff);
    end
end
[vs inds]=sort(mean(ORNr));

figure;imagesc(odorpairstd(inds,inds))
h=colorbar;
title(h,'\sigma')
xlabel('odorant')
ylabel('odorant')
set(gca,'FontSize',15)

figure; plot(mean(ORNr),mean(odorpairstd),'k*','LineWidth',3)
xlabel('mean odor-evoked firing rate (hz)')
ylabel('mean response difference \sigma') 
box off
set(gca,'FontSize',15)



%same thing for glomeruli
glompairstd=zeros(norn,norn);
for i=1:norn
    for j=1:norn
        temp1=mean(pnio2(:,((i-1)*nodor+1):(i*nodor)),2);
        temp2=mean(pnio2(:,((j-1)*nodor+1):(j*nodor)),2);
        tempdiff=temp1-temp2;
        
        
%         temp1=mean(ornio2(:,((i-1)*nodor+1):(i*nodor)),2);
%         temp2=mean(ornio2(:,((j-1)*nodor+1):(j*nodor)),2);
%         tempdiff=temp1-temp2;
        glompairstd(i,j)=std(tempdiff);
    end
end
[vs inds]=sort(meanF);
figure;imagesc(glompairstd(inds,inds))
set(gca,'XTick',[1:norn])
set(gca,'XTickLabel',TXT(inds))
xtickangle(30)
set(gca,'YTick',[1:norn])
set(gca,'YTickLabel',TXT(inds))
h=colorbar;
title(h,'\sigma')
ylabel('OR')
xlabel('OR')
set(gca,'FontSize',15)

figure; plot(meanF,mean(glompairstd),'k*','LineWidth',3)
xlabel('mean ORN firing rate (hz)')
ylabel('mean response difference \sigma') 
box off
set(gca,'FontSize',15)