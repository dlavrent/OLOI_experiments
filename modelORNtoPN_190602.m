%% model ORN-PN to estimate effects of noise
% based off Luo et al PNAS 2010
clear all
close all

norn=20;
nodor=110;
minFiringRate=25;
maxFiringRate=300;
meanFiringRate=80;

% generate mean firing rates for each ORN
meanF=rand(norn,1)*(maxFiringRate-minFiringRate)+minFiringRate;
meanO=rand(nodor,1)*(maxFiringRate-minFiringRate)+minFiringRate;

noiselevel=3;

ORNr=zeros(norn,nodor);
PNr=zeros(norn,nodor);
% generate response vector for each orn
% for i=1:norn
%     ORNr(i,:)=randn(1,nodor)*meanF(i)/noiselevel+meanF(i);
%     ORNr(i,ORNr(i,:)<minFiringRate)=minFiringRate;
%     ORNr(i,ORNr(i,:)>maxFiringRate)=maxFiringRate;
% end

for i=1:norn
    for j=1:nodor
    ORNr(i,j)=randn(1,1)*meanF(i)*meanO(j)/noiselevel+meanF(i)*meanO(j);
    ORNr(i,ORNr(i,:)<minFiringRate)=minFiringRate;
    end
end
ORNr=meanFiringRate*ORNr/mean(mean(ORNr));
ORNr(ORNr<minFiringRate)=minFiringRate;
% generate PN responses
rmax=165;
sigmap=12;
m=0.05;
for j=1:nodor
PNr(:,j)=rmax*ORNr(:,j).^(1.5)./(sigmap^(1.5)+ORNr(:,j).^(1.5)+m*sum(ORNr(:,j))^(1.5));
end

[coeff, score, latent, tsqu
    
ared, explained] = pca(ORNr');
[coeffp, scorep, latentp, tsquaredp, explainedp] = pca(PNr');

%% now add uniform noise to each ORN->PN synapse

PNrn=zeros(norn,nodor);
% generate PN responses
rmax=165;
sigmap=12;
m=0.05;
ornpnnoiselevel=50;
ornpnnoise=randn(norn,1)*ornpnnoiselevel;
for j=1:nodor
    temp=ORNr(:,j)+ornpnnoise*meanO(j)/mean(meanO);

    temp(temp<minFiringRate)=minFiringRate;
    PNrn(:,j)=rmax*(temp).^(1.5)./(sigmap^(1.5)+(temp).^(1.5)+m*sum(ORNr(:,j))^(1.5));
end

[coeff, score, latent, tsquared, explained] = pca(PNrn');
[coeffp, scorep, latentp, tsquaredp, explainedp] = pca(PNr');

figure
plot(meanF,abs(coeffp(:,1)),'o','LineWidth',2)
hold on
plot(meanF,abs(coeff(:,1)),'o','LineWidth',2)
legend('no noise','noise')

figure
plot(meanF,abs(coeff(:,1))-abs(coeffp(:,1)),'o','LineWidth',2)