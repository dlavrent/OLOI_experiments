clear all
close all
[num txt raw]=xlsread('JeanneWilson2018_Table2');
glomNames=txt(2,4:end);
% remove widely connected lines: local 1, local 3, and ML
rowstoremove=[1:3 6 27:51];

numPruned=num;
numPruned(rowstoremove,:)=[];

numBinarized=numPruned>0;

% find # convergences for each pair of gloms
convergences=zeros(39,39);
for i=1:39
    for j=1:39
        if j~=i
           convergences(i,j)=sum(numBinarized(:,i).*numBinarized(:,j)); 
        end
    end
end

% create randomly permuted matrix
shuffles=10000;
randConvergences=zeros(39,39,shuffles);
for k=1:shuffles
    temp=zeros(39,39);
    for i=1:39
        for j=1:39
            if j~=i
                temp(i,j)=sum(numBinarized(randperm(size(numBinarized,1)),i).*numBinarized(randperm(size(numBinarized,1)),j));
            end
        end
    end
    randConvergences(:,:,k)=temp;
    if mod(k,100)==0
    disp(['iter ' num2str(k) ' of ' num2str(shuffles)])
    end
end
randConvergencesMean=mean(randConvergences,3);

figure;
subplot(1,3,1)
imagesc(convergences)
title('Observed')
set(gca,'ytick',1:length(glomNames),'yticklabel',string(glomNames),'FontSize',10)
set(gca,'xtick',1:length(glomNames),'xticklabel',string(glomNames),'FontSize',10)
xtickangle(30)
subplot(1,3,2)
imagesc(randConvergencesMean)
title('Shuffled')
set(gca,'ytick',1:length(glomNames),'yticklabel',string(glomNames),'FontSize',10)
set(gca,'xtick',1:length(glomNames),'xticklabel',string(glomNames),'FontSize',10)
xtickangle(30)
subplot(1,3,3)
imagesc((convergences-randConvergencesMean))
title('Difference')
set(gca,'ytick',1:length(glomNames),'yticklabel',string(glomNames),'FontSize',10)
set(gca,'xtick',1:length(glomNames),'xticklabel',string(glomNames),'FontSize',10)
xtickangle(30)

figure;
imagesc((convergences-randConvergencesMean)>1)
title('LHN Convergences Above Expectation')
set(gca,'ytick',1:length(glomNames),'yticklabel',string(glomNames),'FontSize',10)
set(gca,'xtick',1:length(glomNames),'xticklabel',string(glomNames),'FontSize',10)
xtickangle(30)