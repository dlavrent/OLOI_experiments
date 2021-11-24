%% extract data from chou paper, nature neuroscience 2010
% figure 7a
clear all
close all

a=imread('ChouFig7.png');
b=rgb2gray(a);
cropxy=[1127 1391; 1432 1920];
c=b(cropxy(1,2):cropxy(2,2),cropxy(1,1):cropxy(2,1));

thresh=150;

d=double(c>thresh);

numGloms=54;

xwidth=abs(1127-1432);
xstep=xwidth/numGloms;

nCells=161;
ywidth=abs(1390-1921);
ystep=ywidth/nCells;

% sample from cropped image to get original data
originalData=zeros(nCells,numGloms);
for i=1:numGloms
    for j=1:nCells
        originalData(j,i)=d(round(ystep/2+(j-1)*ystep),round(xstep/2+(i-1)*xstep));
    end
end

% for each glomerulus pair, count how many times they are innervated by the
% same cell

pairedInnervation=zeros(numGloms,numGloms);
for i=1:numGloms
    for j=1:numGloms
        if i~=j
            pairedInnervation(i,j)=sum(originalData(:,i).*originalData(:,j));
        end
    end
end


% create randomly permuted matrix
shuffles=1000;
randConvergences=zeros(numGloms,numGloms,shuffles);
for k=1:shuffles
    temp=zeros(numGloms,numGloms);
    for i=1:numGloms
        for j=1:numGloms
            if j~=i
                temp(i,j)=sum(originalData(randperm(size(originalData,1)),i).*originalData(randperm(size(originalData,1)),j));
            end
        end
    end
    randConvergences(:,:,k)=temp;
    if mod(k,100)==0
    disp(['iter ' num2str(k) ' of ' num2str(shuffles)])
    end
end
randConvergencesMean=mean(randConvergences,3);
convergencesAboveExpectation=(pairedInnervation-randConvergencesMean);
figure;
hist(convergencesAboveExpectation(:),100)
xlabel('observed convergences minus shuffled convergences')
ylabel('number of instances')

figure;
subplot(1,3,1)
imagesc(pairedInnervation)
title('Observed')
% set(gca,'ytick',1:length(glomNames),'yticklabel',string(glomNames),'FontSize',10)
% set(gca,'xtick',1:length(glomNames),'xticklabel',string(glomNames),'FontSize',10)
xtickangle(30)
subplot(1,3,2)
imagesc(randConvergencesMean)
title('Shuffled')
% set(gca,'ytick',1:length(glomNames),'yticklabel',string(glomNames),'FontSize',10)
% set(gca,'xtick',1:length(glomNames),'xticklabel',string(glomNames),'FontSize',10)
xtickangle(30)
subplot(1,3,3)
imagesc(convergencesAboveExpectation)
title('Difference')
% set(gca,'ytick',1:length(glomNames),'yticklabel',string(glomNames),'FontSize',10)
% set(gca,'xtick',1:length(glomNames),'xticklabel',string(glomNames),'FontSize',10)
xtickangle(30)

%% fig 3a

clear all
close all

a=imread('ChouFig3a.png');
b=rgb2gray(a);
cropxy=[1042 941; 1516 2371];
c=b(cropxy(1,2):cropxy(2,2),cropxy(1,1):cropxy(2,1));

thresh=150;

d=double(c>thresh);

numGloms=54;

xwidth=abs(1042-1515);
xstep=xwidth/numGloms;

nCells=1431; % real number is 1532
ywidth=abs(942-2372);
ystep=ywidth/nCells;

% sample from cropped image to get original data
originalData=zeros(nCells,numGloms);
for i=1:numGloms
    for j=1:nCells
        originalData(j,i)=d(round(ystep/2+(j-1)*ystep+0.1),round(xstep/2+(i-1)*xstep));
    end
end

% for each glomerulus pair, count how many times they are innervated by the
% same cell

pairedInnervation=zeros(numGloms,numGloms);
for i=1:numGloms
    for j=1:numGloms
        if i~=j
            pairedInnervation(i,j)=sum(originalData(:,i).*originalData(:,j));
        end
    end
end


% create randomly permuted matrix
shuffles=100;
randConvergences=zeros(numGloms,numGloms,shuffles);
for k=1:shuffles
    temp=zeros(numGloms,numGloms);
    for i=1:numGloms
        for j=1:numGloms
            if j~=i
                temp(i,j)=sum(originalData(randperm(size(originalData,1)),i).*originalData(randperm(size(originalData,1)),j));
            end
        end
    end
    randConvergences(:,:,k)=temp;
    if mod(k,10)==0
    disp(['iter ' num2str(k) ' of ' num2str(shuffles)])
    end
end
randConvergencesMean=mean(randConvergences,3);
convergencesAboveExpectation=(pairedInnervation-randConvergencesMean);
figure;
hist(convergencesAboveExpectation(:),100)
xlabel('observed convergences minus shuffled convergences')
ylabel('number of instances')

figure;
subplot(1,3,1)
imagesc(pairedInnervation)
title('Observed')
% set(gca,'ytick',1:length(glomNames),'yticklabel',string(glomNames),'FontSize',10)
% set(gca,'xtick',1:length(glomNames),'xticklabel',string(glomNames),'FontSize',10)
xtickangle(30)
subplot(1,3,2)
imagesc(randConvergencesMean)
title('Shuffled')
% set(gca,'ytick',1:length(glomNames),'yticklabel',string(glomNames),'FontSize',10)
% set(gca,'xtick',1:length(glomNames),'xticklabel',string(glomNames),'FontSize',10)
xtickangle(30)
subplot(1,3,3)
imagesc(convergencesAboveExpectation)
title('Difference')
% set(gca,'ytick',1:length(glomNames),'yticklabel',string(glomNames),'FontSize',10)
% set(gca,'xtick',1:length(glomNames),'xticklabel',string(glomNames),'FontSize',10)
xtickangle(30)

%% supp fig 6a. Line 1.H84-93 gal4

clear all
close all

a=imread('suppfig6a.png');
b=rgb2gray(a);
cropxy=[561   329;  716 1325];
c=b(cropxy(1,2):cropxy(2,2),cropxy(1,1):cropxy(2,1));

thresh=160;

d=double(c>thresh);

numGloms=54;

xwidth=abs(cropxy(1,1)-cropxy(2,1));
xstep=xwidth/numGloms;

nCells=362; 
ywidth=abs(cropxy(1,2)-cropxy(2,2));
ystep=ywidth/nCells;

% sample from cropped image to get original data
originalData=zeros(nCells,numGloms);
for i=1:numGloms
    for j=1:nCells
        originalData(j,i)=d(round(ystep/2+(j-1)*ystep+0.1),round(xstep/2+(i-1)*xstep));
    end
end

% for each glomerulus pair, count how many times they are innervated by the
% same cell

pairedInnervation=zeros(numGloms,numGloms);
for i=1:numGloms
    for j=1:numGloms
        if i~=j
            pairedInnervation(i,j)=sum(originalData(:,i).*originalData(:,j));
        end
    end
end


% create randomly permuted matrix
shuffles=100;
randConvergences=zeros(numGloms,numGloms,shuffles);
for k=1:shuffles
    temp=zeros(numGloms,numGloms);
    for i=1:numGloms
        for j=1:numGloms
            if j~=i
                temp(i,j)=sum(originalData(randperm(size(originalData,1)),i).*originalData(randperm(size(originalData,1)),j));
            end
        end
    end
    randConvergences(:,:,k)=temp;
    if mod(k,10)==0
    disp(['iter ' num2str(k) ' of ' num2str(shuffles)])
    end
end
randConvergencesMean=mean(randConvergences,3);
convergencesAboveExpectation=(pairedInnervation-randConvergencesMean);
figure;
hist(convergencesAboveExpectation(:),100)
xlabel('observed convergences minus shuffled convergences')
ylabel('number of instances')

figure;
subplot(1,3,1)
imagesc(pairedInnervation)
title('Observed')
% set(gca,'ytick',1:length(glomNames),'yticklabel',string(glomNames),'FontSize',10)
% set(gca,'xtick',1:length(glomNames),'xticklabel',string(glomNames),'FontSize',10)
xtickangle(30)
subplot(1,3,2)
imagesc(randConvergencesMean)
title('Shuffled')
% set(gca,'ytick',1:length(glomNames),'yticklabel',string(glomNames),'FontSize',10)
% set(gca,'xtick',1:length(glomNames),'xticklabel',string(glomNames),'FontSize',10)
xtickangle(30)
subplot(1,3,3)
imagesc(convergencesAboveExpectation)
title('Difference')
% set(gca,'ytick',1:length(glomNames),'yticklabel',string(glomNames),'FontSize',10)
% set(gca,'xtick',1:length(glomNames),'xticklabel',string(glomNames),'FontSize',10)
xtickangle(30)


%% supp fig 6e. Line 5. NP3056-gal4

clear all
close all

a=imread('suppfig6e.png');
b=rgb2gray(a);
cropxy=[558 323 ; 712 1324 ];
c=b(cropxy(1,2):cropxy(2,2),cropxy(1,1):cropxy(2,1));

thresh=130;

d=double(c>thresh);

numGloms=54;

xwidth=abs(558-712);
xstep=xwidth/numGloms;

nCells=578; 
ywidth=abs(323-1324);
ystep=ywidth/nCells;

% sample from cropped image to get original data
originalData=zeros(nCells,numGloms);
for i=1:numGloms
    for j=1:nCells
        originalData(j,i)=d(round(ystep/2+(j-1)*ystep+0.1),round(xstep/2+(i-1)*xstep));
    end
end

% for each glomerulus pair, count how many times they are innervated by the
% same cell

pairedInnervation=zeros(numGloms,numGloms);
for i=1:numGloms
    for j=1:numGloms
        if i~=j
            pairedInnervation(i,j)=sum(originalData(:,i).*originalData(:,j));
        end
    end
end


% create randomly permuted matrix
shuffles=100;
randConvergences=zeros(numGloms,numGloms,shuffles);
for k=1:shuffles
    temp=zeros(numGloms,numGloms);
    for i=1:numGloms
        for j=1:numGloms
            if j~=i
                temp(i,j)=sum(originalData(randperm(size(originalData,1)),i).*originalData(randperm(size(originalData,1)),j));
            end
        end
    end
    randConvergences(:,:,k)=temp;
    if mod(k,10)==0
    disp(['iter ' num2str(k) ' of ' num2str(shuffles)])
    end
end
randConvergencesMean=mean(randConvergences,3);
convergencesAboveExpectation=(pairedInnervation-randConvergencesMean);
figure;
hist(convergencesAboveExpectation(:),100)
xlabel('observed convergences minus shuffled convergences')
ylabel('number of instances')

figure;
subplot(1,3,1)
imagesc(pairedInnervation)
title('Observed')
% set(gca,'ytick',1:length(glomNames),'yticklabel',string(glomNames),'FontSize',10)
% set(gca,'xtick',1:length(glomNames),'xticklabel',string(glomNames),'FontSize',10)
xtickangle(30)
subplot(1,3,2)
imagesc(randConvergencesMean)
title('Shuffled')
% set(gca,'ytick',1:length(glomNames),'yticklabel',string(glomNames),'FontSize',10)
% set(gca,'xtick',1:length(glomNames),'xticklabel',string(glomNames),'FontSize',10)
xtickangle(30)
subplot(1,3,3)
imagesc(convergencesAboveExpectation)
title('Difference')
% set(gca,'ytick',1:length(glomNames),'yticklabel',string(glomNames),'FontSize',10)
% set(gca,'xtick',1:length(glomNames),'xticklabel',string(glomNames),'FontSize',10)
xtickangle(30)
