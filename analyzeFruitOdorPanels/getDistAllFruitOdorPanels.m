clear all
close all
realhomedir = '/Users/mattchurgin/Dropbox (Harvard University)/flyimaging/analysis/analyzeFruitOdorPanels';

imagehomedir = '/Users/mattchurgin/Dropbox (Harvard University)/flyimaging/allFruitOdorPanels';
cd(imagehomedir)

odornames=[{'air'},{'3-octanol'},{'1-hexanol'},{'ethyl lactate'},{'N/A'},{'2-heptanone'},{'1-pentanol'},{'N/A'},{'N/A'},{'hexyl acetate'},{'4-methylcyclohexanol'},{'pentyl acetate'},{'1-butanol'}];

folds = dir(pwd); folds=folds(3:end);

averaged = zeros(13,13);
averagedCos = zeros(13,13);
nvols = 0;
nflies = 0;
for i = 1:length(folds) % day
    cd(folds(i).name)
    folds2 = dir(pwd); folds2=folds2(3:end);
    
    for j = 1:length(folds2) % fly
        cd(folds2(j).name)
        currd = zeros(13,13);
        currdcos = zeros(13,13);
        nflies = nflies+1;
        currvols=0;
        folds3 = dir(pwd); folds3=folds3(3:end);
        for k =1:length(folds3) % lobes
            cd(folds3(k).name)
            files = dir(pwd);
            if length(files)~=2
                files=files(3:end);
                for l = 1:length(files)
                    load(files(l).name)
                    
                    grnResponse = reshape(grnResponse,13,length(clusterVolU)*20);
                    
                    %grnResponse = mean(grnResponse(:,:,6:10),3);
                    
                    
                    d = squareform(pdist(grnResponse,'euclidean'));
                    dcos = squareform(pdist(grnResponse,'cosine'));
                    %imagesc(d)
                    %pause
                    currd=currd+d;
                    currdcos=currdcos+dcos;
                    nvols = nvols+1;
                    currvols=currvols+1;
                end
            end
            cd ..
        end
        cd ..
    end
    averaged = averaged + currd/currvols;
    averagedCos = averagedCos + currdcos/currvols;
    cd ..
end
cd(realhomedir)

eyeunwrapped = eye(13);

averaged=averaged/nflies;
averagedCos=averagedCos/nflies;
averaged(find(eyeunwrapped(:)))=NaN;
averagedCos(find(eyeunwrapped(:)))=NaN;

toremove =[5 8 9];
averaged([toremove],:)=NaN;
averaged(:,[toremove])=NaN;
averagedCos([toremove],:)=NaN;
averagedCos(:,[toremove])=NaN;

figure;
imagesc(averaged)
set(gca,'XTick',[1:length(odornames)])
set(gca,'XTickLabel',odornames)
xtickangle(30)
set(gca,'YTick',[1:length(odornames)])
set(gca,'YTickLabel',odornames)
set(gca,'FontSize',15)
hcb = colorbar;
title(hcb,'euclidean distance')

figure;
imagesc(averagedCos)
set(gca,'XTick',[1:length(odornames)])
set(gca,'XTickLabel',odornames)
xtickangle(30)
set(gca,'YTick',[1:length(odornames)])
set(gca,'YTickLabel',odornames)
set(gca,'FontSize',15)
hcb = colorbar;
title(hcb,'cosine distance')

figure
subplot(1,2,1)
histogram(averaged(:),'Normalization','probability')
hold on
plot(ones(1,5)*averaged(2,11),linspace(0,0.3,5),'k--','LineWidth',3)
plot(ones(1,5)*averaged(10,13),linspace(0,0.3,5),'r--','LineWidth',3)
xlabel('euclidean distance')
ylabel('probability')
box off
set(gca,'FontSize',15)
subplot(1,2,2)
histogram(averagedCos(:),'Normalization','probability')
hold on
plot(ones(1,5)*averagedCos(2,11),linspace(0,0.3,5),'k--','LineWidth',3)
plot(ones(1,5)*averagedCos(10,13),linspace(0,0.3,5),'r--','LineWidth',3)
xlabel('cosine distance')
ylabel('probability')
box off
set(gca,'FontSize',15)
%% 

plot(averaged(2,11),averagedCos(2,11),'x','Color',[0.2 0.3 1],'LineWidth',5,'MarkerSize',15)
hold on
plot(averaged(10,13),averagedCos(10,13),'rx','LineWidth',5,'MarkerSize',15)
plot(averaged(1,:),averagedCos(1,:),'x','Color',[0 0.8 0],'LineWidth',5,'MarkerSize',15)
plot(averaged(:),averagedCos(:),'ko')

xlabel('euclidean distance')
ylabel('cosine distance')
legend('OCT vs. MCH','Hexyl Acetate vs. 1-butanol','all vs. air')
legend boxoff
box off
set(gca,'FontSize',15)