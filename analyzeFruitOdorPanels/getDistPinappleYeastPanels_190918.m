clear all
close all
realhomedir = '/Users/mattchurgin/Dropbox (Harvard University)/flyimaging/analysis/analyzeFruitOdorPanels';

homedir = '/Users/mattchurgin/Dropbox (Harvard University)/flyimaging/Yeast7_Pineapple8_OdorPanels';
cd(homedir)

odornames=[{'air'},{'3-octanol'},{'1-hexanol'},{'ethyl lactate'},{'citronella'},{'2-heptanone'},{'1-pentanol'},{'yeast'},{'pineapple'},{'hexyl acetate'},{'4-methylcyclohexanol'},{'pentyl acetate'},{'1-butanol'}];

folds = dir(pwd); folds=folds(3:end);

averaged = zeros(13,13);
nvols = 0;
nflies = 0;
for i = 1:length(folds) % day
    cd(folds(i).name)
    folds2 = dir(pwd); folds2=folds2(3:end);
    
    for j = 1:length(folds2) % fly
        cd(folds2(j).name)
        currd = zeros(13,13);
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
                    
                    
                    d = squareform(pdist(grnResponse,'cosine'));
                    %imagesc(d)
                    %pause
                    currd=currd+d;
                    nvols = nvols+1;
                    currvols=currvols+1;
                end
            end
            cd ..
        end
        cd ..
    end
    averaged = averaged + currd/currvols;
    cd ..
end
cd(realhomedir)

eyeunwrapped = eye(13);

averaged=averaged/nflies;
averaged(find(eyeunwrapped(:)))=NaN;

figure;
imagesc(averaged)
set(gca,'XTick',[1:length(odornames)])
set(gca,'XTickLabel',odornames)
xtickangle(30)
set(gca,'YTick',[1:length(odornames)])
set(gca,'YTickLabel',odornames)
set(gca,'FontSize',15)
hcb = colorbar;
title(hcb,'distance')

%% sort matrix
averaged(find(eyeunwrapped(:)))=0;
tree = linkage(averaged,'average');
leafOrder = optimalleaforder(tree,averaged);

averaged(find(eyeunwrapped(:)))=NaN;
figure;
imagesc(averaged(leafOrder,leafOrder))
set(gca,'XTick',[1:length(odornames)])
set(gca,'XTickLabel',odornames(leafOrder))
xtickangle(30)
set(gca,'YTick',[1:length(odornames)])
set(gca,'YTickLabel',odornames(leafOrder))
set(gca,'FontSize',15)
hcb = colorbar;
title(hcb,'distance')
%% compare distance in coding space to behavior persistence
% Behavior Persistence scores

combos = [{'yeast-pineapple'}];
persistence = [0]; 
codingdist = [averaged(8,9)];

figure
c = hsv(1);
scatter(codingdist, persistence,100,c,'filled')
text(codingdist,persistence+.01,combos,'FontSize',15)
hold on
scatter(averaged(2,11),0.4,100,'k','filled')
text(averaged(2,11),0.4+.01,'OCT-MCH','FontSize',15)
xlabel('coding distance')
ylabel('behavioral persistence (r)')
set(gca,'FontSize',15)
box off
axis([0 3.2 0 0.5])

%% run pca on one volume
clear all
close all
fname ='fly4/leftLobe/clusterResponses_Volumes.mat';
odornames=[{'air'},{'3-octanol'},{'1-hexanol'},{'ethyl lactate'},{'citronella'},{'2-heptanone'},{'1-pentanol'},{'yeast'},{'pineapple'},{'hexyl acetate'},{'4-methylcyclohexanol'},{'pentyl acetate'},{'1-butanol'}];
file =['/Users/mattchurgin/Dropbox (Harvard University)/flyimaging/Yeast7_Pineapple8_OdorPanels/190918_gh146gcamp_nobehavior_age12days/' fname];
load(file)
%
%grnr = reshape(grnResponse,13,length(clusterVolU)*20);
grnr = mean(grnResponse(:,:,6:10),3);
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(grnr);

cm = hsv(13);
figure;
scatter(SCORE(:,1),SCORE(:,2),100,cm,'filled');text(SCORE(:,1)+0.05,SCORE(:,2)+0.03,odornames,'FontSize',15)
xlabel(['PC 1 (' num2str(round(EXPLAINED(1))) '%)'])
ylabel(['PC 2 (' num2str(round(EXPLAINED(2))) '%)'])
box off
set(gca,'FontSize',15)

