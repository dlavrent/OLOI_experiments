%% batch process apply mask to spot images
clear all
close all

load ORN_PN_colors
load analysis_dir_path
homedir=fullfile(analysis_dir_path, 'IHC/expansion_experiment/191120_or13a_brpshort_behaviorAndExpansionImaging_age10days');
spotfolder = 'maskedSpotImages_lowthreshmask';
maskfolder ='masks_lowthresh';
rawimagemaskedfolder='rawimages_masked';
behaviorfolder='behavior';
imagedir = [homedir '/' spotfolder];
maskdir = [homedir '/' maskfolder];
rawimagedir = [homedir '/' rawimagemaskedfolder];
behaviordir = [homedir '/' behaviorfolder];

currfiles = dir(imagedir);
currfiles=currfiles(3:end);
maskedSpots=cell(1,length(currfiles));
masks=cell(1,length(currfiles));
rawimagesmasked=cell(1,length(currfiles));
maskedSpotNames=cell(1,length(currfiles));
behaviorOcc=zeros(1,length(currfiles));
for i = 1:length(currfiles)
    load([imagedir '/' currfiles(i).name])
    
    maskname=currfiles(i).name(1:(end-28));
    load([maskdir '/' maskname '_brpshort_mask_lowthresh.mat'])
    
    load([rawimagedir '/' maskname '__rawimagemasked.mat'])
    
    flyname=currfiles(i).name(1:5);
    load([behaviordir '/' flyname '_behavior.mat'])
    
    maskedSpots{i}=maskedSpotIm;
    maskedSpotNames{i}=currfiles(i).name;
    masks{i}=mask;
    rawimagesmasked{i}=rawim;
    behaviorOcc(i)=occ-preocc;
    disp(['loaded stack ' num2str(i) ' of ' num2str(length(currfiles))])
end

disp('done loading')

%% get mask volume and binarize spots
maskVols=zeros(1,length(maskedSpots));
totalF=zeros(1,length(maskedSpots));
binarizedSpots=cell(1,length(maskedSpots));

for i =1:length(maskedSpots)
    maskVols(i)=sum(masks{i}(:))*0.1*0.1*0.54; % um^3
    binarizedSpots{i}=maskedSpots{i}>0;
    %binarizedSpots{i}=maskedSpots{i}>15;
    
    totalF(i)=sum(sum(sum(rawimagesmasked{i})));
    disp(['processed stack ' num2str(i) ' of ' num2str(length(currfiles))])
end
%%  get spot number and properties
spotStats=cell(1,length(maskedSpots));
for i =1:length(maskedSpots)
   stats = regionprops3(binarizedSpots{i},'volume');
   spotStats{i}=stats;
    disp(['processed spot properties ' num2str(i) ' of ' num2str(length(currfiles))])
end
%% get total spot volume and total number of spots
spotTotal=zeros(1,length(maskedSpots));
totalFDensity=zeros(1,length(maskedSpots));
spotDensity=zeros(1,length(maskedSpots));
spotMedianSize=zeros(1,length(maskedSpots));
for i =1:length(maskedSpots)
    spotTotal(i)=size(spotStats{i},1);
    spotDensity(i)=spotTotal(i)/maskVols(i);
    totalFDensity(i)=totalF(i)/maskVols(i);
    spotMedianSize(i)=prctile(spotStats{i}.Volume,50)*0.1*0.1*0.54;
    
%     vthresh=10;
%     spotTotal(i)=length(find(spotStats{i}.Volume>vthresh));
%     spotDensity(i)=spotTotal(i)/maskVols(i);
%     spotMedianSize(i)=prctile(spotStats{i}.Volume(find(spotStats{i}.Volume>vthresh)),50)*0.1*0.1*0.54;
    disp(['processed spot totals ' num2str(i) ' of ' num2str(length(currfiles))])
end

%% average across left and right lobe
b = zeros(1,length(behaviorOcc)/2);
tf = zeros(1,length(behaviorOcc)/2); % total fluorescence
tfd = zeros(1,length(behaviorOcc)/2); % total fluorescence density
st = zeros(1,length(behaviorOcc)/2); % spot total #
sd = zeros(1,length(behaviorOcc)/2); % spot density (total spots/mask volume)
sm = zeros(1,length(behaviorOcc)/2); % spot median size
for i = 2:2:length(behaviorOcc)
   b(i/2)=behaviorOcc(i); 
   tf(i/2)=nanmean(totalF((i-1):i));
   tfd(i/2)=nanmean(totalFDensity((i-1):i));
   st(i/2)=nanmean(spotTotal((i-1):i));
   sd(i/2)=nanmean(spotDensity((i-1):i));
   sm(i/2)=nanmean(spotMedianSize((i-1):i));
end

%%

figure %1
plot(tf,b,'.','Color',ocolor,'LineWidth',3,'MarkerSize',15)
[r p]=corrcoef(tf,b);
text(min(tf),min(b)+.1,['r = ' num2str(r(1,2),'%02.2f')],'FontSize',15)
text(min(tf),min(b),['p = ' num2str(p(1,2),'%02.2f')],'FontSize',15)
xlabel('Brp-Short fluorescence')
ylabel('preference score')
set(gca,'FontSize',15)
axis square

figure %2
plot(tfd,b,'.','Color',ocolor,'LineWidth',3,'MarkerSize',15)
[r p]=corrcoef(tfd,b);
text(min(tfd),min(b)+.1,['r = ' num2str(r(1,2),'%02.2f')],'FontSize',15)
text(min(tfd),min(b),['p = ' num2str(p(1,2),'%02.2f')],'FontSize',15)
xlabel('Brp-Short fluorescence density')
ylabel('preference score')
set(gca,'FontSize',15)
axis square

figure %3
plot(st,b,'.','Color',ocolor,'LineWidth',3,'MarkerSize',15)
[r p]=corrcoef(st,b);
text(min(st),min(b)+.1,['r = ' num2str(r(1,2),'%02.2f')],'FontSize',15)
text(min(st),min(b),['p = ' num2str(p(1,2),'%02.2f')],'FontSize',15)
xlabel('# of puncta')
ylabel('preference score')
set(gca,'FontSize',15)
axis square

% FIG 3j
figure %4
plot(sd,b,'.','Color',ocolor,'LineWidth',3,'MarkerSize',15)
[r p]=corrcoef(sd,b);
text(min(sd),min(b)+.1,['r = ' num2str(r(1,2),'%02.2f')],'FontSize',15)
text(min(sd),min(b),['p = ' num2str(p(1,2),'%02.2f')],'FontSize',15)
xlabel('puncta density')
ylabel('preference score')
set(gca,'FontSize',15)
axis([0.0025 0.0225 -0.6 0.2])
axis square

% FIG 3k
figure %5
plot(sm,b,'.','Color',ocolor,'LineWidth',3,'MarkerSize',15)
[r p]=corrcoef(sm,b);
text(min(sm),min(b)+.1,['r = ' num2str(r(1,2),'%02.2f')],'FontSize',15)
text(min(sm),min(b),['p = ' num2str(p(1,2),'%02.2f')],'FontSize',15)
xlabel('median puncta volume (um^3)')
ylabel('preference score')
set(gca,'FontSize',15)
axis([1 4 -0.6 0.2])
axis square

% FIG 3i
figure %6
imagesc(1-binarizedSpots{8}(:,:,65))
colormap('gray')