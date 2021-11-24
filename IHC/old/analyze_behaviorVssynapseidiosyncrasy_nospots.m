clear all
close all

homedir=pwd;
ornOrPn=1; % 1 for pn, 0 for orn
densityOrTotal=1; % 0 for density, 1 for total fluorescence
saveData=0; 
if ornOrPn
    processdir='/Users/mattchurgin/Dropbox/flyimaging/analysis/IHC/Gh146Dalpha7/gh146OX10Dalpha7OXZ';
    processdir='C:\Users\mac0456\Dropbox\flyimaging\analysis\IHC\Gh146Dalpha7\gh146OX10Dalpha7OXZ';
    mycolor=[0 0.7 0];
else
    processdir='C:\Users\mac0456\Dropbox\flyimaging\analysis\IHC\OrcoBrpshort\orcoOX5brpshortOXZ';
    mycolor=[0.9 0.0 0.7];
end
currfolders=dir(processdir);
currfolders=currfolders(3:end);

punctadiffL=zeros(1,length(currfolders));
punctadiffR=zeros(1,length(currfolders));
punctadiffM=zeros(1,length(currfolders));

behaviorpre=zeros(1,length(currfolders));
behaviorraw=zeros(1,length(currfolders));

cd(processdir)
for i=1:length(currfolders)
    cd(currfolders(i).name)
    currfiles=dir(pwd);
    currfiles=currfiles(3:end);
    
    load(currfiles(1).name)
    
    if densityOrTotal==0
        punctadiffL(i)=100*(dm2l_volnormed-dc2l_volnormed)/(dm2l_volnormed+dc2l_volnormed);
        punctadiffR(i)=100*(dm2r_volnormed-dc2r_volnormed)/(dm2r_volnormed+dc2r_volnormed);
        punctadiffM(i)=nanmean([punctadiffL(i) punctadiffR(i)]);
        % save puncta density for each glomerulus
        
        dm2lp(i)=dm2l_volnormed;
        dm2rp(i)=dm2r_volnormed;
        dc2lp(i)=dc2l_volnormed;
        dc2rp(i)=dc2r_volnormed;
        dm2lv(i)=sum(sum(sum(dm2Lmask)));
        dm2rv(i)=sum(sum(sum(dm2Rmask)));
        dc2lv(i)=sum(sum(sum(dc2Lmask)));
        dc2rv(i)=sum(sum(sum(dc2Rmask)));
    else
        % calculate total fluorescence rather than density
        dm2lp(i)=dm2l_volnormed*sum(sum(sum(dm2Lmask)));
        dm2rp(i)=dm2r_volnormed*sum(sum(sum(dm2Rmask)));
        dc2lp(i)=dc2l_volnormed*sum(sum(sum(dc2Lmask)));
        dc2rp(i)=dc2r_volnormed*sum(sum(sum(dc2Rmask)));
        punctadiffL(i)=100*(dm2lp(i)-dc2lp(i))/(dm2lp(i)+dc2lp(i));
        punctadiffR(i)=100*(dm2rp(i)- dc2rp(i))/(dm2rp(i)+ dc2rp(i));
        punctadiffM(i)=nanmean([punctadiffL(i) punctadiffR(i)]);
        
        dm2lv(i)=sum(sum(sum(dm2Lmask)));
        dm2rv(i)=sum(sum(sum(dm2Rmask)));
        dc2lv(i)=sum(sum(sum(dc2Lmask)));
        dc2rv(i)=sum(sum(sum(dc2Rmask)));
    end
    
    load(currfiles(2).name)
    behaviorpre(i)=preocc;
    behaviorraw(i)=occ;
    
    disp(['loaded fly ' num2str(i)])
    cd ..
end

disp('loaded all flies')


% find missing lobe data
lmissing=isnan(punctadiffL);
rmissing=isnan(punctadiffR);

missingALobe=(lmissing+rmissing)>0;

% remove flies for which we don't have data from both lobes

dm2lp(missingALobe)=[];
dm2rp(missingALobe)=[];
dc2lp(missingALobe)=[];
dc2rp(missingALobe)=[];
dm2lv(missingALobe)=[];
dm2rv(missingALobe)=[];
dc2lv(missingALobe)=[];
dc2rv(missingALobe)=[];

pL=punctadiffL; pL(missingALobe)=[];
pR=punctadiffR; pR(missingALobe)=[];
pM=punctadiffM; pM(missingALobe)=[];

bRaw=behaviorraw; bRaw(missingALobe)=[];
bPre=behaviorpre; bPre(missingALobe)=[];
b=bRaw-bPre;

numFly=length(b);
cd(homedir)
if saveData==1
    
    
    if densityOrTotal==0
        if ornOrPn==1
            save('pn_data_density.mat','dm2lp','dm2rp','dc2lp','dc2rp','dm2lv','dm2rv','dc2lv','dc2rv','pL','pR','pM','bRaw', ...
                'bPre','b','numFly','densityOrTotal','mycolor')
        else
            save('orn_data_density.mat','dm2lp','dm2rp','dc2lp','dc2rp','dm2lv','dm2rv','dc2lv','dc2rv','pL','pR','pM','bRaw', ...
                'bPre','b','numFly','densityOrTotal','mycolor')
        end
    else
        if ornOrPn==1
            save('pn_data.mat','dm2lp','dm2rp','dc2lp','dc2rp','dm2lv','dm2rv','dc2lv','dc2rv','pL','pR','pM','bRaw', ...
                'bPre','b','numFly','densityOrTotal','mycolor')
        else
            save('orn_data.mat','dm2lp','dm2rp','dc2lp','dc2rp','dm2lv','dm2rv','dc2lv','dc2rv','pL','pR','pM','bRaw', ...
                'bPre','b','numFly','densityOrTotal','mycolor')
        end
    end
end


%% load previously saved data
ornOrPn=1; % 0 for ORN, 1 for PN
densityOrTotal=0;
if densityOrTotal==0
    if ornOrPn==0
        load('orn_data_density.mat')
    else
        load('pn_data_density.mat')
    end
else
    if ornOrPn==0
        load('orn_data.mat')
    else
        load('pn_data.mat')
    end
end


%%  plot left vs right side puncta density and glomerulus volume
figure
plot(dc2lp,dc2rp,'o','Color',[0.05 0.05 0.05],'LineWidth',3)
hold on
plot(dm2lp,dm2rp,'*','Color',[0.65 0.65 0.65],'LineWidth',3)
if ~densityOrTotal
    xlabel('puncta density (left)')
    ylabel('puncta density (right)')
else
    xlabel('total fluorescence (left)')
    ylabel('total fluorescence (right)')
end
legend('DC2','DM2')
legend boxoff
%text(dc2lp(j),dc2rp(j),[num2str((j))],'FontSize',15)
set(gca,'FontSize',15)
box off


figure
plot(dc2lv,dc2rv,'o','Color',[0.05 0.05 0.05],'LineWidth',3)
hold on
plot(dm2lv,dm2rv,'*','Color',[0.65 0.65 0.65],'LineWidth',3)
    xlabel('voxel volume (left)')
    ylabel('voxel volume (right)')
legend('DC2','DM2')
legend boxoff
%text(dc2lp(j),dc2rp(j),[num2str((j))],'FontSize',15)
set(gca,'FontSize',15)
box off

%% plot dm2 vs dc2 volume

figure
plot([dm2lv],[dc2lv],'o','Color',[0.05 0.05 0.05],'LineWidth',3)
hold on
plot(dm2rv,dc2rv,'*','Color',[0.65 0.65 0.65],'LineWidth',3)
    xlabel('voxel volume (DM2)')
    ylabel('voxel volume (DC2)')
legend('left','right')
legend boxoff
%text(dc2lp(j),dc2rp(j),[num2str((j))],'FontSize',15)
set(gca,'FontSize',15)
box off

%% plot total fluorescence vs total volume
figure
plot([dc2lv dc2rv],[dc2lp dc2rp],'o','Color',[0.05 0.05 0.05],'LineWidth',3)
hold on
plot([dm2lv dm2rv],[dm2lp dm2rp],'*','Color',[0.65 0.65 0.65],'LineWidth',3)
    xlabel('voxel volume')
    ylabel('total fluorescence')
legend('DC2','DM2')
legend boxoff
%text(dc2lp(j),dc2rp(j),[num2str((j))],'FontSize',15)
set(gca,'FontSize',15)
box off

%% plot dm2 v dc2 total fluorescence 
figure
plot(dm2lp,dc2rp,'o','Color',[0.05 0.05 0.05],'LineWidth',3)
hold on
plot(dm2rp,dc2rp,'*','Color',[0.65 0.65 0.65],'LineWidth',3)
if ~densityOrTotal
    xlabel('puncta density (DM2)')
    ylabel('puncta density (DC2)')
else
    xlabel('total fluorescence (DM2)')
    ylabel('total fluorescence (DC2)')
end
legend('left','right')
legend boxoff
%text(dc2lp(j),dc2rp(j),[num2str((j))],'FontSize',15)
set(gca,'FontSize',15)
box off


%%  plot left vs right side puncta density difference

figure
hold on
[rp pp]=corrcoef(pL,pR);

plot(pL,pR,'o','Color',mycolor,'LineWidth',3)

if ~densityOrTotal
    xlabel('DM2-DC2 puncta density (left)')
    ylabel('DM2-DC2 puncta density (right)')
else
    xlabel('DM2-DC2 total fluorescence (left)')
    ylabel('DM2-DC2 total fluorescence (right)')
end

%text(pL(j),pR(j),[num2str((j))],'FontSize',15)
set(gca,'FontSize',15)
box off


%% plot puncta density against behavior

[rp pp]=corrcoef(pM,b);
figure;
hold on
for j=1:numFly
    plot(pM(j),b(j),'o','Color',mycolor,'LineWidth',3)
    
    text(pM(j),b(j),[num2str((j))],'FontSize',15)
end
%text(0,0,['r = ' num2str(rp(1,2))],'FontSize',15)
if ~densityOrTotal
    xlabel('DM2-DC2 labelling density (% difference)')
else
    xlabel('DM2-DC2 total fluorescence (% difference)')
end
ylabel('measured preference')
set(gca,'FontSize',15)
box off



figure;
hold on
for j=1:numFly
    plot((dm2lp(j)+dm2rp(j))/2,b(j),'o','Color',mycolor,'LineWidth',3)
    
    %text((dm2lp(j)+dm2rp(j))/2,b(j),[num2str((j))],'FontSize',15)
end
%text(0,0,['r = ' num2str(rp(1,2))],'FontSize',15)
if ~densityOrTotal
    xlabel('DM2 labelling density')
else
    xlabel('DM2 total fluorescence')
end
ylabel('measured preference')
set(gca,'FontSize',15)
box off


figure;
hold on
for j=1:numFly
    plot((dc2lp(j)+dc2rp(j))/2,b(j),'o','Color',mycolor,'LineWidth',3)
    
   % text((dc2lp(j)+dc2rp(j))/2,b(j),[num2str((j))],'FontSize',15)
end
%text(0,0,['r = ' num2str(rp(1,2))],'FontSize',15)
if ~densityOrTotal
    xlabel('DC2 labelling density')
else
    xlabel('DC2 total fluorescence')
end
ylabel('measured preference')
set(gca,'FontSize',15)
box off



% make puncta diff using standardized data for dm2 and dc2
dm2ls=(dm2lp-mean(dm2lp))/std(dm2lp);
dm2rs=(dm2rp-mean(dm2rp))/std(dm2rp);
dc2ls=(dc2lp-mean(dc2lp))/std(dc2lp);
dc2rs=(dc2rp-mean(dc2rp))/std(dc2rp);
pLs=(dm2ls-dc2ls);
pRs=(dm2rs-dc2rs);
pMs=(pLs+pRs)/2;

[rp pp]=corrcoef(pMs,b);
figure;
hold on
for j=1:numFly
    plot(pMs(j),b(j),'o','Color',mycolor,'LineWidth',3)
    
    %text(pM(j),b(j),[num2str((j))],'FontSize',15)
end
%text(0,0,['r = ' num2str(rp(1,2))],'FontSize',15)
if ~densityOrTotal
    xlabel('DM2-DC2 labelling density')
else
    xlabel('DM2-DC2 total fluorescence')
end
ylabel('measured preference')
set(gca,'FontSize',15)
box off

%%   make correlation matrix

allData=[dm2lv; dm2rv; dm2lp; dm2rp; dc2lv; dc2rv; dc2lp; dc2rp; pL; pR; pM; bPre; bRaw; b];
varNames=[{'DM2Lv'},{'DM2Rv'},{'DM2Lf'},{'DM2Rf'},{'DC2Lv'},{'DC2Rv'},{'DC2Lf'},{'DC2Rf'},{'fdiffL'},{'fdiffR'},{'fdiffAv'},{'preOdorB'},{'odorB'},{'B'}];
dataCov=allData*allData';

[dataCorr pCorr]=corrcoef(allData');


figure
imagesc(dataCorr)
set(gca,'XTick',1:size(dataCorr,1))
set(gca,'XTickLabel',varNames)
xtickangle(30)
set(gca,'YTick',1:size(dataCorr,2))
set(gca,'YTickLabel',varNames)
ytickangle(0)
colorbar
set(gca,'FontSize',15)

figure
imagesc(abs(dataCorr))
set(gca,'XTick',1:size(dataCorr,1))
set(gca,'XTickLabel',varNames)
xtickangle(30)
set(gca,'YTick',1:size(dataCorr,2))
set(gca,'YTickLabel',varNames)
ytickangle(0)
colorbar
set(gca,'FontSize',15)



figure
imagesc(pCorr,[0 0.2])
set(gca,'XTick',1:size(dataCorr,1))
set(gca,'XTickLabel',varNames)
xtickangle(30)
set(gca,'YTick',1:size(dataCorr,2))
set(gca,'YTickLabel',varNames)
ytickangle(0)
colorbar
set(gca,'FontSize',15)
%% use bootstrap resampling to estimate correlation and p-value
iters=25000;
rb=zeros(1,iters-1);
pb=zeros(1,iters-1);
for i = 1:iters
    for j=1:numFly
        temp=rand(1);
        temp2=1+round(temp*(numFly-1));
        bootstrapb(j)=b(temp2);
        bootstrapp(j)=pM(temp2);
    end
    
    [r p]=corrcoef(bootstrapp,bootstrapb);
    rb(i)=r(1,2);
    pb(i)=p(1,2);
    
    clear bootstrapp bootstrapb
    if mod(i,5000)==0
        disp(['iter ' num2str(i) ' of ' num2str(iters)])
    end
end

figure;
histogram(rb,50,'Normalization','probability')
xlabel('bootstrapped r')
ylabel('frequency')
set(gca,'FontSize',15)
box off
% 
% 
% figure;
% histogram(pb,50)
% xlabel('bootstrapped p')
% ylabel('counts')
% set(gca,'FontSize',15)
% box off

hold on
text(0,0,['bootstrapped p = ' num2str(round(1000*(1-sum(rb<0)/length(rb)))/1000)],'FontSize',15)

%% use bootstrap resampling to estimate correlation and p-value for just dm2 or dc2 behavior prediction
iters=25000;
rdm2=zeros(1,iters-1);
pdm2=zeros(1,iters-1);
rdc2=zeros(1,iters-1);
pdc2=zeros(1,iters-1);
for i = 1:iters
    for j=1:numFly
        temp=rand(1);
        temp2=1+round(temp*(numFly-1));
        bootstrapb(j)=b(temp2);
        bootstrapdc2(j)=(dc2lp(temp2)+dc2rp(temp2))/2;
        bootstrapdm2(j)=(dm2lp(temp2)+dm2rp(temp2))/2;
    end
    
    [r p]=corrcoef(bootstrapdm2,bootstrapb);
    rdm2(i)=r(1,2);
    pdm2(i)=p(1,2);
    [r p]=corrcoef(bootstrapdc2,bootstrapb);
    rdc2(i)=r(1,2);
    pdc2(i)=p(1,2);
    
    if mod(i,5000)==0
        disp(['iter ' num2str(i) ' of ' num2str(iters)])
    end
end

figure;
histogram(rdm2,50,'Normalization','probability')
xlabel('bootstrapped r (DM2)')
ylabel('frequency')
set(gca,'FontSize',15)
box off
hold on
text(0,0,['bootstrapped p = ' num2str(round(1000*(1-sum(rdm2<0)/length(rdm2)))/1000)],'FontSize',15)

figure;
histogram(rdc2,50,'Normalization','probability')
xlabel('bootstrapped r (DC2)')
ylabel('frequency')
set(gca,'FontSize',15)
box off
hold on
text(0,0,['bootstrapped p = ' num2str(round(1000*(1-sum(rdc2<0)/length(rdc2)))/1000)],'FontSize',15)
