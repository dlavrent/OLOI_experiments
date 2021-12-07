clear all
close all

foldernames{1} = '191009_oct_vs_air_persistence';
foldernames{2} = '191010_oct_vs_air_persistence';
prefix = 'oct_air';


nbatches{1} = 3;
nbatches{2} = 4;
z=[];
zv=[];
zp=[];
t=[];
tv=[];
tp=[];
for i =1:2
    cd(foldernames{i})
    load([prefix '_twentyfourHourPositions.mat'])
    [sorted sortedInds]=sort(randomizedPositions,2);
    
    velocityThresh = 0.2;
    zerohour=zeros(nbatches{i},15);
    zerohourvelocity=zeros(nbatches{i},15);
    twentyFour=zeros(nbatches{i},15);
    twentyFourvelocity=zeros(nbatches{i},15);
    
    zcounter=1;
    tcounter=1;
    fnames = dir(pwd);
    for i = 1:length(fnames)
        if strfind(fnames(i).name,prefix)
            if strfind(fnames(i).name,'processed')
                load(fnames(i).name)
                if strfind(fnames(i).name, '00h')
                    temp=0;
                    for j=1:length(flyTracks.tunnelActive)
                        if flyTracks.tunnelActive(j)==1
                            temp=temp+1;
                            zerohour(zcounter,j)=flyTracks.occupancy(temp);
                            zerohourvelocity(zcounter,j)=prctile(flyTracks.velocity(:,temp),50);
                            zerohourp(zcounter,j)=flyTracks.preOdorOccupancy(temp);
                        else
                            zerohour(zcounter,j)=NaN;
                            zerohourvelocity(zcounter,j)=NaN;
                            zerohourp(zcounter,j)=NaN;
                        end
                    end
                    zcounter=zcounter+1;
                elseif strfind(fnames(i).name, '24h')
                    temp=0;
                    for j=1:length(flyTracks.tunnelActive)
                        if flyTracks.tunnelActive(j)==1
                            temp=temp+1;
                            twentyFour(tcounter,j)=flyTracks.occupancy(temp);
                            twentyFourvelocity(tcounter,j)=prctile(flyTracks.velocity(:,temp),50);
                            twentyFourp(tcounter,j)=flyTracks.preOdorOccupancy(temp);
                        else
                            twentyFour(tcounter,j)=NaN;
                            twentyFourvelocity(tcounter,j)=NaN;
                            twentyFourp(tcounter,j)=NaN;
                        end
                    end
                    
                    twentyFour(tcounter,:)=twentyFour(tcounter,sortedInds(tcounter,:));
                    twentyFourvelocity(tcounter,:)=twentyFourvelocity(tcounter,sortedInds(tcounter,:));
                    twentyFourp(tcounter,:)=twentyFourp(tcounter,sortedInds(tcounter,:));
                    tcounter=tcounter+1;
                end
            end
        end
    end
    cd ..
    z = [z; zerohour];
    zv = [zv; zerohourvelocity];
    zp = [zp; zerohourp];
    t = [t; twentyFour];
    tv = [tv; twentyFourvelocity];
    tp = [tp; twentyFourp];
end

%%
z = z(:);
zv = zv(:);
zp = zp(:);
t = t(:);
tv = tv(:);
tp = tp(:);

todelete = [];
todelete = [todelete transpose(find(isnan(t)))];
todelete = [todelete transpose(find(isnan(z)))];
todelete = [todelete transpose(find(tv<velocityThresh))];
todelete = [todelete transpose(find(zv<velocityThresh))];

z(todelete)=[];
zp(todelete)=[];
tp(todelete)=[];
zv(todelete)=[];
tv(todelete)=[];
t(todelete)=[];


[r p] = corrcoef(z-zp,t-tp);
figure;
plot(z-zp,t-tp,'k.','LineWidth',2,'MarkerSize',15)
text(0, 0, ['r = ' num2str(r(1,2),'%02.3f')],'FontSize',15)
text(0.1, 0, ['p = ' num2str(p(1,2),'%02.3f')],'FontSize',15)
xlabel('odor preference (t = 0 h)')
ylabel('odor preference (t = 24 h)')
axis square
axis([-1 .3 -1 .3])
set(gca,'FontSize',15)
set(gca,'xtick','')
set(gca,'ytick','')

% 
% [rp pp] = corrcoef(zp,tp);
% figure;
% plot(zp,tp,'ko','LineWidth',2)
% xlabel('pre-odor occupancy (t = 0 h)')
% ylabel('pre-odor occupancy (t = 24 h)')
% text(0.5, 0.5, ['r = ' num2str(rp(1,2))],'FontSize',15)
% box off
% set(gca,'FontSize',15)
%% bootstrap
iters = 10000;
rboot=zeros(1,iters);
for i=1:iters
    bsample = zeros(1,length(z));
    for j = 1:length(z)
        bsample(j)=round(rand(1)*(length(z)-1))+1;
    end
    zsample=z(bsample);
    tsample=t(bsample);
    r = corrcoef(zsample,tsample);
    rboot(i)=r(1,2);
end

figure;
histogram(rboot)