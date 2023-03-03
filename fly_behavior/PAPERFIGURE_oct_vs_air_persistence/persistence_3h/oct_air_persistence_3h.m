clear all
close all

load analysis_dir_path
data_folder = fullfile(analysis_dir_path, 'fly_behavior\PAPERFIGURE_oct_vs_air_persistence\persistence_3h');
cd(data_folder)

prefix = 'oct_air';
nbatches = 4;
load([prefix '_threeHourPositions.mat'])
[sorted sortedInds]=sort(randomizedPositions,2);

velocityThresh = 0.2;
zerohour=zeros(nbatches,15);
zerohourvelocity=zeros(nbatches,15);
threehour=zeros(nbatches,15);
threehourvelocity=zeros(nbatches,15);

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
            elseif strfind(fnames(i).name, '03h')
                temp=0;
                for j=1:length(flyTracks.tunnelActive)
                    if flyTracks.tunnelActive(j)==1
                        temp=temp+1;
                        threehour(tcounter,j)=flyTracks.occupancy(temp);
                        threehourvelocity(tcounter,j)=prctile(flyTracks.velocity(:,temp),50);
                        threehourp(tcounter,j)=flyTracks.preOdorOccupancy(temp);
                    else
                        threehour(tcounter,j)=NaN;
                        threehourvelocity(tcounter,j)=NaN;
                        threehourp(tcounter,j)=NaN;
                    end
                end
                
                threehour(tcounter,:)=threehour(tcounter,sortedInds(tcounter,:));
                threehourvelocity(tcounter,:)=threehourvelocity(tcounter,sortedInds(tcounter,:));
                threehourp(tcounter,:)=threehourp(tcounter,sortedInds(tcounter,:));
                tcounter=tcounter+1;
            end
        end
    end
end

z = zerohour(:);
zv = zerohourvelocity(:);
zp = zerohourp(:);
t = threehour(:);
tv = threehourvelocity(:);
tp = threehourp(:);

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
% SUP FIG OCT-AIR persistence 3h
figure; %1
plot(z-zp,t-tp,'k.','LineWidth',3)
text(0, 0, ['r = ' num2str(r(1,2))],'FontSize',15)
xlabel('occupancy (t = 0 h)')
ylabel('occupancy (t = 3 h)')
box on
axis square
set(gca,'xtick','')
set(gca,'ytick','')
axis([-.8 .5 -.8 .5])
set(gca,'FontSize',15)
% 
% [rp pp] = corrcoef(zp,tp);
% figure;
% plot(zp,tp,'ko','LineWidth',2)
% xlabel('pre-odor occupancy (t = 0 h)')
% ylabel('pre-odor occupancy (t = 3 h)')
% text(0.5, 0.5, ['r = ' num2str(rp(1,2))],'FontSize',15)
% box off
% set(gca,'FontSize',15)

