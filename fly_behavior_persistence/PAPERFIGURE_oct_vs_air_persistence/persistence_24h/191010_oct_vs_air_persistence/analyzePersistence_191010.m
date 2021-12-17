clear all
close all

homedir='/Users/mattchurgin/Dropbox (Harvard University)/fly_behavior/191010_oct_vs_air_persistence';
cd(homedir)

prefix = 'oct_air';
nbatches = 4;
load([prefix '_twentyfourHourPositions.mat'])
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
            elseif strfind(fnames(i).name, '24h')
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

[r p] = corrcoef(z,t);
figure;
plot(z,t,'ko','LineWidth',2)
text(0.5, 0.5, ['r = ' num2str(r(1,2))],'FontSize',15)
xlabel('occupancy (t = 0 h)')
ylabel('occupancy (t = 24 h)')
box off
set(gca,'FontSize',15)

[r p] = corrcoef(z-zp,t-tp);
figure;
plot(z-zp,t-tp,'ko','LineWidth',2)
text(0, 0, ['r = ' num2str(r(1,2))],'FontSize',15)
xlabel('occupancy (t = 0 h)')
ylabel('occupancy (t = 24 h)')
box off
set(gca,'FontSize',15)

% 
% [rp pp] = corrcoef(zp,tp);
% figure;
% plot(zp,tp,'ko','LineWidth',2)
% xlabel('pre-odor occupancy (t = 0 h)')
% ylabel('pre-odor occupancy (t = 24 h)')
% text(0.5, 0.5, ['r = ' num2str(rp(1,2))],'FontSize',15)
% box off
% set(gca,'FontSize',15)

%% plot effect of speed and preference
figure
[rp pp] = corrcoef(zv,z);
[tp asdf] = corrcoef(tv,t);
plot(zv,z,'ko','LineWidth',3)
hold on
plot(tv,t,'r*','LineWidth',3)
legend('day 1','day 2')
text(0.5, 0.85, ['day 1 r = ' num2str(rp(1,2))],'FontSize',15)
text(0.5, 0.5, ['day 2 r = ' num2str(tp(1,2))],'FontSize',15)
xlabel('median velocity')
ylabel('preference')
legend boxoff
box off
set(gca,'FontSize',15)
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