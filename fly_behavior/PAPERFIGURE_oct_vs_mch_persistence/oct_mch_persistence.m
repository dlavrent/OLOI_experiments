clear all 
close all
filesep='/';
foldernames{1}='190123_isokh11_f_0_3_24_hour_pref';
foldernames{2}='190124_isokh11_f_0_3_24_hour_pref';
foldernames{3}='190125_isokh11_f_0_3_24_hour_pref';

hour0{1}=load([foldernames{1} filesep 'flies_01_to_10_00_hour-processed.mat']);
hour3{1}=load([foldernames{1} filesep 'flies_01_to_10_03_hour-processed.mat']);
hour24{1}=load([foldernames{1} filesep 'flies_01_to_10_24_hour-processed.mat']);

hour0{2}=load([foldernames{2} filesep 'flies_01_to_15_00_hour-processed.mat']);
hour3{2}=load([foldernames{2} filesep 'flies_01_to_15_03_hour-processed.mat']);
hour24{2}=load([foldernames{2} filesep 'flies_01_to_15_24_hour-processed.mat']);

hour0{3}=load([foldernames{3} filesep 'flies_01_to_15_00_hour-processed.mat']);
hour3{3}=load([foldernames{3} filesep 'flies_01_to_15_03_hour-processed.mat']);
hour24{3}=load([foldernames{3} filesep 'flies_01_to_15_24_hour-processed.mat']);

hour0{4}=load([foldernames{3} filesep 'flies_16_to_27_00_hour-processed.mat']);
hour3{4}=load([foldernames{3} filesep 'flies_16_to_27_03_hour-processed.mat']);
hour24{4}=load([foldernames{3} filesep 'flies_16_to_27_24_hour-processed.mat']);

numExp=length(hour0);


%% experiments 3 hour persistence
hour0all=[];
hour3all=[];
hour0allPreodor=[];
hour3allPreodor=[];
for i = 1:numExp
    hour0all=[hour0all hour0{i}.flyTracks.occupancy-hour0{i}.flyTracks.preOdorOccupancy];
    hour3all=[hour3all hour3{i}.flyTracks.occupancy-hour3{i}.flyTracks.preOdorOccupancy];
    hour0allPreodor=[hour0allPreodor hour0{i}.flyTracks.preOdorOccupancy];
    hour3allPreodor=[hour3allPreodor hour3{i}.flyTracks.preOdorOccupancy];
end

% remove NaN
missing1=find(isnan(hour0all));
missing2=find(isnan(hour3all));

hour0all([missing1 missing2])=[];
hour3all([missing1 missing2])=[];

hour0allPreodor([missing1 missing2])=[];
hour3allPreodor([missing1 missing2])=[];

% SUP FIG 1d
figure %1
plot(hour0all,hour3all,'k.','LineWidth',3)
[r p1]=corrcoef(hour0all,hour3all);
text(0.1,0.5,['r = ' num2str(r(1,2),'%02.3f')],'FontSize',15)
text(0.4,0.5,['p = ' num2str(p1(1,2),'%02.3f')],'FontSize',15)
xlabel('odor preference (t = 0 h)')
ylabel('odor preference (t = 3 h)')
box on
set(gca,'xtick','')
set(gca,'ytick','')
axis square
axis([-.8 0.5 -.8 0.5])
set(gca,'FontSize',15)

% figure
% plot(hour0allPreodor,hour3allPreodor,'k.','LineWidth',3)
% r=corrcoef(hour0allPreodor,hour3allPreodor);
% text(0.1,0.5,['r = ' num2str(r(1,2))],'FontSize',15)
% xlabel('pre-odor preference (t = 0 h)')
% ylabel('pre-odor preference (t = 3 h)')
% box off
% set(gca,'FontSize',15)
% 


%% 24 hours persistence
hour0all=[];
hour24all=[];
hour0allPreodor=[];
hour24allPreodor=[];
for i = 1:2
    hour0all=[hour0all hour0{i}.flyTracks.occupancy-hour0{i}.flyTracks.preOdorOccupancy];
    hour24all=[hour24all hour24{i}.flyTracks.occupancy-hour24{i}.flyTracks.preOdorOccupancy];
    hour0allPreodor=[hour0allPreodor hour0{i}.flyTracks.preOdorOccupancy];
    hour24allPreodor=[hour24allPreodor hour24{i}.flyTracks.preOdorOccupancy];
end
% manually add 3rd day's experiments because the first and last fly from
% that day was lost at 24 hour time point
i=3;
hour0all=[hour0all hour0{i}.flyTracks.occupancy];
hour24all=[hour24all NaN hour24{i}.flyTracks.occupancy];
hour0allPreodor=[hour0allPreodor hour0{i}.flyTracks.preOdorOccupancy];
hour24allPreodor=[hour24allPreodor NaN hour24{i}.flyTracks.preOdorOccupancy];

i=4;
hour0all=[hour0all hour0{i}.flyTracks.occupancy];
hour24all=[hour24all hour24{i}.flyTracks.occupancy NaN];
hour0allPreodor=[hour0allPreodor hour0{i}.flyTracks.preOdorOccupancy];
hour24allPreodor=[hour24allPreodor hour24{i}.flyTracks.preOdorOccupancy NaN];

% remove NaN
missing1=find(isnan(hour0all));
missing2=find(isnan(hour24all));

hour0all([missing1 missing2])=[];
hour24all([missing1 missing2])=[];

hour0allPreodor([missing1 missing2])=[];
hour24allPreodor([missing1 missing2])=[];

% SUP FIG 1e
figure %2
plot(hour0all,hour24all,'k.','LineWidth',2,'MarkerSize',15)
[r p1]=corrcoef(hour0all,hour24all);
text(0.1,0.5,['r = ' num2str(r(1,2),'%02.3f')],'FontSize',15)
text(0.4,0.5,['p = ' num2str(p1(1,2),'%02.3f')],'FontSize',15)
xlabel('odor preference (t = 0 h)')
ylabel('odor preference (t = 24 h)')
box on
set(gca,'xtick','')
set(gca,'ytick','')
axis square
axis([-1 1 -1 1])
set(gca,'FontSize',15)
% 
% figure
% plot(hour0allPreodor,hour24allPreodor,'*','LineWidth',3)
% [r p2]=corrcoef(hour0allPreodor,hour24allPreodor);
% text(0.1,0.5,['r = ' num2str(r(1,2))],'FontSize',15)
% xlabel('pre-odor preference (t = 0 h)')
% ylabel('pre-odor preference (t = 24 h)')
% box off
% set(gca,'FontSize',15)
