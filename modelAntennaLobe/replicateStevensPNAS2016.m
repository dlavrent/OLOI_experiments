clear all
close all

fdir='C:\Users\mac0456\Dropbox\flyimaging\analysis\modelAntennaLobe';
fname='hallemAndCarlson2004_tableS1data.xlsx';
[NUM,TXT,RAW]=xlsread([fdir '\' fname]);
ORNr=NUM(1:110,:)'+NUM(111,:)';  % add background firing rate for each OR

% sort responses for each odor
sortedResponses=sort(ORNr);


rates=[-50:300];
cdfresponses=zeros(length(rates),110);
for j=1:110
   for i=1:length(rates)
       cdfresponses(i,j)=sum(ORNr(:,j)<rates(i))/length(ORNr(:,j));
   end
end

figure
hold on
for i=1:110
    plot(rates,cdfresponses(:,i),'k')
end

%% scale rates to all have mean of 100


rates=[-50:300];
cdfresponses=zeros(length(rates),110);
for j=1:110
   for i=1:length(rates)
       temp=ORNr(:,j)*100/mean(ORNr(:,j));
       cdfresponses(i,j)=sum(temp<rates(i))/length(temp);
   end
end

figure
hold on
for i=1:110
    plot(rates,cdfresponses(:,i),'k')
end

%% sample from exponential distribution with mean of 100 hz
lambda=1/100; % yields exponential distribution with mean of 100 hz

 % % testing individual sample from exponential distribution
% for i=1:10000
%     u=rand(1);
% exponentialsample(i)=(-1/lambda)*log(1-u);
% end

% generate 110 collections of norn random samples
norns=24;
sampledResponses=zeros(110,norns);
for j=1:110
    for i=1:norns
        u=rand(1);
        sampledResponses(j,i)=(-1/lambda)*log(1-u);
    end
end

rates=[-50:300];
cdfsampledresponses=zeros(length(rates),110);
for j=1:110
   for i=1:length(rates)
       temp=sampledResponses(j,:);
       cdfsampledresponses(i,j)=sum(temp<rates(i))/length(temp);
   end
end

figure
hold on
for i=1:110
    plot(rates,cdfsampledresponses(:,i),'k')
end