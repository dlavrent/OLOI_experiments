% model of null odor preference distribution
clear all
close all


meanDecisions=30; % hard code # decisions
odorAprob=0.4;

nflies=1000;

for j=1:nflies
    flydecisions(j) = binornd(meanDecisions,odorAprob);
end

nbins=meanDecisions/2;
xbins=linspace(0,1,nbins);
figure;
[yf xf]=hist(flydecisions/meanDecisions,xbins);
bar(xf,yf/sum(yf))
xlabel('preference')
ylabel('probability')
box off
set(gca,'FontSize',20)

%% model number of decisions as variable

stdDecisions=10; % normal distribution std of # decisions

nflydecision=zeros(1,nflies);
flydecisions2=zeros(1,nflies);
for j=1:nflies
    nflydecision(j)=max(round(stdDecisions*randn(1)+meanDecisions),1);
    flydecisions2(j) = binornd(nflydecision(j),odorAprob);
end

figure;
[yf2 xf2]=hist(flydecisions2./nflydecision,xbins);
bar(xf2,yf2/sum(yf2))
xlabel('preference')
ylabel('probability')
box off
set(gca,'FontSize',20)


figure;
plot(xf,yf/sum(yf),'LineWidth',3)
hold on
plot(xf2,yf2/sum(yf2),'o-','LineWidth',3)
legend('fixed decisions','variable decisions')
xlabel('preference')
ylabel('probability')
box off
legend boxoff
set(gca,'FontSize',20)
