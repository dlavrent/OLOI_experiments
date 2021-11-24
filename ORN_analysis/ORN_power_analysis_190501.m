%% ORN power analysis
clear all
close all

nflies=35;
iters=10000;
trueCorr=sqrt(0.2);

for i=1:iters
   r=mvnrnd(zeros(nflies,2),[1 trueCorr; trueCorr 1]);
   rr=corrcoef(r);
   
   rp(i)=rr(1,2);
end

possiblecorrs=[-1:0.01:1];
for i=1:length(possiblecorrs)
   rcdf(i)=sum(rp<possiblecorrs(i))/length(rp); 
   
   if i>1
    rprob(i)=(sum(rp<possiblecorrs(i))-sum(rp<possiblecorrs(i-1)))/length(rp); 
   end
end

nullhpoint=find(rcdf>0.05);
nullhpoint=nullhpoint(1)-1;

figure;
plot(possiblecorrs,rcdf,'k','LineWidth',3)
title(['true r = ' num2str(trueCorr) ', n = ' num2str(nflies) ' flies'])
xlabel('measured r')
ylabel('CDF')
set(gca,'FontSize',15)
box off

figure;
plot(possiblecorrs,rcdf,'k','LineWidth',3)
hold on
plot(possiblecorrs(nullhpoint),rcdf(nullhpoint),'m*','LineWidth',3,'MarkerSize',10)
title(['true r = ' num2str(trueCorr) ', n = ' num2str(nflies) ' flies'])
xlabel('measured r')
ylabel('CDF')
set(gca,'FontSize',15)
box off


figure;
plot(possiblecorrs,rprob,'k','LineWidth',3)
title(['true r = ' num2str(trueCorr) ', n = ' num2str(nflies) ' flies'])
xlabel('measured r')
ylabel('probability')
set(gca,'FontSize',15)
box off
