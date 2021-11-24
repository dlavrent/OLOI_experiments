clear all
close all

samplesizes=10:10:100;
r = 0.4;
N = 10000;
alpha = 0.05;
for n = samplesizes
    out = sampleCorrelated(r,n,N);
    pow(n)=sum(out.sampleStats(:,2)<alpha)/N;
    n

end

%
r = 0.3;
for n = samplesizes
    out = sampleCorrelated(r,n,N);
    pow2(n)=sum(out.sampleStats(:,2)<alpha)/N;
    n

end

%%
figure
plot(samplesizes,pow(samplesizes),'k','LineWidth',3)
hold on
plot(samplesizes,pow2(samplesizes),'m','LineWidth',3)
plot(0:100,0.8*ones(length(0:100)),'r--','LineWidth',2)
legend('true r = 0.4','true r = 0.3')
legend boxoff
xlabel('sample size')
ylabel('power')
box off
set(gca,'FontSize',15)