clear all
close all

load summary_200226
t=tunnelAssignment;
p=z;

load summary_200227
t=[t; tunnelAssignment];
p=[p; z];

tunnelAssignment=t;
z=p;

genotypes=[{'Or22a/+'},{'Or13a/+'},{'Or13a/LRP4-HA'},{'LRP4-HA/+'}, ...
    {'Or22a/LRP4-HA'}, {'TNT/+'}, {'Or22a/TNT'}, {'Or13a/TNT'}, {'IR-108629/+'}, ...
    {'Or13a/IR-108629'}, {'Or22a/IR-29900'}, {'Or13a/IR-29900'}, {'IR-29900/+'}, ...
    {'Or13a/Brp-Short'}, {'Or22a/Brp-Short'}, {'Brp-Short/+'}];

%% plot all together
figure
hold on
for i = 1:16
    plot([(i)*ones(1,length(find(tunnelAssignment==i)))'],[z(tunnelAssignment==i)],'o')
    plot([i],[mean(z(tunnelAssignment==i))],'k*-','LineWidth',3)
end
set(gca,'XTick',[1:length(genotypes)])
%set(gca,'XTickLabel',(genotypes))
set(gca,'XTickLabel',genotypes)
ylabel('preference')
set(gca,'FontSize',10)
axis([0 17 -1 1])


%%
fs=15;
ymin=-0.6;
ymax=0.5;

figure
subplot(1,2,1)
hold on
gens=[1 5 4];
c=0;
for i = gens
    c=c+1;
    plot([(c)*ones(1,length(find(tunnelAssignment==i)))'],[z(tunnelAssignment==i)],'o')
    plot([c],[mean(z(tunnelAssignment==i))],'k*-','LineWidth',3)
end
set(gca,'XTick',[1:length(gens)])
%set(gca,'XTickLabel',(genotypes))
set(gca,'XTickLabel',genotypes(gens))
ylabel('preference')
xtickangle(30)
set(gca,'FontSize',fs)
axis([0 length(gens)+1 ymin ymax])
subplot(1,2,2)
hold on
gens=[2 3 4];
c=0;
for i = gens
    c=c+1;
    plot([(c)*ones(1,length(find(tunnelAssignment==i)))'],[z(tunnelAssignment==i)],'o')
    plot([c],[mean(z(tunnelAssignment==i))],'k*-','LineWidth',3)
end
set(gca,'XTick',[1:length(gens)])
%set(gca,'XTickLabel',(genotypes))
set(gca,'XTickLabel',genotypes(gens))
ylabel('preference')
xtickangle(30)
set(gca,'FontSize',fs)
axis([0 length(gens)+1  ymin ymax])

figure
subplot(1,2,1)
hold on
gens=[1 7 6];
c=0;
for i = gens
    c=c+1;
    plot([(c)*ones(1,length(find(tunnelAssignment==i)))'],[z(tunnelAssignment==i)],'o')
    plot([c],[mean(z(tunnelAssignment==i))],'k*-','LineWidth',3)
end
set(gca,'XTick',[1:length(gens)])
%set(gca,'XTickLabel',(genotypes))
set(gca,'XTickLabel',genotypes(gens))
ylabel('preference')
xtickangle(30)
set(gca,'FontSize',fs)
axis([0 length(gens)+1  ymin ymax])
subplot(1,2,2)
hold on
gens=[2 8 6];
c=0;
for i = gens
    c=c+1;
    plot([(c)*ones(1,length(find(tunnelAssignment==i)))'],[z(tunnelAssignment==i)],'o')
    plot([c],[mean(z(tunnelAssignment==i))],'k*-','LineWidth',3)
end
set(gca,'XTick',[1:length(gens)])
%set(gca,'XTickLabel',(genotypes))
set(gca,'XTickLabel',genotypes(gens))
ylabel('preference')
xtickangle(30)
set(gca,'FontSize',fs)
axis([0 length(gens)+1  ymin ymax])


figure
hold on
c=0;
gens=[2 10 9];
for i = gens
    c=c+1;
    plot([(c)*ones(1,length(find(tunnelAssignment==i)))'],[z(tunnelAssignment==i)],'o')
    plot([c],[mean(z(tunnelAssignment==i))],'k*-','LineWidth',3)
end
set(gca,'XTick',[1:length(gens)])
%set(gca,'XTickLabel',(genotypes))
set(gca,'XTickLabel',genotypes(gens))
ylabel('preference')
xtickangle(30)
set(gca,'FontSize',fs)
axis([0 length(gens)+1  ymin ymax])


figure
subplot(1,2,1)
hold on
gens=[1 11 13];
c=0;
for i = gens
    c=c+1;
    plot([(c)*ones(1,length(find(tunnelAssignment==i)))'],[z(tunnelAssignment==i)],'o')
    plot([c],[mean(z(tunnelAssignment==i))],'k*-','LineWidth',3)
end
set(gca,'XTick',[1:length(gens)])
%set(gca,'XTickLabel',(genotypes))
set(gca,'XTickLabel',genotypes(gens))
ylabel('preference')
set(gca,'FontSize',fs)
xtickangle(30)
axis([0 length(gens)+1 ymin ymax])
subplot(1,2,2)
hold on
gens=[2 12 13];
c=0;
for i = gens
    c=c+1;
    plot([(c)*ones(1,length(find(tunnelAssignment==i)))'],[z(tunnelAssignment==i)],'o')
    plot([c],[mean(z(tunnelAssignment==i))],'k*-','LineWidth',3)
end
set(gca,'XTick',[1:length(gens)])
%set(gca,'XTickLabel',(genotypes))
set(gca,'XTickLabel',genotypes(gens))
ylabel('preference')
xtickangle(30)
set(gca,'FontSize',fs)
axis([0 length(gens)+1  ymin ymax])

figure
subplot(1,2,1)
hold on
gens=[1 15 16];
c=0;
for i = gens
    c=c+1;
    plot([(c)*ones(1,length(find(tunnelAssignment==i)))'],[z(tunnelAssignment==i)],'o')
    plot([c],[mean(z(tunnelAssignment==i))],'k*-','LineWidth',3)
end
set(gca,'XTick',[1:length(gens)])
%set(gca,'XTickLabel',(genotypes))
set(gca,'XTickLabel',genotypes(gens))
ylabel('preference')
set(gca,'FontSize',fs)
xtickangle(30)
axis([0 length(gens)+1  ymin ymax])
subplot(1,2,2)
hold on
gens=[2 14 16];
c=0;
for i = gens
    c=c+1;
    plot([(c)*ones(1,length(find(tunnelAssignment==i)))'],[z(tunnelAssignment==i)],'o')
    plot([c],[mean(z(tunnelAssignment==i))],'k*-','LineWidth',3)
end
set(gca,'XTick',[1:length(gens)])
%set(gca,'XTickLabel',(genotypes))
set(gca,'XTickLabel',genotypes(gens))
ylabel('preference')
xtickangle(30)
set(gca,'FontSize',fs)
axis([0 length(gens)+1  ymin ymax])





