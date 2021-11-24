%% modulate measured r
clear all
close all
rs=linspace(0,1,50);
truer=0.4;
Niters=1000;
for i=1:length(rs)
temp=sampleCorrelated(rs(i),27,Niters);
out(i)=sum(temp.sampleStats(:,1)>truer);
disp(['iter ' num2str(i) ' of ' num2str(length(rs))])
end

figure
plot(rs,out/Niters,'k','LineWidth',3)
xlabel('?true? r of Ca++-behavior model')
ylabel('p-value')
set(gca,'FontSize',15)
box off

%% modulate measured r
clear all
close all
rs=[5:50];
truer=0.4;
Niters=1000;
for i=1:length(rs)
temp=sampleCorrelated(truer,rs(i),Niters);
out(i)=sum(temp.sampleStats(:,1)>0.1732);
disp(['iter ' num2str(i) ' of ' num2str(length(rs))])
end

figure
plot(rs,out/Niters,'k','LineWidth',3)
xlabel('n')
ylabel('p-value')
set(gca,'FontSize',15)
box off

