clear all
close all

load 191009
zz=z;
zzv=zv;
zzp=zp;
tt=t;
ttv=tv;
ttp=tp;
load 191010
zz=[zz; z];
tt=[tt; t];
zzv=[zzv; zv];
ttv=[ttv; tv];
zzp=[zzp; zp];
ttp=[ttp; tp];

z=zz;
t=tt;
zp=zzp;
tp=ttp;
zv=zzv;
tv=ttv;

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

%%
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