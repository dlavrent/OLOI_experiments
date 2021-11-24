clear all
close all

% from hallem and carlson 2004
fruitnames=[{'apple'},{'apricot'},{'banana'},{'cherry'},{'mango'},{'peach'},{'pineapple'},{'raspberry'},{'strawberry'}];
fruit=[40 210 126 28 73 248 68 6 130 76 158 250 -13 20 165 15 109 158 40 87 203 112 49 135; ...
    10 165 71 5 56 139 30 -3 40 47 48 115 -16 16 146 5 94 78 28 26 55 70 34 42; ...
    55 203 255 117 126 245 77 10 89 112 247 238 -26 26 158 26 209 213 28 80 213 97 19 175; ...
    5 206 42 43 47 130 25 -8 83 68 51 94 -15 91 65 9 170 83 26 25 23 68 9 4; ...
    15 143 75 49 104 243 39 -2 60 51 95 108 -10 12 106 9 151 90 28 162 94 67 21 98; ...
    16 159 84 1 47 161 27 -5 78 52 49 106 -13 26 135 7 108 76 31 24 31 60 8 20; ...
    13 136 133 30 38 234 27 -3 51 59 82 144 06 14 151 9 125 84 29 85 92 67 7 23; ...
    19 181 109 42 63 234 57 -3 78 44 59 114 -23 24 140 13 119 97 42 30 195 80 -6 98; ...
    35 197 157 116 47 237 67 -3 90 51 126 183 -19 27 172 18 171 149 27 123 196 77 -9 32];

fruitcorr=corr(fruit');

figure
imagesc(fruitcorr,[0 1])
set(gca,'XTick',[1:length(fruitnames)])
set(gca,'XTickLabel',fruitnames)
xtickangle(30)
set(gca,'YTick',[1:length(fruitnames)])
set(gca,'YTickLabel',fruitnames)
set(gca,'FontSize',15)
h=colorbar;
title('response correlation (24 ORNs)')
title(h,'r')




fdir='/Users/mattchurgin/Dropbox (Harvard University)/flyimaging/analysis/modelAntennaLobe';
fname='hallemAndCarlson2004_tableS1data.xlsx';
[NUM,TXT,RAW]=xlsread([fdir '/' fname]);

ORNr=NUM(1:110,:)'+NUM(111,:)';  % add background firing rate for each OR

otherodorcorr=corr(ORNr);

figure
imagesc(otherodorcorr,[0 1])
set(gca,'FontSize',15)
xlabel('odor')
ylabel('odor')
h=colorbar;
title('response correlation (24 ORNs)')
title(h,'r')


figure
imagesc(corr(ORNr'))
set(gca,'XTick',[1:length(TXT)])
set(gca,'XTickLabel',TXT)
xtickangle(30)
set(gca,'YTick',[1:length(TXT)])
set(gca,'YTickLabel',TXT)
set(gca,'FontSize',15)
h=colorbar;
title('response correlation (24 ORNs)')
title(h,'r')


% calculate euclidean distance between fruit odors
fruitdist = squareform(pdist(fruit));
allodordist = squareform(pdist(ORNr'));

figure
imagesc(fruitdist)
set(gca,'XTick',[1:length(fruitnames)])
set(gca,'XTickLabel',fruitnames)
xtickangle(30)
set(gca,'YTick',[1:length(fruitnames)])
set(gca,'YTickLabel',fruitnames)
set(gca,'FontSize',15)
h=colorbar;
title('coding distance (24 ORNs)')
title(h,'distance')

%%
%plot mean distances
cmap = hsv(9);
mf = mean(fruitdist);
figure
for i = 1:9
bar(i,mf(i),'FaceColor',cmap(i,:))
hold on
end
errorbar(1:9,mf,std(fruitdist)/sqrt(9),'k.','LineWidth',3)
set(gca,'XTick',[1:length(fruitnames)])
set(gca,'XTickLabel',fruitnames)
xtickangle(30)
ylabel('inter-odor coding distance')
set(gca,'Fontsize',15)