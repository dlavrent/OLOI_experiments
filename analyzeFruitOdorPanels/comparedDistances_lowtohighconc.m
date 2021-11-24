clear all 
close all
odornames=[{'air'},{'3-octanol'},{'1-hexanol'},{'ethyl lactate'},{'balsamic vinegar'},{'2-heptanone'},{'1-pentanol'},{'apple cider vinegar'},{'pineapple'},{'hexyl acetate'},{'4-methylcyclohexanol'},{'pentyl acetate'},{'1-butanol'}];
load lowConc_distances
ld = averaged;

load highConc_distances
hd = averaged;

r = hd./ld;
imagesc(r)
set(gca,'XTick',[1:length(odornames)])
set(gca,'XTickLabel',odornames)
xtickangle(30)
set(gca,'YTick',[1:length(odornames)])
set(gca,'YTickLabel',odornames)
set(gca,'FontSize',15)
hcb = colorbar;
title(hcb,'distance ratio')
