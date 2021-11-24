function out=odorModel2d(m,z)

odornames=[{'air'},{'3-octanol'},{'1-hexanol'},{'ethyl lactate'},{'citronella'},{'2-heptanone'},{'1-pentanol'},{'ethanol'},{'geranyl acetate'},{'hexyl actetate'},{'4-methylcyclohexanol'},{'pentyl acetate'},{'1-butanol'}];
numX=size(m,2);
numY=size(m,1);

xBins=linspace(0,1,numX+1); xDel=xBins(2)-xBins(1);
yBins=1-linspace(0,1,numY+1); yDel=yBins(1)-yBins(2);

nyidalur=interp1([1 128 129 256],[1 0 1; .2559 .248 .2559; .248 .2559 .248; 0 1 0],1:256);

m=m*128+128;
m=round(m);
m(m<1)=1;
m(m>256)=256;

z=abs(z);

figure;
hold on;

for i=1:numX
    for j=1:numY
        rectangle('Position',[xBins(i) yBins(j) xDel yDel], 'FaceColor', [nyidalur(m(j,i),:) min([1 z(j,i)/2])],'LineStyle','none');
        
    end
end

axis([min(xBins) max(xBins) min(yBins)+yDel max(yBins)+yDel])
set(gca,'XTick',[xBins]+xDel/2)
set(gca,'XTickLabel',[{'D'},{'DC2'},{'DL5'},{'DM1'},{'DM2'},{'DM3'}])
xtickangle(30)
set(gca,'YTick',[yBins(end:-1:1)]+1.5*yDel)
set(gca,'YTickLabel',string(odornames(end:-1:1)))
ytickangle(30)
set(gca,'FontSize',12)
%colorbar