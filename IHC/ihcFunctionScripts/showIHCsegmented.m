function [] = showIHCsegmented(filename)

load(filename)
figure
view(3);
axis tight
camlight
lighting gouraud

hold on
mycmap=distinguishable_colors(length(glomeruliSegs));

for i=[1 2 3 5]
    if i == 1
        currcolor = 'y';
    elseif i==2
        currcolor = 'b';
    elseif i==3
        currcolor = 'g';
    elseif i==4
        currcolor = 'p';
    elseif i==5
        currcolor = 'r';
    end
    p2=patch(isosurface(glomeruliSegs{1,i}),'FaceColor',currcolor,'EdgeColor','none','FaceAlpha',0.3);
    isonormals(glomeruliSegs{1,i},p2)
    
    ctrL=regionprops(glomeruliSegs{1,i},'Centroid');
    ctrL=ctrL.Centroid;
    %text(ctrL(1),ctrL(2),ctrL(3),[glomeruliToSegment{i} 'L'],'FontSize',15,'FontWeight','Bold')
    
    p2=patch(isosurface(glomeruliSegs{2,i}),'FaceColor',currcolor,'EdgeColor','none','FaceAlpha',0.3);
    isonormals(glomeruliSegs{2,i},p2)
    
    ctrR=regionprops(glomeruliSegs{2,i},'Centroid');
    ctrR=ctrR.Centroid;
    %text(ctrR(1),ctrR(2),ctrR(3),[glomeruliToSegment{i} 'R'],'FontSize',15,'FontWeight','Bold')
    
end