function [] = labelAutoSegmentation_orcobrpshort(filename)
% loads filename and user selects clusters corresponding to glomeruli of
% interest
% assumes only labels and labels 2 were saved by batch processing script
load(filename)
savename=[filename(1:end-4) '_labelled'];
close all

glomeruliToSegment=[{'DM2'},{'DM3'},{'DC2'},{'DL5'},{'DM1'}];
glomeruliSegs=cell(2,length(glomeruliToSegment));

% randomize labels to help visualize
labels2orig=labels2;
labels2=zeros(size(labels2orig));
rp=randperm(max(labels2orig(:)));
for i=1:max(labels2orig(:))
    labels2(labels2orig==rp(i))=i;
end

% use two options for segmentation: one with low and one with high
% over-clustering
segmentation1=labels;
segmentation2=labels2;

minCtrst=-10;
ctrst=[minCtrst 1000];

for glomsToSegment=1:length(glomeruliToSegment)
    for leftright=[1 2]
        
        currentSeg=zeros(size(maskedNnorm));
        
        fig = gcf;
        ax = axes('Parent', fig);
        
        k=1;
        was_a_key=1;
        while k
            
            if was_a_key
                text(0,0,['z slice: ' num2str(k)],'FontSize',20)
                
                subplot(1,3,1)
                imagesc(maskedNnorm(:,:,k),ctrst)
                if leftright==1
                    title(['select glomerulus: ' glomeruliToSegment{glomsToSegment} ' left'])
                else
                    title(['select glomerulus: ' glomeruliToSegment{glomsToSegment} ' right'])
                end
                
                subplot(1,3,2)
                imagesc(segmentation1(:,:,k))
                set(gca,'tag',num2str(1))
                title('labels')
                
                subplot(1,3,3)
                imagesc(segmentation2(:,:,k))
                set(gca,'tag',num2str(2))
                title('labels2')
                
                drawnow
            end
            
            was_a_key = waitforbuttonpress;
            if was_a_key && (strcmp(get(fig, 'CurrentKey'), 'downarrow') || strcmp(get(fig, 'CurrentKey'), 'leftarrow'))
                if k>1
                    k = k - 1;
                else
                    k=1;
                end
            elseif was_a_key && (strcmp(get(fig, 'CurrentKey'), 'uparrow') || strcmp(get(fig, 'CurrentKey'), 'rightarrow'))
                if k<size(n,3)
                    k = k + 1;
                else
                    k=size(n,3);
                end
            elseif was_a_key  && get(gcf,'currentcharacter')==27 % 27 = escape key
                k=0;
                disp('current label recorded')
            elseif ~was_a_key
                [xi,yi] = ginput(1); xi=round(xi);yi=round(yi);
                mousept = get(gca,'currentPoint');
                xs = mousept(1,1);
                ys = mousept(1,2);
                subplotnum=get(gca,'tag');
                
                if str2num(subplotnum)==1
                    selectedRegion=segmentation1(yi,xi,k);
                    currentSeg(segmentation1==selectedRegion)=1;
                    segmentation1(segmentation1==selectedRegion)=minCtrst;
                elseif str2num(subplotnum)==2
                    selectedRegion=segmentation2(yi,xi,k);
                    currentSeg(segmentation2==selectedRegion)=1;
                    segmentation2(segmentation2==selectedRegion)=minCtrst;
                end
            end
        end
        
%         figure
%         view(3);
%         axis tight
%         camlight
%         lighting gouraud
%         
%         p2=patch(isosurface(currentSeg),'FaceColor',[0.8 0 0.2],'EdgeColor','none','FaceAlpha',0.3);
%         isonormals(currentSeg,p2)
        
        glomeruliSegs{leftright,glomsToSegment}=currentSeg;
        
    end
end

figure
view(3);
axis tight
camlight
lighting gouraud

hold on
mycmap=distinguishable_colors(length(glomeruliSegs));

for i=1:length(glomeruliSegs)
    p2=patch(isosurface(glomeruliSegs{1,i}),'FaceColor',mycmap(i,:),'EdgeColor','none','FaceAlpha',0.3);
    isonormals(glomeruliSegs{1,i},p2)
    
    
    ctrL=regionprops(glomeruliSegs{1,i},'Centroid');
    ctrL=ctrL.Centroid;
    text(ctrL(1),ctrL(2),ctrL(3),[glomeruliToSegment{i} 'L'],'FontSize',15,'FontWeight','Bold')
    
    p2=patch(isosurface(glomeruliSegs{2,i}),'FaceColor',mycmap(i,:),'EdgeColor','none','FaceAlpha',0.3);
    isonormals(glomeruliSegs{2,i},p2)
    
    ctrR=regionprops(glomeruliSegs{2,i},'Centroid');
    ctrR=ctrR.Centroid;
    text(ctrR(1),ctrR(2),ctrR(3),[glomeruliToSegment{i} 'R'],'FontSize',15,'FontWeight','Bold')
    
end

save(savename,'glomeruliSegs','glomeruliToSegment')
disp('saved labels!')