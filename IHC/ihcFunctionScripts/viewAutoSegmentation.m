close all



% randomize labels to help visualize
%
% labels1orig=labels;
% labels=zeros(size(labels1orig));
% rp=randperm(max(labels1orig(:)));
% for i=1:max(labels1orig(:))
%    labels(labels1orig==rp(i))=i;
% end

labels2orig=labels2;
labels2=zeros(size(labels2orig));
rp=randperm(max(labels2orig(:)));
for i=1:max(labels2orig(:))
    labels2(labels2orig==rp(i))=i;
end
%
% labels3orig=labels3;
% labels3=zeros(size(labels3orig));
% rp=randperm(max(labels3orig(:)));
% for i=1:max(labels3orig(:))
%    labels3(labels3orig==rp(i))=i;
% end

% labels4orig=labels4;
% labels4=zeros(size(labels4orig));
% rp=randperm(max(labels4orig(:)));
% for i=1:max(labels4orig(:))
%     labels4(labels4orig==rp(i))=i;
% end
% 
% labels5orig=labels5;
% labels5=zeros(size(labels5orig));
% rp=randperm(max(labels5orig(:)));
% for i=1:max(labels5orig(:))
%     labels5(labels5orig==rp(i))=i;
% end



segmentation1=labels2;
segmentation2=labels;

fig = gcf;
ax = axes('Parent', fig);
k=1;
was_a_key=1;
ctrst=[0 1000];
while k
    
    if was_a_key
        text(0,0,['z slice: ' num2str(k)],'FontSize',20)
        
        
        subplot(1,3,1)
        imagesc(maskedNnorm(:,:,k),ctrst)
        title('n')
        
        subplot(1,3,2)
        imagesc(segmentation1(:,:,k))
        set(gca,'tag',num2str(1))
        
        subplot(1,3,3)
        imagesc(segmentation2(:,:,k))
        set(gca,'tag',num2str(2))
   
        
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
    elseif was_a_key  && strcmp(get(fig, 'CurrentKey'), 'x')
        k=0;
        disp('exit')
    end
end