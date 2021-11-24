figure

for zslice=1:size(labels,3)
     subplot(2,3,1)
    imagesc(n(:,:,zslice))
     subplot(2,3,2)
    imagesc(b(:,:,zslice))
    
    subplot(2,3,3)
    imagesc(labels(:,:,zslice))
    
    subplot(2,3,4)
    imagesc(labels2(:,:,zslice))
    
        subplot(2,3,5)
    imagesc(labels3(:,:,zslice))
    
        subplot(2,3,6)
    imagesc(labels4(:,:,zslice))
    drawnow
    pause(0.1)
end