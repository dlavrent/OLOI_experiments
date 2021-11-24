close all
figure
for i = 1:size(l,3)
    imagesc(l(:,:,i),[0 max(max(max(l)))])
    pause(0.01)
end
%%
%close all
%figure
for i = 1:size(maskedSpotIm,3)
    imagesc(maskedSpotIm(:,:,i))
    pause(0.05)
end

%%
close all
figure
for i = 1:size(mask,3)
    imagesc(mask(:,:,i))
    pause(0.05)
end
%%
close all
figure
for i = 1:size(rawim,3)
    imagesc(rawim(:,:,i))
    pause(0.05)
end

%%
close all
figure
for i = 1:size(rawim,3)
    subplot(1,2,1)
    imagesc(rawim(:,:,i))
    subplot(1,2,2)
    imagesc(maskedSpotIm(:,:,i))
    
    drawnow
    pause(0.05)
end