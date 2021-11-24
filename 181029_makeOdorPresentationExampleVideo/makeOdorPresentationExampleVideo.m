clear all

load 181012_fly2_rL_volumes_3.mat

%%
close all
figure

ctrst=[0 400];
imageLength=25; % seconds
zDepth=70; % um
nRepeats=10; % times to repeat video
% Prepare the new file.
    vidObj = VideoWriter('odorPres_ethyllactate.avi');
    vidObj.FrameRate = 10;
    open(vidObj);
set(gca,'nextplot','replacechildren');

 set(gcf, 'Position', [100, 50, 1000, 400])
for repeat=1:nRepeats
    for i=1:size(greenChannel,4)
        for j=1:10
            subplot(2,5,j)
            imagesc(greenChannel(:,:,j,i),ctrst)
            colormap(hot)
            if j==1
               text(10,10,[num2str((i-1)*imageLength/size(greenChannel,4),'%4.2f') ' seconds'],'FontSize',10,'Color',[1 1 1])
               if i<6
               text(-50,-30,['Odor Off'],'FontSize',15)
               elseif i>=6 && i<10
                   text(-50,-30,['Odor On'],'FontSize',15)
               else
                   text(-50,-30,['Odor Off'],'FontSize',15)
               end
            end
            title(['z = ' num2str(round(zDepth-(j-1)*zDepth/(size(greenChannel,3)))) ' um'])
            set(gca,'XTick','')
            set(gca,'YTick','')
        end
        drawnow
        currFrame = getframe(gcf);
        writeVideo(vidObj,currFrame);
        pause(0.25)
    end
end
close(vidObj);