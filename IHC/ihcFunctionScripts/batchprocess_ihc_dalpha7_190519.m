
clear all
close all

startDir='E:\IHC_and_behavior';
cd(startDir)

saveDest='C:\Users\mac0456\Dropbox\flyimaging\analysis\IHC\Gh146Dalpha7\autoSegmentation';


folderNames{1}='190919_gh146dalpha7_behaviorIHC_age10days';
flyNums{1}=[8 12 15];


for currFolderNames=1:length(folderNames)
    disp(['starting to analyze folder ' folderNames{currFolderNames}])
    cd(folderNames{currFolderNames})
    for currFlyNum=flyNums{currFolderNames}
        close all
        disp(['analyzing folder ' folderNames{currFolderNames} '. fly ' num2str(currFlyNum)])
        % 1. load images
        try
            filePrefix=['fly' num2str(currFlyNum) '_'];
            nc82name=[filePrefix 'nc82.tif'];
            brpname=[filePrefix 'dalpha7.tif'];
            mcd8name=[filePrefix 'myr.tif'];
            
            % get image info
            nc82info=imfinfo(nc82name);
        catch % some fly names have 2 decimal places
            filePrefix=['fly' num2str(currFlyNum,'%02d') '_'];
            nc82name=[filePrefix 'nc82.tif'];
            brpname=[filePrefix 'dalpha7.tif'];
            mcd8name=[filePrefix 'myr.tif'];
            
            % get image info
            nc82info=imfinfo(nc82name);
        end
        
        nc82l=length(nc82info);
        
        nc82im=zeros(nc82info(1).Height,nc82info(1).Width,nc82l);
        brpim=zeros(nc82info(1).Height,nc82info(1).Width,nc82l);
        mcd8im=zeros(nc82info(1).Height,nc82info(1).Width,nc82l);
        for i=1:nc82l
            
            nc82im(:,:,i)=imread(nc82name,'Index',i);
            brpim(:,:,i)=imread(brpname,'Index',i);
            mcd8im(:,:,i)=imread(mcd8name,'Index',i);
            
            if mod(i,10)==0
                disp(['loaded ' num2str(i)])
            end
        end
        disp('images loaded')
        
        
        
        % 2. median filter and normalize images
        
        zs=1:size(nc82im,3);
        for xystep=[13]
            tic
            
            n=nc82im(:,:,zs);
            b=brpim(:,:,zs);
            m=mcd8im(:,:,zs);
            
            %xystep=9; % default 13
            
            xystart=(xystep-1)/2;
            n=medfilt3(n, [xystep xystep 1]); n=n(xystart:xystep:end,xystart:xystep:end,:);
            disp('median filtering step 1 complete')
            b=medfilt3(b, [xystep xystep 1]); b=b(xystart:xystep:end,xystart:xystep:end,:);
            disp('median filtering step 2 complete')
            m=medfilt3(m, [xystep xystep 1]); m=m(xystart:xystep:end,xystart:xystep:end,:);
            disp('median filtering step 3 complete')
            
            
            gs=30;
            x=round([-10:10]*13/xystep); % depends on xystep (smaller range for larger xystep). correction factor is round(13/xystep)
            y=x;
            z=-1:1;
            [xx yy zz]=meshgrid(x,y,z);
            gau=1/(length(x)*gs^2)*exp(-(xx.^2+yy.^2 +zz.^2)/(2*gs^2));
            disp('normalizing nc82 image')
            
            nNorm=n./(convn(n,gau,'same'));
            bNorm=b./(convn(b,gau,'same'));
            mNorm=m./(convn(m,gau,'same'));
            
            lmask=ones(size(nNorm));
            disp(['done normalizing. time elapsed: ' num2str(toc/60) ' minutes'])
            
            % 3.  mask by expression pattern of driver lines
            
            driver=m.*b;
            dm=driver;
            driverThresh=75;
            dm(dm<prctile(dm(:),driverThresh))=0;
            dm(dm>0)=1;
            dm=logical(dm);
            
            
            % 4.threshold and watershed
            maskedNnorm=nNorm.*dm.*mNorm;
            
            
            %thresh=prctile(maskedNnorm(:),90);
            
            asdf=maskedNnorm(:);
            maskednormthresh=55;
            thresh=prctile(asdf(asdf>0),maskednormthresh);
            
            im=maskedNnorm<thresh;
            im=im.*lmask;
            imb=bwdist(im);
            imb=imb.*lmask;
            
            % Watershed transform
            minSuppressionThreshold=1;
            preL=max(imb(:))-imb;
            preL2=imhmin(preL,minSuppressionThreshold);
            L = watershed(preL2);
            labels=bwlabeln(L&~im);
            labels=labels.*lmask;
            
            disp('finished watershed')
            
            % 6. save initial watershed results
            
            warning('off')
            areathresh=300;
            areathreshhigh=600000;
            
            lbls=regionprops(labels,'all');
            todelete=[];
            clusterVols=cell(1,length(lbls));
            clusterInfo=cell(1,length(lbls));
            bb=zeros(length(lbls),7);
            for j=1:length(lbls)
                if lbls(j).Area>areathresh && lbls(j).Area<areathreshhigh
                    clusterVols{j}=(labels==j);
                    clusterInfo{j}=lbls(j);
                    
                    bb(j,:)=[clusterInfo{j}.BoundingBox clusterInfo{j}.FilledArea];
                else
                    labels(labels==j)=0;
                    todelete=[todelete j];
                end
                if mod(j,100)==0
                    disp(['post-processed ' num2str(j)])
                end
            end
            u=unique(labels(:));
            for i=1:length(u)
                labels(labels(:)==u(i))=i;
            end
            clusterVols(todelete)=[];
            clusterInfo(todelete)=[];
            bb(todelete,:)=[];
            
            
            % 7. apply 3d convolution plus watershed to split merged gloms
            splitInfo=cell(1,length(clusterVols));
            labels2=zeros(size(clusterVols{1}));
            clear clusterVols2 clusterInfo2
            newcno=1;
            for i=1:length(clusterVols)
                try
                ballr=round(nthroot(clusterInfo{i}.Area,3)/4);
                tempBall= strel('sphere',ballr);
                tempconv=convn(maskedNnorm.*clusterVols{i},tempBall.Neighborhood/sum(sum(sum(tempBall.Neighborhood))),'same');
                tempconv(isnan(tempconv))=0;
                
                tcb=tempconv<prctile(tempconv(tempconv~=0),50);
                tcb2=bwdist(tcb);
                
                % Watershed transform
                minSuppressionThreshold=0.25;
                preL=max(tcb2(:))-tcb2;
                preL2=imhmin(preL,minSuppressionThreshold);
                L = watershed(preL2);
                labelm=bwlabeln(L&~tcb);
                
                r1=regionprops(labelm);
                splitInfo{i}=r1;
                
                if length(splitInfo{i})>1
                    p=clusterInfo{i}.PixelList;
                    pidx=clusterInfo{i}.PixelIdxList;
                    clear tempd
                    for kkk=1:length(splitInfo{i})
                        tempc(kkk,:)=splitInfo{i}(kkk).Centroid;
                        tempd(kkk,:)=sqrt(sum((p-tempc(kkk,:)).^2,2));
                    end
                    
                    [v ind]=min(tempd);
                    
                    for kkk=1:length(splitInfo{i})
                        newv1=logical(zeros(size(clusterVols{i})));
                        newv1(pidx(ind==kkk))=1;
                        clusterVols2{newcno}=newv1;
                        clusterInfo2{newcno}=regionprops(newv1,'all');
                        labels2(clusterVols2{newcno}(:)==1)=newcno;
                        newcno=newcno+1;
                    end
                else
                    clusterVols2{newcno}=clusterVols{i};
                    clusterInfo2{newcno}=clusterInfo{i};
                    labels2(clusterVols2{newcno}(:)==1)=newcno;
                    newcno=newcno+1;
                end
                
                disp(['cluster num: ' num2str(i) '. new cluster num: ' num2str(newcno)])
                catch
                end
            end
            
            % 7. perform pca on pixels
            
            % mask out pixels below 50th percentile of combined b and m
            combined=b.*m;
            tmask=double(combined>prctile(combined(:),50));
            tmask(tmask==0)=NaN;
            
            nm=n.*tmask;
            bm=b.*tmask;
            mm=m.*tmask;
            
            nm(nm==0)=NaN;
            mm(mm==0)=NaN;
            bm(bm==0)=NaN;
            
            % set up data matrix
            % make xyz mesh
            [yyy xxx zzz]= ndgrid(1:size(nm,1),1:size(nm,2),1:size(nm,3));
            
            % unwrap images for pca
            data = [bm(:)'; mm(:)'; transpose(bm(:).*mm(:)); xxx(:)'; yyy(:)'; zzz(:)']';
            data=(data-nanmean(data))./nanstd(data);
            
            interactionweight=0.75;
            xyweight=1;
            zweight=xyweight;
            
            data(:,3)=interactionweight*data(:,3);
            data(:,(end-2):(end-1))=xyweight*data(:,(end-2):(end-1));
            data(:,end)=zweight*data(:,end);
            
            
            % disp('finished preprocessing')
            
            [coeff, score, latent, tsquared, explained] = pca(data);
            
            % 8. split clusters based on pixel similarity
            labels3=labels;
            
            nomoresplits=1;
            twosplits=0;
            
            numwhileloops=0;
            while numwhileloops<2%nomoresplits
                %while nomoresplits
                numwhileloops=numwhileloops+1;
                splitthistime=0;
                for ws=1:max(labels3(:))
                    temp=score(labels3(:)==ws,:);
                    if sum(any(isfinite(temp)'))>2
                        
                        [ktemp ctemp sumtemp]=kmeans(temp,1);
                        [ktemp2 ctemp2 sumtemp2]=kmeans(temp,2);
                        
                        % determine if one or two clusters is a better fit for current cluster
                        onekmeansum(ws)=sumtemp/length(ktemp);
                        twokmeansum(ws)=sumtemp2(1)/sum(ktemp2==1) + sumtemp2(2)/sum(ktemp2==2);
                        
                        
                        if twokmeansum(ws)<onekmeansum(ws)
                            splitthistime=1;
                            pixels=find(labels3==ws);
                            labels3(pixels(ktemp2==1))=ws;
                            labels3(pixels(ktemp2==2))=max(labels3(:))+1;
                            twosplits=twosplits+1;
                        end
                    end
                    %disp(['cluster ' num2str(ws)])
                end
                
                if splitthistime==0
                    nomoresplits=0;
                end
                disp(['while loop iteration #' num2str(numwhileloops)])
            end
            
            warning('off')
            areathresh=500;
            areathreshhigh=600000;
            
            lbls=regionprops(labels3,'all');
            todelete=[];
            clusterVols3=cell(1,length(lbls));
            clusterInfo3=cell(1,length(lbls));
            bb=zeros(length(lbls),7);
            for j=1:length(lbls)
                if lbls(j).Area>areathresh && lbls(j).Area<areathreshhigh
                    clusterVols3{j}=(labels3==j);
                    clusterInfo3{j}=lbls(j);
                    
                    bb(j,:)=[clusterInfo3{j}.BoundingBox clusterInfo3{j}.FilledArea];
                else
                    labels3(labels3==j)=0;
                    todelete=[todelete j];
                end
                if mod(j,100)==0
                    disp(['post-processed ' num2str(j)])
                end
            end
            u=unique(labels3(:));
            for i=1:length(u)
                labels3(labels3(:)==u(i))=i;
            end
            clusterVols3(todelete)=[];
            clusterInfo3(todelete)=[];
            bb(todelete,:)=[];
            
            %
            % labelsEdited=newlabels;
            % for i=1:length(todelete)
            %     labelsEdited(labelsEdited==todelete(i))=0;
            % end
            
            % 9. use 3d convolution plus watershed
            splitInfo=cell(1,length(clusterVols3));
            labels4=zeros(size(clusterVols3{1}));
            clear clusterVols4 clusterInfo4
            newcno=1;
            for i=1:length(clusterVols3)
                try
                ballr=round(nthroot(clusterInfo3{i}.Area,3)/4);
                tempBall= strel('sphere',ballr);
                tempconv=convn(maskedNnorm.*clusterVols3{i},tempBall.Neighborhood/sum(sum(sum(tempBall.Neighborhood))),'same');
                tempconv(isnan(tempconv))=0;
                
                tcb=tempconv<prctile(tempconv(tempconv~=0),50);
                tcb2=bwdist(tcb);
                
                % Watershed transform
                minSuppressionThreshold=0.25;
                preL=max(tcb2(:))-tcb2;
                preL2=imhmin(preL,minSuppressionThreshold);
                L = watershed(preL2);
                labelm=bwlabeln(L&~tcb);
                
                %labelm=labels.*clusterVols2{i};
                
                r1=regionprops(labelm);
                splitInfo{i}=r1;
                
                if length(splitInfo{i})>1
                    p=clusterInfo3{i}.PixelList;
                    pidx=clusterInfo3{i}.PixelIdxList;
                    clear tempd
                    for kkk=1:length(splitInfo{i})
                        tempc(kkk,:)=splitInfo{i}(kkk).Centroid;
                        tempd(kkk,:)=sqrt(sum((p-tempc(kkk,:)).^2,2));
                    end
                    
                    [v ind]=min(tempd);
                    
                    for kkk=1:length(splitInfo{i})
                        newv1=logical(zeros(size(clusterVols3{i})));
                        newv1(pidx(ind==kkk))=1;
                        clusterVols4{newcno}=newv1;
                        clusterInfo4{newcno}=regionprops(newv1,'all');
                        labels4(clusterVols4{newcno}(:)==1)=newcno;
                        newcno=newcno+1;
                    end
                else
                    clusterVols4{newcno}=clusterVols3{i};
                    clusterInfo4{newcno}=clusterInfo3{i};
                    labels4(clusterVols4{newcno}(:)==1)=newcno;
                    newcno=newcno+1;
                end
                
                disp(['cluster num: ' num2str(i) '. new cluster num: ' num2str(newcno)])
                catch
                end
            end
            
            
            splitInfo=cell(1,length(clusterVols3));
            labels5=zeros(size(clusterVols3{1}));
            clear clusterVols5 clusterInfo5
            newcno=1;
            for i=1:length(clusterVols3)
                try
                ballr=round(nthroot(clusterInfo3{i}.Area,3)/8);
                tempBall= strel('sphere',ballr);
                tempconv=convn(maskedNnorm.*clusterVols3{i},tempBall.Neighborhood/sum(sum(sum(tempBall.Neighborhood))),'same');
                tempconv(isnan(tempconv))=0;
                
                tcb=tempconv<prctile(tempconv(tempconv~=0),50);
                tcb2=bwdist(tcb);
                
                % Watershed transform
                minSuppressionThreshold=0.25;
                preL=max(tcb2(:))-tcb2;
                preL2=imhmin(preL,minSuppressionThreshold);
                L = watershed(preL2);
                labelm=bwlabeln(L&~tcb);
                
                %labelm=labels.*clusterVols2{i};
                
                r1=regionprops(labelm);
                splitInfo{i}=r1;
                
                if length(splitInfo{i})>1
                    p=clusterInfo3{i}.PixelList;
                    pidx=clusterInfo3{i}.PixelIdxList;
                    clear tempd
                    for kkk=1:length(splitInfo{i})
                        tempc(kkk,:)=splitInfo{i}(kkk).Centroid;
                        tempd(kkk,:)=sqrt(sum((p-tempc(kkk,:)).^2,2));
                    end
                    
                    [v ind]=min(tempd);
                    
                    for kkk=1:length(splitInfo{i})
                        newv1=logical(zeros(size(clusterVols3{i})));
                        newv1(pidx(ind==kkk))=1;
                        clusterVols5{newcno}=newv1;
                        clusterInfo5{newcno}=regionprops(newv1,'all');
                        labels5(clusterVols5{newcno}(:)==1)=newcno;
                        newcno=newcno+1;
                    end
                else
                    clusterVols5{newcno}=clusterVols3{i};
                    clusterInfo5{newcno}=clusterInfo3{i};
                    labels5(clusterVols5{newcno}(:)==1)=newcno;
                    newcno=newcno+1;
                end
                
                disp(['cluster num: ' num2str(i) '. new cluster num: ' num2str(newcno)])
                catch
                end
            end
            
            if ~exist([saveDest '\' folderNames{currFolderNames}])
                mkdir([saveDest '\' folderNames{currFolderNames}])
            end
            save([saveDest '\' folderNames{currFolderNames} '\' filePrefix '_autoseg_xystep' num2str(xystep) '.mat'],'clusterVols','clusterInfo','clusterVols2','clusterInfo2','clusterVols3','clusterInfo3','clusterVols4','clusterInfo4','clusterVols5','clusterInfo5','n','b','m','labels','labels2','labels3','labels4','labels5','maskedNnorm','thresh')
            
            disp(['successfully saved data for fly ' num2str(currFlyNum) ' in folder ' folderNames{currFolderNames}])
        end
    end
    cd(startDir)
end


%% reprocess and overwrite existing files. only recalculates labels2 and labels4

clear all
close all

startDir='C:\Users\mac0456\Dropbox\flyimaging\analysis\IHC\Gh146Dalpha7\autoSegmentation';
cd(startDir)

saveDest='C:\Users\mac0456\Dropbox\flyimaging\analysis\IHC\Gh146Dalpha7\autoSegmentation';
xystep=13;

folderNames{1}='190301_gh146_dalpha7OXZ';
flyNums{1}=[];%8 18 19 23 26 28 30];
folderNames{2}='190304_gh146_dalpha7OXZ';
flyNums{2}=[]%;1 3 17 18 30 40 43 44];
folderNames{3}='190507_gh146dalpha7_behaviorPlusIHC';
flyNums{3}=[2 4 5 8 9 10 12 14 15];
folderNames{4}='190214_gh146_Dalpha7gfp_myr_behaviorANDIHC';
flyNums{4}=[2 5 8 23 28];

for currFolderNames=1:length(folderNames)
    disp(['starting to analyze folder ' folderNames{currFolderNames}])
    cd(folderNames{currFolderNames})
    for currFlyNum=flyNums{currFolderNames}
        close all
        disp(['analyzing folder ' folderNames{currFolderNames} '. fly ' num2str(currFlyNum)])
        
        try
            filePrefix=['fly' num2str(currFlyNum) '_'];
            currfname=[saveDest '\' folderNames{currFolderNames} '\' filePrefix '_autoseg_xystep' num2str(xystep) '.mat'];
            load(currfname)
        catch % some fly names have 2 decimal places
            filePrefix=['fly' num2str(currFlyNum,'%02d') '_'];
            currfname=[saveDest '\' folderNames{currFolderNames} '\' filePrefix '_autoseg_xystep' num2str(xystep) '.mat'];
            load(currfname)
        end
        
        
        % apply 3d convolution plus watershed to split merged gloms
        splitInfo=cell(1,length(clusterVols));
        labels2=zeros(size(clusterVols{1}));
        clear clusterVols2 clusterInfo2
        newcno=1;
        for i=1:length(clusterVols)
            try
            ballr=round(nthroot(clusterInfo{i}.Area,3)/4);
            tempBall= strel('sphere',ballr);
            tempconv=convn(maskedNnorm.*clusterVols{i},tempBall.Neighborhood/sum(sum(sum(tempBall.Neighborhood))),'same');
            tempconv(isnan(tempconv))=0;
            
            tcb=tempconv<prctile(tempconv(tempconv~=0),50);
            tcb2=bwdist(tcb);
            
            % Watershed transform
            minSuppressionThreshold=0.25;
            preL=max(tcb2(:))-tcb2;
            preL2=imhmin(preL,minSuppressionThreshold);
            L = watershed(preL2);
            labelm=bwlabeln(L&~tcb);
            
            r1=regionprops(labelm);
            splitInfo{i}=r1;
            warning('off')
            if length(splitInfo{i})>1
                p=clusterInfo{i}.PixelList;
                pidx=clusterInfo{i}.PixelIdxList;
                clear tempd
                for kkk=1:length(splitInfo{i})
                    tempc(kkk,:)=splitInfo{i}(kkk).Centroid;
                    tempd(kkk,:)=sqrt(sum((p-tempc(kkk,:)).^2,2));
                end
                
                [v ind]=min(tempd);
                
                for kkk=1:length(splitInfo{i})
                    newv1=logical(zeros(size(clusterVols{i})));
                    newv1(pidx(ind==kkk))=1;
                    clusterVols2{newcno}=newv1;
                    clusterInfo2{newcno}=regionprops(newv1,'all');
                    labels2(clusterVols2{newcno}(:)==1)=newcno;
                    newcno=newcno+1;
                end
            else
                clusterVols2{newcno}=clusterVols{i};
                clusterInfo2{newcno}=clusterInfo{i};
                labels2(clusterVols2{newcno}(:)==1)=newcno;
                newcno=newcno+1;
            end
            
            disp(['cluster num: ' num2str(i) '. new cluster num: ' num2str(newcno)])
            catch
            end
        end
        
        
        % 9. use 3d convolution plus watershed
        splitInfo=cell(1,length(clusterVols3));
        labels4=zeros(size(clusterVols3{1}));
        clear clusterVols4 clusterInfo4
        newcno=1;
        for i=1:length(clusterVols3)
            try
            ballr=round(nthroot(clusterInfo3{i}.Area,3)/4);
            tempBall= strel('sphere',ballr);
            tempconv=convn(maskedNnorm.*clusterVols3{i},tempBall.Neighborhood/sum(sum(sum(tempBall.Neighborhood))),'same');
            tempconv(isnan(tempconv))=0;
            
            tcb=tempconv<prctile(tempconv(tempconv~=0),50);
            tcb2=bwdist(tcb);
            
            % Watershed transform
            minSuppressionThreshold=0.25;
            preL=max(tcb2(:))-tcb2;
            preL2=imhmin(preL,minSuppressionThreshold);
            L = watershed(preL2);
            labelm=bwlabeln(L&~tcb);
            
            %labelm=labels.*clusterVols2{i};
            
            r1=regionprops(labelm);
            splitInfo{i}=r1;
            
            if length(splitInfo{i})>1
                p=clusterInfo3{i}.PixelList;
                pidx=clusterInfo3{i}.PixelIdxList;
                clear tempd
                for kkk=1:length(splitInfo{i})
                    tempc(kkk,:)=splitInfo{i}(kkk).Centroid;
                    tempd(kkk,:)=sqrt(sum((p-tempc(kkk,:)).^2,2));
                end
                
                [v ind]=min(tempd);
                
                for kkk=1:length(splitInfo{i})
                    newv1=logical(zeros(size(clusterVols3{i})));
                    newv1(pidx(ind==kkk))=1;
                    clusterVols4{newcno}=newv1;
                    clusterInfo4{newcno}=regionprops(newv1,'all');
                    labels4(clusterVols4{newcno}(:)==1)=newcno;
                    newcno=newcno+1;
                end
            else
                clusterVols4{newcno}=clusterVols3{i};
                clusterInfo4{newcno}=clusterInfo3{i};
                labels4(clusterVols4{newcno}(:)==1)=newcno;
                newcno=newcno+1;
            end
            disp(['cluster num: ' num2str(i) '. new cluster num: ' num2str(newcno)])
            catch
            end
        end
        
        splitInfo=cell(1,length(clusterVols3));
        labels5=zeros(size(clusterVols3{1}));
        clear clusterVols5 clusterInfo5
        newcno=1;
        for i=1:length(clusterVols3)
            try
            ballr=round(nthroot(clusterInfo3{i}.Area,3)/8);
            tempBall= strel('sphere',ballr);
            tempconv=convn(maskedNnorm.*clusterVols3{i},tempBall.Neighborhood/sum(sum(sum(tempBall.Neighborhood))),'same');
            tempconv(isnan(tempconv))=0;
            
            tcb=tempconv<prctile(tempconv(tempconv~=0),50);
            tcb2=bwdist(tcb);
            
            % Watershed transform
            minSuppressionThreshold=0.25;
            preL=max(tcb2(:))-tcb2;
            preL2=imhmin(preL,minSuppressionThreshold);
            L = watershed(preL2);
            labelm=bwlabeln(L&~tcb);
            
            %labelm=labels.*clusterVols2{i};
            
            r1=regionprops(labelm);
            splitInfo{i}=r1;
            
            if length(splitInfo{i})>1
                p=clusterInfo3{i}.PixelList;
                pidx=clusterInfo3{i}.PixelIdxList;
                clear tempd
                for kkk=1:length(splitInfo{i})
                    tempc(kkk,:)=splitInfo{i}(kkk).Centroid;
                    tempd(kkk,:)=sqrt(sum((p-tempc(kkk,:)).^2,2));
                end
                
                [v ind]=min(tempd);
                
                for kkk=1:length(splitInfo{i})
                    newv1=logical(zeros(size(clusterVols3{i})));
                    newv1(pidx(ind==kkk))=1;
                    clusterVols5{newcno}=newv1;
                    clusterInfo5{newcno}=regionprops(newv1,'all');
                    labels5(clusterVols5{newcno}(:)==1)=newcno;
                    newcno=newcno+1;
                end
            else
                clusterVols5{newcno}=clusterVols3{i};
                clusterInfo5{newcno}=clusterInfo3{i};
                labels5(clusterVols5{newcno}(:)==1)=newcno;
                newcno=newcno+1;
            end
            
            disp(['cluster num: ' num2str(i) '. new cluster num: ' num2str(newcno)])
            catch
            end
        end
        
        if ~exist([saveDest '\' folderNames{currFolderNames}])
            mkdir([saveDest '\' folderNames{currFolderNames}])
        end
        save([saveDest '\' folderNames{currFolderNames} '\' filePrefix '_autoseg_xystep' num2str(xystep) '.mat'],'clusterVols','clusterInfo','clusterVols2','clusterInfo2','clusterVols3','clusterInfo3','clusterVols4','clusterInfo4','clusterVols5','clusterInfo5','n','b','m','labels','labels2','labels3','labels4','labels5','maskedNnorm','thresh')
        
        disp(['successfully saved data for fly ' num2str(currFlyNum) ' in folder ' folderNames{currFolderNames}])
    end
    cd(startDir)
end

