%% 1. load images
clear all
close all

filePrefix='fly1_';
nc82name=[filePrefix 'nc82.tif'];
brpname=[filePrefix 'brpshort.tif'];
mcd8name=[filePrefix 'mcd8.tif'];

% get image info
nc82info=imfinfo(nc82name);

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



%% 2. median filter and normalize images

zs=1:size(nc82im,3);
tic

n=nc82im(:,:,zs);
b=brpim(:,:,zs);
m=mcd8im(:,:,zs);

xystep=13;
xystart=(xystep-1)/2;
n=medfilt3(n, [xystep xystep 1]); n=n(xystart:xystep:end,xystart:xystep:end,:);
disp('median filtering step 1 complete')
b=medfilt3(b, [xystep xystep 1]); b=b(xystart:xystep:end,xystart:xystep:end,:);
disp('median filtering step 2 complete')
m=medfilt3(m, [xystep xystep 1]); m=m(xystart:xystep:end,xystart:xystep:end,:);
disp('median filtering step 3 complete')


gs=30;
x=[-10:10]; % depends on xystep (smaller range for larger xystep)
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

%% 3.  mask by expression pattern of driver lines

driver=m.*b;
dm=driver;
dm(dm<prctile(dm(:),85))=0;
dm(dm>0)=1;
dm=logical(dm);


% 4.threshold and watershed
maskedNnorm=nNorm.*dm.*mNorm;
thresh=prctile(maskedNnorm(:),90);
asdf=maskedNnorm(:);
thresh=prctile(asdf(asdf>0),50);

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

%% 5. check watershed
close all
figure;

for i=1:size(nNorm,3)
    subplot(1,2,1)
    imagesc(maskedNnorm(:,:,i),[0 50])
    subplot(1,2,2)
    imagesc(labels(:,:,i))
    
    drawnow
    pause(0.1)
end

% %% watershed hyperparameter search
% n=nc82im(:,:,zs);
% figure
% imagesc(mean(n,3))
% pmask=roipoly();
% pmask=double(pmask);
% close all
% 
% 
% xystep=5:2:15;
% minSuppressionThreshold=[1:0.5:3];
% thresh=[10:5:40];
% convLength=[5:10:25];
% 
% nclusters=zeros(length(xystep),length(minSuppressionThreshold),length(thresh),length(convLength));
% 
% for hh=1:length(xystep)
%     n=nc82im(:,:,zs);
%     
%     n=medfilt3(n, [xystep(hh) xystep(hh) 3]); n=n(xystep(hh):xystep(hh):end,xystep(hh):xystep(hh):end,:);
%     disp('median filtering complete')
%     
%     lmask=pmask(xystep(hh):xystep(hh):end,xystep(hh):xystep(hh):end);
%     
%     for ii=1:length(minSuppressionThreshold)
%         for jj=1:length(thresh)
%             for kk=1:length(convLength)
%                 
%                 gs=30;
%                 x=[-convLength(kk):convLength(kk)]; % depends on xystep (smaller range for larger xystep)
%                 y=x;
%                 z=-1:1;
%                 [xx yy zz]=meshgrid(x,y,z);
%                 gau=1/(length(x)*gs^2)*exp(-(xx.^2+yy.^2 +zz.^2)/(2*gs^2));
%                 
%                 nNorm=n./(convn(n,gau,'same'));
%                 im=nNorm<thresh(jj);
%                 im=im.*lmask;
%                 imb=bwdist(im);
%                 imb=imb.*lmask;
%                 
%                 % Watershed transform
%                 preL=max(imb(:))-imb;
%                 preL2=imhmin(preL,minSuppressionThreshold(ii));
%                 L = watershed(preL2);
%                 labels=bwlabeln(L&~im);
%                 labels=labels.*lmask;
%                 
%                 areathresh=100;
%                 areathreshhigh=600000;
%                 
%                 lbls=regionprops(labels,'all');
%                 todelete=[];
%                 clusterVols2=cell(1,length(lbls));
%                 clusterInfo2=cell(1,length(lbls));
%                 bb=zeros(length(lbls),7);
%                 for j=1:length(lbls)
%                     if lbls(j).Area>areathresh && lbls(j).Area<areathreshhigh
%                         clusterVols2{j}=(labels==j);
%                         clusterInfo2{j}=lbls(j);
%                         
%                         bb(j,:)=[clusterInfo2{j}.BoundingBox clusterInfo2{j}.FilledArea];
%                     else
%                         todelete=[todelete j];
%                     end
%                     
%                 end
%                 clusterVols2(todelete)=[];
%                 clusterInfo2(todelete)=[];
%                 bb(todelete,:)=[];
%                 bb=(bb-mean(bb))./std(bb);
%                 
%                 [coeff, score, latent, tsquared, explained] = pca(bb);
%                 try
%                     todelete2=find(abs(score(:,1))>(mean(abs(score(:,1)))+2*std(abs(score(:,1)))));
%                     
%                     clusterVols2(todelete2)=[];
%                     clusterInfo2(todelete2)=[];
%                     bb(todelete2,:)=[];
%                     
%                     
%                     nclusters(hh,ii,jj,kk)=length(clusterVols2);
%                 catch
%                     nclusters(hh,ii,jj,kk)=NaN;
%                 end
%                 disp([num2str(hh) ' ' num2str(ii) ' ' num2str(jj) ' ' num2str(kk)])
%             end
%         end
%     end
% end
% disp('done')




%% 6. save initial watershed results

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


%% 7. apply 3d convolution plus watershed to split merged gloms
splitInfo=cell(1,length(clusterVols));
labels2=zeros(size(clusterVols{1}));
clear clusterVols2 clusterInfo2
newcno=1;
for i=1:length(clusterVols)
    ballr=round(nthroot(clusterInfo{i}.Area,3)/2);
   tempBall= strel('sphere',ballr);
   tempconv=convn(clusterVols{i},tempBall.Neighborhood,'same');
   
   
   tcb=tempconv<prctile(tempconv(tempconv~=0),50);
   tcb2=bwdist(tcb);
   
   % Watershed transform
   minSuppressionThreshold=0.15;
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
end
% 



% 
% figure
% 
% for zslice=1:length(zs)
%     
%     subplot(1,2,1)
%    % imagesc(labels(:,:,zslice))
%      imagesc(maskedNnorm(:,:,zslice),[0 1000])
%     subplot(1,2,2)
%     imagesc(labels2(:,:,zslice),[0 200])
%     drawnow
%     pause(0.5)
% end
% % 
% 
% 
% 
% 
% figure
% view(3);
% axis tight
% camlight
% lighting gouraud
% 
% hold on
% 
% rng('shuffle')
% mycmap=hsv(length(clusterVols));
% mycmap=mycmap(randperm(length(mycmap)),:);
% 
% randv=randperm(length(clusterVols));
% 
% for i=randv
%     p2=patch(isosurface(clusterVols{i}),'FaceColor',mycmap(i,:),'EdgeColor','none','FaceAlpha',0.8);
%     isonormals(clusterVols{i},p2)
%     %text(clusterInfo2{i}.Centroid(1),clusterInfo2{i}.Centroid(2),clusterInfo2{i}.Centroid(3),num2str(i),'Color',[0 0 0],'FontSize',25)
%     drawnow
% end
% title('original watershed')
% 
% figure
% view(3);
% axis tight
% camlight
% lighting gouraud
% 
% hold on
% 
% rng('default')
% mycmap=hsv(length(clusterVols2));
% mycmap=mycmap(randperm(length(mycmap)),:);
% 
% randv=randperm(length(clusterVols2));
% 
% for i=randv
%     p2=patch(isosurface(clusterVols2{i}),'FaceColor',mycmap(i,:),'EdgeColor','none','FaceAlpha',0.8);
%     isonormals(clusterVols2{i},p2)
%     %text(clusterInfo2{i}.Centroid(1),clusterInfo2{i}.Centroid(2),clusterInfo2{i}.Centroid(3),num2str(i),'Color',[0 0 0],'FontSize',25)
%     drawnow
% end
% title('watershed+3dconv+watershed')
% disp('done plotting')
%% 7. perform pca on pixels

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

% unwrap images for kmeans
%data = [n(:)'; b(:)'; m(:)'; xxx(:)'; yyy(:)'; zzz(:)']';
data = [bm(:)'; mm(:)'; transpose(bm(:).*mm(:)); xxx(:)'; yyy(:)'; zzz(:)']';


data=(data-nanmean(data))./nanstd(data);

interactionweight=0.75;
xyweight=1;
zweight=xyweight;

data(:,3)=interactionweight*data(:,3);
data(:,(end-2):(end-1))=xyweight*data(:,(end-2):(end-1));
data(:,end)=zweight*data(:,end);


disp('finished preprocessing')

% run kmeans and pca

% [IDX, C] = kmeans(gpuArray(data),100,'Display','iter','MaxIter',500);
% 
% kout=reshape(gather(IDX),size(n));
% cout=gather(C);
% disp('done')


[coeff, score, latent, tsquared, explained] = pca(data);

disp('finished pca')
%Y = tsne(data,'Perplexity',20,'Verbose',2);

% %%  plot pcs and raw data for pixels within each cluster IDed by watershed
% numclustersfound=length(unique(labelsEdited))-1;
% uniqueclusters=unique(labelsEdited);
% figure
% myc=hsv(numclustersfound);
% randv=randperm(numclustersfound);
% for i= 1:numclustersfound
%     temp=find(labelsEdited==uniqueclusters(i+1));
%     plot3(score(temp,1),score(temp,2),score(temp,3),'o','Color',myc(randv(i),:))
%     hold on
% end
% xlabel('pc 1')
% ylabel('pc 2')
% zlabel('pc 3')
% 
% 
% figure
% 
% for i= 1:numclustersfound
%     temp=find(labelsEdited==uniqueclusters(i+1));
%     plot(data(temp,1),data(temp,2),'o','Color',myc(randv(i),:))
%     hold on
% end

%% 8. split clusters based on pixel similarity
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

figure

for zslice=1:size(labels,3)
    
    subplot(2,2,1)
    imagesc(labels(:,:,zslice))
    
    subplot(2,2,2)
    imagesc(labels2(:,:,zslice))
    
        subplot(2,2,3)
    imagesc(labels3(:,:,zslice))
    
        subplot(2,2,4)
    imagesc(labels4(:,:,zslice))
    drawnow
    pause(0.1)
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
%% use 3d convolution plus watershed
splitInfo=cell(1,length(clusterVols3));
labels4=zeros(size(clusterVols3{1}));
clear clusterVols4 clusterInfo4
newcno=1;
for i=1:length(clusterVols3)
    ballr=round(nthroot(clusterInfo3{i}.Area,3)/2);
   tempBall= strel('sphere',ballr);
   tempconv=convn(clusterVols3{i},tempBall.Neighborhood,'same');
   
   
   tcb=tempconv<prctile(tempconv(tempconv~=0),50);
   tcb2=bwdist(tcb);
   
   % Watershed transform
   minSuppressionThreshold=0.15;
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
end
% 



% 
% figure
% 
% for zslice=1:length(zs)
%     
%     subplot(1,3,1)
%     imagesc(labels(:,:,zslice))
%     
%     subplot(1,3,2)
%     imagesc(newlabels(:,:,zslice),[0 100])
%     subplot(1,3,3)
%     imagesc(finalLabels(:,:,zslice),[0 200])
%     drawnow
%     pause(0.1)
% end
% % 
% 


% 
% figure
% view(3);
% axis tight
% camlight
% lighting gouraud
% 
% hold on
% 
% rng('shuffle')
% mycmap=hsv(length(clusterVols));
% mycmap=mycmap(randperm(length(mycmap)),:);
% 
% randv=randperm(length(clusterVols));
% 
% for i=randv
%     p2=patch(isosurface(clusterVols{i}),'FaceColor',mycmap(i,:),'EdgeColor','none','FaceAlpha',0.8);
%     isonormals(clusterVols{i},p2)
%     %text(clusterInfo2{i}.Centroid(1),clusterInfo2{i}.Centroid(2),clusterInfo2{i}.Centroid(3),num2str(i),'Color',[0 0 0],'FontSize',25)
%     drawnow
% end
% title('original watershed')
% 
% figure
% view(3);
% axis tight
% camlight
% lighting gouraud
% 
% hold on
% 
% rng('shuffle')
% mycmap=hsv(length(clusterVols3));
% mycmap=mycmap(randperm(length(mycmap)),:);
% 
% randv=randperm(length(clusterVols3));
% 
% for i=randv
%     p2=patch(isosurface(clusterVols3{i}),'FaceColor',mycmap(i,:),'EdgeColor','none','FaceAlpha',0.8);
%     isonormals(clusterVols3{i},p2)
%     %text(clusterInfo2{i}.Centroid(1),clusterInfo2{i}.Centroid(2),clusterInfo2{i}.Centroid(3),num2str(i),'Color',[0 0 0],'FontSize',25)
%     drawnow
% end
% title('watershed+split')
% 
% figure
% view(3);
% axis tight
% camlight
% lighting gouraud
% 
% hold on
% 
% rng('default')
% mycmap=hsv(length(clusterVols4));
% mycmap=mycmap(randperm(length(mycmap)),:);
% 
% randv=randperm(length(clusterVols4));
% 
% for i=randv
%     p2=patch(isosurface(clusterVols4{i}),'FaceColor',mycmap(i,:),'EdgeColor','none','FaceAlpha',0.8);
%     isonormals(clusterVols4{i},p2)
%     %text(clusterInfo2{i}.Centroid(1),clusterInfo2{i}.Centroid(2),clusterInfo2{i}.Centroid(3),num2str(i),'Color',[0 0 0],'FontSize',25)
%     drawnow
% end
% title('watershed+split+watershed')
% disp('done plotting')

