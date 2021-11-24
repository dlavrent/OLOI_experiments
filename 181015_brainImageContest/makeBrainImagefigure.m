clear all
close all

f1='clusterResponses_Volumes_181012fly6Left.mat';
f2='clusterResponses_Volumes_181012fly6right.mat';

glomstonotkeep1=[4 17 24 29 30 35 40];
glomstonotkeep2=[20 30];
neworder1=[28 22 2 11 1 12 30 13 17 3 20 33 23 32 26 9 25 29 19 14 27];
neworder2=[13 19 9 11 22 14 23 15 21 5 27 28 10 6 12 1 25 26 4 8 3];

load(f1)
pgloms=1:length(clusterVolU);
tokeep=setxor(pgloms,glomstonotkeep1);

clusterVolU=clusterVolU(tokeep);
clusterInfoU=clusterInfoU(tokeep);

pgloms=1:length(clusterVolU);
no1=[neworder1];% setxor(pgloms,neworder1)];
no1(19)=[];
clusterVolU1=clusterVolU(no1);
clusterInfoU1=clusterInfoU(no1);

%clusterVolU(neworder1)=[];
%clusterInfoU(neworder1)=[];

showClusters(clusterVolU1,clusterInfoU1)
set(gcf, 'color', 'none');
 set(gca, 'color', 'none');
axis off

load(f2)
pgloms=1:length(clusterVolU);
tokeep=setxor(pgloms,glomstonotkeep2);

clusterVolU=clusterVolU(tokeep);
clusterInfoU=clusterInfoU(tokeep);

pgloms=1:length(clusterVolU);
no2=[neworder2];% setxor(pgloms,neworder2)];
no2(19)=[];
clusterVolU2=clusterVolU(no2);
clusterInfoU2=clusterInfoU(no2);

showClusters(clusterVolU2,clusterInfoU2)
set(gcf, 'color', 'none');
 set(gca, 'color', 'none');
 axis off
