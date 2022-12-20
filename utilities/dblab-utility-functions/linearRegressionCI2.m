function out=linearRegressionCI2(X,Y,fitOrder,figureBool,xMin,xMax)

bsReps=2000;
xPoints=100;
%fitOrder=3;

%figureBool=1;
lineColor=[0 0 0];
ciColor=[0.8 0.8 0.8];

linFits=zeros(bsReps,xPoints);
nPoints=length(X);
%xMin=min(X);
%xMax=max(X);
%xVals=linspace(xMin-(xMax-xMin)*0.05,xMax+(xMax-xMin)*0.05,xPoints);
xVals=linspace(xMin,xMax,xPoints);

h=waitbar(0,'bootstrapping progress');
for i=1:bsReps
   
    which=ceil(rand(nPoints,1)*nPoints);
    
    Xtemp=X(which);
    Ytemp=Y(which);
    
    p=polyfit(Xtemp,Ytemp,fitOrder);
    
    linFits(i,:)=polyval(p,xVals);
    
    waitbar(i/bsReps,h);
end

close(h);

out.fits=linFits;
out.pOverall=polyfit(X,Y,fitOrder);
out.xVals=xVals;

if figureBool==1
   
    figure;
    hold on;
    
    areaBar(out.xVals,polyval(out.pOverall,out.xVals),2*std(out.fits),lineColor,ciColor);
    scatter(X,Y,320,'k.');
end

