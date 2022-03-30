function [X,fval,CDFpts]=fitGEV(Sizes,minSizeOffset,maxSizeOffset)
%fit skew normal distribution pdf to min, max, and mean and return
%location, scale, and shape
%Edward Tekwa Oct 11, 2019

%convert to log scale
meanLogS=log10(Sizes(2));
minLogS=log10(Sizes(1));
maxLogS=log10(Sizes(3));

%fit distribution to observed min, mean, and max body sizes
X0=[0,((maxLogS-minLogS)/5),meanLogS]; %initial guesses for k (shape), sigma (scale), mu (location)
Xmin=[-Inf,eps,-Inf];
Xmax=[Inf,Inf,Inf];
options = optimset('MaxFunEvals',500,'MaxIter',500,'Display','off','TolFun',1e-20,'TolX',1e-20);
[X,fval,exitflag,output] = fminsearchbnd(@(params) GEV3pts(params,meanLogS,minLogS,maxLogS,minSizeOffset,maxSizeOffset),X0,Xmin,Xmax,options);

%reference CDFs, should get close to 0.005 (min), something (mean), and
%0.995 (max)
CDF005=gevcdf(minLogS,X(1),X(2),X(3));
CDF500=gevcdf(meanLogS,X(1),X(2),X(3));
CDF995=gevcdf(maxLogS,X(1),X(2),X(3));

X;
CDFpts=[CDF005 CDF500 CDF995];