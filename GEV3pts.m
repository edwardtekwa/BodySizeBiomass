function fval = GEV3pts(params,meanLogS,minLogS,maxLogS,minSizeOffset,maxSizeOffset)
    k=params(1); %shape (~skew)
    s=params(2); %scale (~standard deviation)
    u=params(3); %location (~mean)
    
    pd=makedist('gev','k',k,'sigma',s,'mu',u);
    td=truncate(pd,minLogS-minSizeOffset,maxLogS+maxSizeOffset); %truncate GEV distribution
    modelMinLogS=icdf(td,0.0005);
    modelMaxLogS=icdf(td,0.9995);
    modelMeanLogS=median(pd); %modelled size mode
    
    %compare modelled to observed
    F(1)=modelMinLogS-minLogS;
    F(2)=modelMaxLogS-maxLogS;
    F(3)=modelMeanLogS-meanLogS;
    
    fval=sum(F.^2); %sum of squares
end