%bodySizeGEVPowerLawsBootstrapped.m
%Edward Tekwa March 2, 2020

%run bodySizeGEVPlot_indivGroup.m first!

%Figure 4 or 5
BA=1; %specify regression type: (1) biomass - Figure 5, or (2) abundance - Figure 4
numOrders=3; %three habitat realm groups for regression analyses: all, terrestrial, marine (rows in figure)
NotTypeCodes={[];[3 4 5];[1 2 5]}; %exclusion codes: 1. terrestrial producer, 2. terrestrial consumer, 3. marine producer, 4. marine consumer, 5. subterranean
LogOffset=-log10(10^-5); %add this to the log10(biomass) to ensure positive values for regression (0.0002 10^-5, 10^-10)
lowerPerc=2.5; %2.5 : lower bound of AICc or BIC
upperPerc=97.5; %97.5: upper bound of AICc or BIC

fitObjects1=nan(2,5,numOrders);
fitObjects2=nan(2,5,numOrders);
fitObjects3=nan(2,5,numOrders);
allLogBS=[];
allLogBiomass=[];
allClass=[];

LinEsts=nan(2,2,numOrders);
Gauss1Ests=nan(2,3,numOrders);
Gauss2Ests=nan(2,6,numOrders);
Gauss3Ests=nan(2,9,numOrders);
Gauss4Ests=nan(2,12,numOrders);
ModelStats=nan(5,11,numOrders);
modelsSelected=zeros(1,numOrders);


scrsz = get(0,'ScreenSize');
set(0, 'DefaultAxesFontSize', 16)
set(0,'defaultaxeslinewidth',1)
set(0,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
%figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2 scrsz(4)]);
figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2 scrsz(4)/2.5]);

for k=1:numOrders
    %select data subset to exclude from plot:
    NotTypes=NotTypeCodes{k};
    typeidx=[];
    if BA==1
        plot_orderedBootSortedBiomassMatrix=All_orderedBootSortedBiomassMatrix;
    elseif BA==2
        plot_orderedBootSortedBiomassMatrix=All_orderedBootSortedBiomassMatrix./(10.^BSbins);
    end
    plotBSGroups=sortedBSGroups(:,1:3);
    if ~isempty(NotTypes)
        for j=1:length(NotTypes)
            typeidx = [typeidx; find(sortedBSGroups(:,5)==NotTypes(j))];
        end
        plot_orderedBootSortedBiomassMatrix(:,:,typeidx)=0;
        plotBSGroups(typeidx,:)=NaN;
    end
    BootSortedBiomassMatrix=sum(plot_orderedBootSortedBiomassMatrix,3);
    orderedBootSortedBiomassMatrix=sort(BootSortedBiomassMatrix); %this generates the empirical confidence bounds
    
    %begin bootstrapped linear and Gaussian mixture regressions (entire size range is a single correlated block,
    %pseudoreplicates are drawn from cumulative biomass range sets)
    %set up x-variable for regressions
    LinVals=nan(numBoots,2); %linear regression coefficients from all bootstraps
    Gauss1Vals=nan(numBoots,3); %Gaussian (1) regression coefficients from all bootstraps
    Gauss2Vals=nan(numBoots,6); %Gaussian (2) regression coefficients from all bootstraps
    Gauss3Vals=nan(numBoots,9); %Gaussian (3) regression coefficients from all bootstraps
    Gauss4Vals=nan(numBoots,12); %Gaussian (4) regression coefficients from all bootstraps
    LinPred=nan(numBoots,size(BootSortedBiomassMatrix,2)); %linear predictions
    Gauss1Pred=nan(numBoots,size(BootSortedBiomassMatrix,2)); %Gaussian (1) predictions
    Gauss2Pred=nan(numBoots,size(BootSortedBiomassMatrix,2)); %Gaussian (2) predictions
    Gauss3Pred=nan(numBoots,size(BootSortedBiomassMatrix,2)); %Gaussian (3) predictions
    Gauss4Pred=nan(numBoots,size(BootSortedBiomassMatrix,2)); %Gaussian (4) predictions
    LinR2=nan(numBoots,1);
    LinSSE=nan(numBoots,1);
    relAICcLin=nan(numBoots,1);
    relBICLin=nan(numBoots,1);
    Gauss1R2=nan(numBoots,1);
    Gauss2R2=nan(numBoots,1);
    Gauss3R2=nan(numBoots,1);
    Gauss4R2=nan(numBoots,1);
    relAICcGauss1=nan(numBoots,1);
    relAICcGauss2=nan(numBoots,1);
    relAICcGauss3=nan(numBoots,1);
    relAICcGauss4=nan(numBoots,1);
    relBICGauss1=nan(numBoots,1);
    relBICGauss2=nan(numBoots,1);
    relBICGauss3=nan(numBoots,1);
    relBICGauss4=nan(numBoots,1);
    relAICcGauss1Diff=nan(numBoots,1);
    relAICcGauss2Diff=nan(numBoots,1);
    relAICcGauss3Diff=nan(numBoots,1);
    relAICcGauss4Diff=nan(numBoots,1);
    relBICGauss1Diff=nan(numBoots,1);
    relBICGauss2Diff=nan(numBoots,1);
    relBICGauss3Diff=nan(numBoots,1);
    relBICGauss4Diff=nan(numBoots,1);
    
    TimesModelSelected=zeros(1,5);
    
    for boot=1:numBoots
        logBiomass=log10(BootSortedBiomassMatrix(boot,:)); %set up biomass data for many-point regression (1):
        logBiomass(logBiomass==-Inf)=NaN; %take out 0 biomass points
        logBiomassFull=logBiomass;
        logBiomassFull(logBiomassFull<-LogOffset+1)=NaN; %take out biomass <0.0002 GtC
        RealBiomassPowerLawBoot=fitlm(BSbins,logBiomassFull); %with x-value being all pseudo datapoints across log body size
        LinVals(boot,:)=RealBiomassPowerLawBoot.Coefficients.Estimate'; %record coefficients
        LinPred(boot,:)=predict(RealBiomassPowerLawBoot, BSbins'); %record model predictions
        LinR2(boot,1)=RealBiomassPowerLawBoot.Rsquared.Ordinary;
        LinSSE(boot,1)=RealBiomassPowerLawBoot.SSE;
        relAICcLin(boot,1)=2*3*(1+(3+1)/(length(BSbins)-3-1))+length(BSbins)*log(RealBiomassPowerLawBoot.SSE);
        relBICLin(boot,1)=3*log(length(BSbins))+length(BSbins)*log(RealBiomassPowerLawBoot.SSE);
        
        TransformedYBoot=logBiomassFull(~isnan(logBiomassFull))+LogOffset; %add LogOffset orders of magnitude to avoid negative values
        BSbinsSubset=BSbins(~isnan(logBiomassFull));
        [Gaussian1Boot_fitobj, Gaussian1Boot_gof] = fit(BSbinsSubset',TransformedYBoot','gauss1','Lower',[0 -Inf 0]); %fit 1 Gaussian on linear biomass scale
        [Gaussian2Boot_fitobj, Gaussian2Boot_gof] = fit(BSbinsSubset',TransformedYBoot','gauss2','Lower',[0 -Inf 0 0 -Inf 0]); %fit 2 Gaussian mixtures on linear biomass scale
        [Gaussian3Boot_fitobj, Gaussian3Boot_gof] = fit(BSbinsSubset',TransformedYBoot','gauss3','Lower',[0 -Inf 0 0 -Inf 0 0 -Inf 0]); %fit 3 Gaussian mixtures on linear biomass scale
        [Gaussian4Boot_fitobj, Gaussian4Boot_gof] = fit(BSbinsSubset',TransformedYBoot','gauss4','Lower',[0 -Inf 0 0 -Inf 0 0 -Inf 0 0 -Inf 0]); %fit 4 Gaussian mixtures on linear biomass scale
        Gauss1Vals(boot,:)=coeffvalues(Gaussian1Boot_fitobj); %record coefficients
        Gauss2Vals(boot,:)=coeffvalues(Gaussian2Boot_fitobj); %record coefficients
        Gauss3Vals(boot,:)=coeffvalues(Gaussian3Boot_fitobj); %record coefficients
        Gauss4Vals(boot,:)=coeffvalues(Gaussian4Boot_fitobj); %record coefficients
        Gauss1Pred(boot,:) = Gaussian1Boot_fitobj(BSbins'); %record model predictions
        Gauss2Pred(boot,:) = Gaussian2Boot_fitobj(BSbins'); %record model predictions
        Gauss3Pred(boot,:) = Gaussian3Boot_fitobj(BSbins'); %record model predictions
        Gauss4Pred(boot,:) = Gaussian4Boot_fitobj(BSbins'); %record model predictions
        Gauss1R2(boot,1)=Gaussian1Boot_gof.rsquare;
        Gauss1SSE(boot,1)=Gaussian1Boot_gof.sse;
        Gauss2R2(boot,1)=Gaussian2Boot_gof.rsquare;
        Gauss2SSE(boot,1)=Gaussian2Boot_gof.sse;
        Gauss3R2(boot,1)=Gaussian3Boot_gof.rsquare;
        Gauss3SSE(boot,1)=Gaussian3Boot_gof.sse;
        Gauss4R2(boot,1)=Gaussian4Boot_gof.rsquare;
        Gauss4SSE(boot,1)=Gaussian4Boot_gof.sse;
        relAICcGauss1(boot,1)=2*4*(1+(4+1)/(length(BSbins)-4-1))+length(BSbins)*log(Gaussian1Boot_gof.sse);
        relAICcGauss2(boot,1)=2*7*(1+(7+1)/(length(BSbins)-7-1))+length(BSbins)*log(Gaussian2Boot_gof.sse);
        relAICcGauss3(boot,1)=2*10*(1+(10+1)/(length(BSbins)-10-1))+length(BSbins)*log(Gaussian3Boot_gof.sse);
        relAICcGauss4(boot,1)=2*13*(1+(13+1)/(length(BSbins)-13-1))+length(BSbins)*log(Gaussian4Boot_gof.sse);
        relBICGauss1(boot,1)=4*log(length(BSbins))+length(BSbins)*log(Gaussian1Boot_gof.sse);
        relBICGauss2(boot,1)=7*log(length(BSbins))+length(BSbins)*log(Gaussian2Boot_gof.sse);
        relBICGauss3(boot,1)=10*log(length(BSbins))+length(BSbins)*log(Gaussian3Boot_gof.sse);
        relBICGauss4(boot,1)=13*log(length(BSbins))+length(BSbins)*log(Gaussian4Boot_gof.sse);
        
        %get pairwise AIC BIC differences compare to next simplest model):
        relAICcGauss1Diff(boot,1)=relAICcGauss1(boot,1)-relAICcLin(boot,1);
        relAICcGauss2Diff(boot,1)=relAICcGauss2(boot,1)-relAICcGauss1(boot,1);
        relAICcGauss3Diff(boot,1)=relAICcGauss3(boot,1)-relAICcGauss2(boot,1);
        relAICcGauss4Diff(boot,1)=relAICcGauss4(boot,1)-relAICcGauss3(boot,1);
        
        relBICGauss1Diff(boot,1)=relBICGauss1(boot,1)-relBICLin(boot,1);
        relBICGauss2Diff(boot,1)=relBICGauss2(boot,1)-relBICGauss1(boot,1);
        relBICGauss3Diff(boot,1)=relBICGauss3(boot,1)-relBICGauss2(boot,1);
        relBICGauss4Diff(boot,1)=relBICGauss4(boot,1)-relBICGauss3(boot,1);
        
        [temp,bootModelSelected]=min([relAICcLin(boot,1) relAICcGauss1(boot,1) relAICcGauss2(boot,1) relAICcGauss3(boot,1) relAICcGauss4(boot,1)]); %select Gaussian mixture model within bootstrap
        TimesModelSelected(bootModelSelected)=TimesModelSelected(bootModelSelected)+1; %add to tally of which model is selected

    end
    
    %get frequency of a model being selected based on minAICc or minBIC 
    
    meanLinearPred=mean(LinPred); %get mean regression curve
    lowerLinearPred=prctile(LinPred,2.5); %get lower bound
    upperLinearPred=prctile(LinPred,97.5); %get upper bound
    meanLinVals=mean(LinVals); %mean coefficients
    stdLinVals=std(LinVals); %standard deviation of coefficients
    meanLinR2=mean(LinR2); %mean R2
    stdLinR2=std(LinR2); %standard deviation R2
    mean_relAICcLin=mean(relAICcLin); %mean AICc
    std_relAICcLin=std(relAICcLin);
    lower_relAICcLin=prctile(relAICcLin,lowerPerc); %AICc CI
    upper_relAICcLin=prctile(relAICcLin,upperPerc); %AICc CI
    mean_relBICLin=mean(relBICLin); %mean AICc
    lower_relBICLin=prctile(relBICLin,lowerPerc); %AICc CI
    upper_relBICLin=prctile(relBICLin,upperPerc); %AICc CI
    LinEsts(:,:,k)=[meanLinVals;stdLinVals];
    pAICcLin=0; %probability next simplest model is better
    pBICLin=0;
    ModelStats(1,:,k)=[meanLinR2 stdLinR2 mean_relAICcLin lower_relAICcLin upper_relAICcLin pAICcLin TimesModelSelected(1)/numBoots mean_relBICLin lower_relBICLin upper_relBICLin pBICLin];
    

    meanGauss1Pred=mean(Gauss1Pred); %get mean regression curve
    lowerGauss1Pred=prctile(Gauss1Pred,2.5); %get lower bound
    upperGauss1Pred=prctile(Gauss1Pred,97.5); %get upper bound
    meanGauss1Vals=mean(Gauss1Vals); %mean coefficients
    stdGauss1Vals=std(Gauss1Vals); %standard deviation of coefficients
    meanGauss1R2=mean(Gauss1R2); %mean R2
    stdGauss1R2=std(Gauss1R2); %standard deviation R2
    mean_relAICcGauss1=mean(relAICcGauss1); %mean AICc
    std_relAICcGauss1=std(relAICcGauss1); 
    lower_relAICcGauss1=prctile(relAICcGauss1,lowerPerc); %AICc CI
    upper_relAICcGauss1=prctile(relAICcGauss1,upperPerc); %AICc CI
    mean_relBICGauss1=mean(relBICGauss1); %mean AICc
    lower_relBICGauss1=prctile(relBICGauss1,lowerPerc); %AICc CI
    upper_relBICGauss1=prctile(relBICGauss1,upperPerc); %AICc CI
    Gauss1Ests(:,:,k)=[meanGauss1Vals;stdGauss1Vals];
    pAICcGauss1=sum(relAICcGauss1Diff>0)/numBoots; %probability next simplest model is better
    pBICGauss1=sum(relBICGauss1Diff>0)/numBoots;
    ModelStats(2,:,k)=[meanGauss1R2 stdGauss1R2 mean_relAICcGauss1 lower_relAICcGauss1 upper_relAICcGauss1 pAICcGauss1 TimesModelSelected(2)/numBoots mean_relBICGauss1 lower_relBICGauss1 upper_relBICGauss1 pBICGauss1];
    
    
    meanGauss2Pred=mean(Gauss2Pred); %get mean regression curve
    lowerGauss2Pred=prctile(Gauss2Pred,2.5); %get lower bound
    upperGauss2Pred=prctile(Gauss2Pred,97.5); %get upper bound
    meanGauss2Vals=mean(Gauss2Vals); %mean coefficients
    stdGauss2Vals=std(Gauss2Vals); %standard deviation of coefficients
    meanGauss2R2=mean(Gauss2R2); %mean R2
    stdGauss2R2=std(Gauss2R2); %standard deviation R2
    mean_relAICcGauss2=mean(relAICcGauss2); %mean AICc
    std_relAICcGauss2=std(relAICcGauss2); 
    lower_relAICcGauss2=prctile(relAICcGauss2,lowerPerc); %AICc CI
    upper_relAICcGauss2=prctile(relAICcGauss2,upperPerc); %AICc CI
    mean_relBICGauss2=mean(relBICGauss2); %mean AICc
    lower_relBICGauss2=prctile(relBICGauss2,lowerPerc); %AICc CI
    upper_relBICGauss2=prctile(relBICGauss2,upperPerc); %AICc CI
    Gauss2Ests(:,:,k)=[meanGauss2Vals;stdGauss2Vals];
    pAICcGauss2=sum(relAICcGauss2Diff>0)/numBoots; %probability next simplest model is better
    pBICGauss2=sum(relBICGauss2Diff>0)/numBoots;
    ModelStats(3,:,k)=[meanGauss2R2 stdGauss2R2 mean_relAICcGauss2 lower_relAICcGauss2 upper_relAICcGauss2 pAICcGauss2 TimesModelSelected(3)/numBoots mean_relBICGauss2 lower_relBICGauss2 upper_relBICGauss2 pBICGauss2];
    
    
    meanGauss3Pred=mean(Gauss3Pred); %get mean regression curve
    lowerGauss3Pred=prctile(Gauss3Pred,2.5); %get lower bound
    upperGauss3Pred=prctile(Gauss3Pred,97.5); %get upper bound
    meanGauss3Vals=mean(Gauss3Vals); %mean coefficients
    stdGauss3Vals=std(Gauss3Vals); %standard deviation of coefficients
    meanGauss3R2=mean(Gauss3R2); %mean R2
    stdGauss3R2=std(Gauss3R2); %standard deviation R2
    mean_relAICcGauss3=mean(relAICcGauss3); %mean AICc
    std_relAICcGauss3=std(relAICcGauss3); 
    lower_relAICcGauss3=prctile(relAICcGauss3,lowerPerc); %AICc CI
    upper_relAICcGauss3=prctile(relAICcGauss3,upperPerc); %AICc CI
    mean_relBICGauss3=mean(relBICGauss3); %mean AICc
    lower_relBICGauss3=prctile(relBICGauss3,lowerPerc); %AICc CI
    upper_relBICGauss3=prctile(relBICGauss3,upperPerc); %AICc CI
    Gauss3Ests(:,:,k)=[meanGauss3Vals;stdGauss3Vals];
    pAICcGauss3=sum(relAICcGauss3Diff>0)/numBoots; %probability next simplest model is better
    pBICGauss3=sum(relBICGauss3Diff>0)/numBoots;
    ModelStats(4,:,k)=[meanGauss3R2 stdGauss3R2 mean_relAICcGauss3 lower_relAICcGauss3 upper_relAICcGauss3 pAICcGauss3 TimesModelSelected(4)/numBoots mean_relBICGauss3 lower_relBICGauss3 upper_relBICGauss3 pBICGauss3];
    
    
    meanGauss4Pred=mean(Gauss4Pred); %get mean regression curve
    lowerGauss4Pred=prctile(Gauss4Pred,2.5); %get lower bound
    upperGauss4Pred=prctile(Gauss4Pred,97.5); %get upper bound
    meanGauss4Vals=mean(Gauss4Vals); %mean coefficients
    stdGauss4Vals=std(Gauss4Vals); %standard deviation of coefficients
    meanGauss4R2=mean(Gauss4R2); %mean R2
    stdGauss4R2=std(Gauss4R2); %standard deviation R2
    mean_relAICcGauss4=mean(relAICcGauss4); %mean AICc
    std_relAICcGauss4=std(relAICcGauss4); 
    lower_relAICcGauss4=prctile(relAICcGauss4,lowerPerc); %AICc CI
    upper_relAICcGauss4=prctile(relAICcGauss4,upperPerc); %AICc CI
    mean_relBICGauss4=mean(relBICGauss4); %mean AICc
    lower_relBICGauss4=prctile(relBICGauss4,lowerPerc); %AICc CI
    upper_relBICGauss4=prctile(relBICGauss4,upperPerc); %AICc CI
    Gauss4Ests(:,:,k)=[meanGauss4Vals;stdGauss4Vals];
    pAICcGauss4=sum(relAICcGauss4Diff>0)/numBoots; %probability next simplest model is better
    pBICGauss4=sum(relBICGauss4Diff>0)/numBoots;
    ModelStats(5,:,k)=[meanGauss4R2 stdGauss4R2 mean_relAICcGauss4 lower_relAICcGauss4 upper_relAICcGauss4 pAICcGauss4 TimesModelSelected(5)/numBoots mean_relBICGauss4 lower_relBICGauss4 upper_relBICGauss4 pBICGauss4];
    
    %Gaussian model selection based on greatest loss in BIC score with an
    %incremental model complexity (with base score being from the linear
    %model)
    %[temp modelSelected]=min(diff([mean_relBICLin mean_relBICGauss1 mean_relBICGauss2 mean_relBICGauss3 mean_relBICGauss4]));
    %[temp modelSelected]=min(diff([mean_relAICcLin mean_relAICcGauss1 mean_relAICcGauss2 mean_relAICcGauss3 mean_relAICcGauss4]));
    
    [temp,modelSelected]=min([mean_relAICcGauss1 mean_relAICcGauss2 mean_relAICcGauss3 mean_relAICcGauss4]); %select Gaussian mixture model
    %[temp,modelSelected]=min([mean_relBICGauss1 mean_relBICGauss2 mean_relBICGauss3 mean_relBICGauss4]); %select Gaussian mixture model based on mean BIC
    
    %     %Gaussian model selection based on significantly lowest BIC score:
    %     if modelSelected==4 %&& upper_relBICGauss4<lower_relBICGauss3 %if Gauss4 is significantly better than Gauss3, choose Gauss4
    %         lowerGaussPred=lowerGauss4Pred;
    %         upperGaussPred=upperGauss4Pred;
    %         GaussR2=meanGauss4R2;
    %         mean_relBICGauss=mean_relBICGauss4;
    %         lower_relBICGauss=lower_relBICGauss4;
    %         upper_relBICGauss=upper_relBICGauss4;
    %     elseif modelSelected>=3 %&& upper_relBICGauss3<lower_relBICGauss2 %if Gauss3 is significantly better than Gauss2, choose Gauss3
    %         modelSelected=3;
    %         lowerGaussPred=lowerGauss3Pred;
    %         upperGaussPred=upperGauss3Pred;
    %         GaussR2=meanGauss3R2;
    %         mean_relBICGauss=mean_relBICGauss3;
    %         lower_relBICGauss=lower_relBICGauss3;
    %         upper_relBICGauss=upper_relBICGauss3;
    %     elseif modelSelected>=2 %&& upper_relBICGauss2<lower_relBICGauss1 %if Gauss2 is significantly better than Gauss1, choose Gauss2
    %         modelSelected=2;
    %         lowerGaussPred=lowerGauss2Pred;
    %         upperGaussPred=upperGauss2Pred;
    %         GaussR2=meanGauss2R2;
    %         mean_relBICGauss=mean_relBICGauss2;
    %         lower_relBICGauss=lower_relBICGauss2;
    %         upper_relBICGauss=upper_relBICGauss2;
    %     else
    %         modelSelected=1;
    %         lowerGaussPred=lowerGauss1Pred;
    %         upperGaussPred=upperGauss1Pred;
    %         GaussR2=meanGauss1R2;
    %         mean_relBICGauss=mean_relBICGauss1;
    %         lower_relBICGauss=lower_relBICGauss1;
    %         upper_relBICGauss=upper_relBICGauss1;
    %     end
    
    %Gaussian model selection based on significantly lowest AICc score:
%     if modelSelected==4 %&& upper_relAICcGauss4<lower_relAICcGauss3 %if Gauss4 is significantly better than Gauss3, choose Gauss4
%         lowerGaussPred=lowerGauss4Pred;
%         upperGaussPred=upperGauss4Pred;
%         GaussR2=meanGauss4R2;
%         mean_relAICcGauss=mean_relAICcGauss4;
%         lower_relAICcGauss=lower_relAICcGauss4;
%         upper_relAICcGauss=upper_relAICcGauss4;
%     elseif modelSelected>=3 %&& upper_relAICcGauss3<lower_relAICcGauss2 %if Gauss3 is significantly better than Gauss2, choose Gauss3
%         modelSelected=3;
%         lowerGaussPred=lowerGauss3Pred;
%         upperGaussPred=upperGauss3Pred;
%         GaussR2=meanGauss3R2;
%         mean_relAICcGauss=mean_relAICcGauss3;
%         lower_relAICcGauss=lower_relAICcGauss3;
%         upper_relAICcGauss=upper_relAICcGauss3;
%     elseif modelSelected>=2 %&& upper_relAICcGauss2<lower_relAICcGauss1 %if Gauss2 is significantly better than Gauss1, choose Gauss2
%         modelSelected=2;
%         lowerGaussPred=lowerGauss2Pred;
%         upperGaussPred=upperGauss2Pred;
%         GaussR2=meanGauss2R2;
%         mean_relAICcGauss=mean_relAICcGauss2;
%         lower_relAICcGauss=lower_relAICcGauss2;
%         upper_relAICcGauss=upper_relAICcGauss2;
%     else
%         modelSelected=1;
%         lowerGaussPred=lowerGauss1Pred;
%         upperGaussPred=upperGauss1Pred;
%         GaussR2=meanGauss1R2;
%         mean_relAICcGauss=mean_relAICcGauss1;
%         lower_relAICcGauss=lower_relAICcGauss1;
%         upper_relAICcGauss=upper_relAICcGauss1;
%     end
    %Gaussian model selection based on lowest AICc score that is lower than next simplest model:
    if modelSelected==4 && pAICcGauss4<=0.05 %if Gauss4 is significantly better than Gauss3, choose Gauss4
        lowerGaussPred=lowerGauss4Pred;
        upperGaussPred=upperGauss4Pred;
        GaussR2=meanGauss4R2;
        mean_relAICcGauss=mean_relAICcGauss4;
        lower_relAICcGauss=lower_relAICcGauss4;
        upper_relAICcGauss=upper_relAICcGauss4;
    elseif modelSelected>=3 && pAICcGauss3<=0.05 %if Gauss3 is significantly better than Gauss2, choose Gauss3
        modelSelected=3;
        lowerGaussPred=lowerGauss3Pred;
        upperGaussPred=upperGauss3Pred;
        GaussR2=meanGauss3R2;
        mean_relAICcGauss=mean_relAICcGauss3;
        lower_relAICcGauss=lower_relAICcGauss3;
        upper_relAICcGauss=upper_relAICcGauss3;
    elseif modelSelected>=2 && pAICcGauss2<=0.05 %if Gauss2 is significantly better than Gauss1, choose Gauss2
        modelSelected=2;
        lowerGaussPred=lowerGauss2Pred;
        upperGaussPred=upperGauss2Pred;
        GaussR2=meanGauss2R2;
        mean_relAICcGauss=mean_relAICcGauss2;
        lower_relAICcGauss=lower_relAICcGauss2;
        upper_relAICcGauss=upper_relAICcGauss2;
    else
        modelSelected=1;
        lowerGaussPred=lowerGauss1Pred;
        upperGaussPred=upperGauss1Pred;
        GaussR2=meanGauss1R2;
        mean_relAICcGauss=mean_relAICcGauss1;
        lower_relAICcGauss=lower_relAICcGauss1;
        upper_relAICcGauss=upper_relAICcGauss1;
    end
    modelsSelected(k)=modelSelected;
    
    %subplot(2,3,k)
    subplot(3,5,(k-1)*5+1)
    %yyaxis left
    hold on
    
    %plot raw spectrum
    %CIcolor=[191 0 191]./255; %purple
    CIcolor=[200 200 200]./255; %grey [220 220 220]./25
    plot(BSbins,log10(prctile(orderedBootSortedBiomassMatrix,2.5)'),'-','Color',CIcolor,'LineWidth',2) %lower bound of raw data
    plot(BSbins,log10(prctile(orderedBootSortedBiomassMatrix,97.5)'),'-','Color',CIcolor,'LineWidth',2) %upper bound of raw data
    x2 = [BSbins, fliplr(BSbins)];
    inBetween = [log10(prctile(orderedBootSortedBiomassMatrix,2.5)), fliplr(log10(prctile(orderedBootSortedBiomassMatrix,97.5)))];
    fill(x2, inBetween,CIcolor,'EdgeColor','none');
    %patch(x2,inBetween);
    
    %plot confidence bounds of linear model
    LinColor=[119 172 48]/255;
    plot(BSbins,lowerLinearPred,'-b','LineWidth',2);
    plot(BSbins,upperLinearPred,'-b','LineWidth',2);

    legend 'off'
    %xlim(BSlims)
    if BA==1
        ylim(log10([0.0002 3000])) %for log body size-log biomass 200 [1E-7 2000]
        text(-13,-2.5,{['\beta=' num2str(meanLinVals(2),2) '\pm' num2str(num2str(stdLinVals(2),2))];['R^2=' num2str(meanLinR2,2) '\pm' num2str(num2str(stdLinR2,2))];['AICc=' num2str(mean_relAICcLin,2) '\pm' num2str(num2str(std_relAICcLin,2))]},'fontsize',16,'Color','b')

    elseif BA==2
        ylim([-13 18]) %for log body size-log abundance
        text(-13,-2.5,{['\alpha=' num2str(meanLinVals(2),2) '\pm' num2str(num2str(stdLinVals(2),2))];['R^2=' num2str(meanLinR2,2) '\pm' num2str(num2str(stdLinR2,2))];['AICc=' num2str(mean_relAICcLin,2) '\pm' num2str(num2str(std_relAICcLin,2))]},'fontsize',16,'Color','b')

    end
    yt=yticks;
    if BA==2
        yyaxis right
        ylim([-13 18]-0.454)
        yticks(yt);
        yticklabels(round(yt-0.454,1)); %normalized biomass
    end
    
    subplot(3,5,(k-1)*5+2)
    hold on
    %ylim(log10([0.0002 3000])) %for log body size-log biomass 200 [1E-7 2000]
    plot(BSbins,log10(prctile(orderedBootSortedBiomassMatrix,2.5)'),'-','Color',CIcolor,'LineWidth',2) %lower bound of raw data
    plot(BSbins,log10(prctile(orderedBootSortedBiomassMatrix,97.5)'),'-','Color',CIcolor,'LineWidth',2) %upper bound of raw data
    x2 = [BSbins, fliplr(BSbins)];
    inBetween = [log10(prctile(orderedBootSortedBiomassMatrix,2.5)), fliplr(log10(prctile(orderedBootSortedBiomassMatrix,97.5)))];
    fill(x2, inBetween,CIcolor,'EdgeColor','none');
        %plot confidence bounds of Gaussian mixture model
    plot(BSbins,lowerGauss1Pred-LogOffset,'b','LineWidth',2);
    plot(BSbins,upperGauss1Pred-LogOffset,'b','LineWidth',2);
    text(-13,-2.5,{['R^2=' num2str(meanGauss1R2,2) '\pm' num2str(num2str(stdGauss1R2,2))];['AICc=' num2str(mean_relAICcGauss1,2) '\pm' num2str(num2str(std_relAICcGauss1,2))]},'fontsize',16,'Color','b')
    if BA==1
        ylim(log10([0.0002 3000])) %for log body size-log biomass 200 [1E-7 2000]
    elseif BA==2
        ylim([-13 18]) %for log body size-log abundance
    end
    yt=yticks;
    if BA==2
        yyaxis right
        ylim([-13 18]-0.454)
        yticks(yt);
        yticklabels(round(yt-0.454,1)); %normalized biomass
    end
    
    subplot(3,5,(k-1)*5+3)
    hold on
    %ylim(log10([0.0002 3000])) %for log body size-log biomass 200 [1E-7 2000]
    plot(BSbins,log10(prctile(orderedBootSortedBiomassMatrix,2.5)'),'-','Color',CIcolor,'LineWidth',2) %lower bound of raw data
    plot(BSbins,log10(prctile(orderedBootSortedBiomassMatrix,97.5)'),'-','Color',CIcolor,'LineWidth',2) %upper bound of raw data
    x2 = [BSbins, fliplr(BSbins)];
    inBetween = [log10(prctile(orderedBootSortedBiomassMatrix,2.5)), fliplr(log10(prctile(orderedBootSortedBiomassMatrix,97.5)))];
    fill(x2, inBetween,CIcolor,'EdgeColor','none');
        %plot confidence bounds of Gaussian mixture model
    plot(BSbins,lowerGauss2Pred-LogOffset,'b','LineWidth',2);
    plot(BSbins,upperGauss2Pred-LogOffset,'b','LineWidth',2);
    text(-13,-2.5,{['R^2=' num2str(meanGauss2R2,2) '\pm' num2str(num2str(stdGauss2R2,2))];['AICc=' num2str(mean_relAICcGauss2,2) '\pm' num2str(num2str(std_relAICcGauss2,2))]},'fontsize',16,'Color','b')
    if BA==1
        ylim(log10([0.0002 3000])) %for log body size-log biomass 200 [1E-7 2000]
    elseif BA==2
        ylim([-13 18]) %for log body size-log abundance
    end
    yt=yticks;
    if BA==2
        yyaxis right
        ylim([-13 18]-0.454)
        yticks(yt);
        yticklabels(round(yt-0.454,1)); %normalized biomass
    end
    
    
    subplot(3,5,(k-1)*5+4)
    hold on
    %ylim(log10([0.0002 3000])) %for log body size-log biomass 200 [1E-7 2000]
    plot(BSbins,log10(prctile(orderedBootSortedBiomassMatrix,2.5)'),'-','Color',CIcolor,'LineWidth',2) %lower bound of raw data
    plot(BSbins,log10(prctile(orderedBootSortedBiomassMatrix,97.5)'),'-','Color',CIcolor,'LineWidth',2) %upper bound of raw data
    x2 = [BSbins, fliplr(BSbins)];
    inBetween = [log10(prctile(orderedBootSortedBiomassMatrix,2.5)), fliplr(log10(prctile(orderedBootSortedBiomassMatrix,97.5)))];
    fill(x2, inBetween,CIcolor,'EdgeColor','none');
        %plot confidence bounds of Gaussian mixture model
    plot(BSbins,lowerGauss3Pred-LogOffset,'b','LineWidth',2);
    plot(BSbins,upperGauss3Pred-LogOffset,'b','LineWidth',2);
    text(-13,-2.5,{['R^2=' num2str(meanGauss3R2,2) '\pm' num2str(num2str(stdGauss3R2,2))];['AICc=' num2str(mean_relAICcGauss3,2) '\pm' num2str(num2str(std_relAICcGauss3,2))]},'fontsize',16,'Color','b')
    if BA==1
        ylim(log10([0.0002 3000])) %for log body size-log biomass 200 [1E-7 2000]
    elseif BA==2
        ylim([-13 18]) %for log body size-log abundance
    end
    yt=yticks;
    if BA==2
        yyaxis right
        ylim([-13 18]-0.454)
        yticks(yt);
        yticklabels(round(yt-0.454,1)); %normalized biomass
    end
    
    subplot(3,5,(k-1)*5+5)
    hold on
    %ylim(log10([0.0002 3000])) %for log body size-log biomass 200 [1E-7 2000]
    plot(BSbins,log10(prctile(orderedBootSortedBiomassMatrix,2.5)'),'-','Color',CIcolor,'LineWidth',2) %lower bound of raw data
    plot(BSbins,log10(prctile(orderedBootSortedBiomassMatrix,97.5)'),'-','Color',CIcolor,'LineWidth',2) %upper bound of raw data
    x2 = [BSbins, fliplr(BSbins)];
    inBetween = [log10(prctile(orderedBootSortedBiomassMatrix,2.5)), fliplr(log10(prctile(orderedBootSortedBiomassMatrix,97.5)))];
    fill(x2, inBetween,CIcolor,'EdgeColor','none');
        %plot confidence bounds of Gaussian mixture model
    plot(BSbins,lowerGauss4Pred-LogOffset,'b','LineWidth',2);
    plot(BSbins,upperGauss4Pred-LogOffset,'b','LineWidth',2);
    text(-13,-2.5,{['R^2=' num2str(meanGauss4R2,2) '\pm' num2str(num2str(stdGauss4R2,2))];['AICc=' num2str(mean_relAICcGauss4,2) '\pm' num2str(num2str(std_relAICcGauss4,2))]},'fontsize',16,'Color','b')
    if BA==1
        ylim(log10([0.0002 3000])) %for log body size-log biomass 200 [1E-7 2000]
    elseif BA==2
        ylim([-13 18]) %for log body size-log abundance
    end
    yt=yticks;
    if BA==2
        yyaxis right
        ylim([-13 18]-0.454)
        yticks(yt);
        yticklabels(round(yt-0.454,1)); %normalized biomass
    end
end

% modelsSelected
% [bootSelectFreq bootSelectedModel]=max(ModelStats(:,7,:));