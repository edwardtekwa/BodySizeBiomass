load('bodySizeGEVPlot All 2_0cutoff.mat'); %default data
%load('bodySizeGEVPlot All ramet 2_0cutoff.mat'); %use ramet instead of genet for body size
%load('bodySizeGEVPlot All no skeleton 2_0cutoff.mat'); %exclude skeleton mass and subterranean microbes

scrsz = get(0,'ScreenSize');
set(0, 'DefaultAxesFontSize', 18)

CIcolor=[140 140 140]./255;


minSizeOffset=2; %log10 offset to reported minimum size for biomass distribution truncation
maxSizeOffset=0; %log10 offset to reported maximum size for biomass distribution truncation

BSlims=[-18 11]; %in log10(gC) scale
res=40; %40 or 1 points per log size bin
numBoots=100; %number of biomass bootstraps (1000)
numBins=(BSlims(2)-BSlims(1))*res+1; %log10(BS)=[-14 8] %44
BSbins=[BSlims(1):(BSlims(2)-BSlims(1)-1)/((BSlims(2)-BSlims(1)-1)*res):BSlims(2)];
BiomassMatrix=zeros(size(BSAllGroups,1),numBins); %col: body size bin, row: organismal groups
sortedBiomassMatrix=zeros(size(BSAllGroups,1),numBins);
BiomassPerSMatrix=zeros(size(BSAllGroups,1),numBins); %col: body size bin, row: organismal groups
AbundanceMatrix=zeros(size(BSAllGroups,1),numBins); %col: body size bin, row: organismal groups
Taxa=zeros(1,numBins); %number of taxa(groups) in each bin
% Xs=zeros(size(BSGroups,1),3);
% fvals=zeros(size(BSGroups,1),1);
% CDFs=zeros(size(BSGroups,1),3);

%sort groups(taxa) from least to most biomass
[~,idx] = sort(BSAllGroups(:,4)); %omit if using original biomass ranked order
sortedBSGroups = BSAllGroups(idx,:);
sortedGroups = AllGroups(idx);
sortedGroupFoldUncert = BSAllGroups(idx,6);
sortedGroupBiomass = BSAllGroups(idx,4);

bootGroupBiomass=[]; %generate bootstrap group biomass samples
for boot=1:numBoots
    bootGroupBiomass=[bootGroupBiomass sortedGroupBiomass.*2.^(normrnd(0,log2(sortedGroupFoldUncert)/1.96))];
end

for s=1:size(sortedBSGroups,1) %fit skew normal distributions to each taxon, record distribution parameters in Xs
    sortedGroups{s}
    [X,fval,CDFpts]=fitGEV([sortedBSGroups(s,1) sortedBSGroups(s,3) sortedBSGroups(s,2)],minSizeOffset,maxSizeOffset); %fit truncated GEV biomass distribution
    Xs(s,:)=X;
    fvals(s)=fval;
    CDFs(s,:)=CDFpts;
    pd=makedist('gev','k',X(1),'sigma',X(2),'mu',X(3));
    td=truncate(pd,log10(sortedBSGroups(s,1))-minSizeOffset,log10(sortedBSGroups(s,2))+maxSizeOffset); %truncate GEV distribution
    for bi=1:numBins %record biomass density of species s at size bi
        LogSizeCenter=BSbins(bi); %center of log size bin
        sortedBiomassMatrix(s,bi)=pdf(td,LogSizeCenter)*sortedBSGroups(s,4);
    end
end

All_orderedBootSortedBiomassMatrix=zeros(numBoots,numBins,36);
for s=1:36
    BootSortedBiomassMatrix=zeros(numBoots,numBins); %store boostrapped biomass by size bin
    Taxa=zeros(numBoots,numBins);
    pd=makedist('gev','k',Xs(s,1),'sigma',Xs(s,2),'mu',Xs(s,3));
    td0=truncate(pd,log10(sortedBSGroups(s,1))-minSizeOffset,log10(sortedBSGroups(s,2))+maxSizeOffset); %truncate GEV distribution
    for boot=1:numBoots %option 3: using skew normal or normal with random size mode
        randMeanSize=10^random(td0); %size output in arithmetic scale
        [Xboot,fvalboot,CDFptsboot]=fitGEV([sortedBSGroups(s,1) randMeanSize sortedBSGroups(s,2)],minSizeOffset,maxSizeOffset);
        pd=makedist('gev','k',Xboot(1),'sigma',Xboot(2),'mu',Xboot(3));
        td=truncate(pd,log10(sortedBSGroups(s,1))-minSizeOffset,log10(sortedBSGroups(s,2))+maxSizeOffset); %truncate GEV distribution
        for bi=1:numBins %construct bootstrapped cumulative biomass density function one bin at a time
            LogSizeCenter=BSbins(bi); %center of log size bin
            biomassProbDistr=pdf(td,LogSizeCenter);
            BootSortedBiomassMatrix(boot,bi)=BootSortedBiomassMatrix(boot,bi)+biomassProbDistr*bootGroupBiomass(s,boot);
            if biomassProbDistr*bootGroupBiomass(s,boot)>eps %threshold in GtC for group to count in bin (eps or 0.0002)
                Taxa(boot,bi)=Taxa(boot,bi)+1;
            end
            
        end
    end
    orderedBootSortedBiomassMatrix=sort(BootSortedBiomassMatrix);
    All_orderedBootSortedBiomassMatrix(:,:,s)=orderedBootSortedBiomassMatrix;
end


% %reorder plotBiomassMatrix to match original colours if incorrect
% [~,idxOriginal] = sort(BSAllGroupsOriginal(:,4));
% for s=1:36
%    newPos(s)=find(idx==idxOriginal(s)); 
% end
% newsortedBSGroups=sortedBSGroups(newPos,:);
% plotBiomassMatrix = plotBiomassMatrix(newPos,:);
% 
% %reorder All_orderedBootSortedBiomassMatrix to match current group biomass
% %rank if incorrect
% All_orderedBootSortedBiomassMatrix=All_orderedBootSortedBiomassMatrix(:,:,newPos);
% sortedGroups=sortedGroups(newPos);
% sortedBSGroups=sortedBSGroups(newPos,:);


%Figure 2: size biomass spectra per biological group
figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/1.6 scrsz(3)/2]);
tickLocations=[1:(numBins-1)/((BSlims(2)-BSlims(1))/2):numBins];
for s=1:36
    subplot(6,6,s)
    orderedBootSortedBiomassMatrix=All_orderedBootSortedBiomassMatrix(:,:,s);
    set(gca, 'YScale', 'log')
    hold on
    plot(prctile(sum(All_orderedBootSortedBiomassMatrix,3),50)','LineWidth',1,'Color',CIcolor)
    plot(prctile(orderedBootSortedBiomassMatrix,2.5)','k--','LineWidth',1)
    plot(prctile(orderedBootSortedBiomassMatrix,50)','k','LineWidth',2)
    plot(prctile(orderedBootSortedBiomassMatrix,97.5)','k--','LineWidth',1)
    xlim ([1 numBins])
    xticks(tickLocations(1:2:end)) %40
    xticklabels(BSbins(tickLocations(1:2:end)))
    ylim ([0.000001 500])
    yticks([1E-6 1E-4 0.01 1 100])
    title(sortedGroups(s))
end

%extra figure: plotting mean versus median biomass
figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/4 scrsz(3)/5]);
set(gca, 'YScale', 'log')
hold on
plot(mean(sum(All_orderedBootSortedBiomassMatrix,3)),'k','LineWidth',2)
plot(prctile(sum(All_orderedBootSortedBiomassMatrix,3),2.5)','k--','LineWidth',1)
plot(prctile(sum(All_orderedBootSortedBiomassMatrix,3),50)','k','LineWidth',4)
plot(prctile(sum(All_orderedBootSortedBiomassMatrix,3),97.5)','k--','LineWidth',1)
xlim ([1 numBins])
xticks(tickLocations(1:2:end)) %40
xticklabels(BSbins(tickLocations(1:2:end)))
ylim ([0.001 1000])
yticks([0.01 1 100])
xlabel 'body size [log_{10} g]'
ylabel 'biomass [Gt]'

set(0,'defaultAxesColorOrder',[[0 0 0]; [77 190 238]./255]);
ColorMapChoice='colorcube'; %jet (overall) or parula (within habitat/trophic level)


%Figure 1 or 3
%select data subset to exclude from plot:
NotTypes=[]; %1. terrestrial producer, 2. terrestrial consumer, 3. marine producer, 4. marine consumer, 5. subterranean
typeidx=[];
plotBiomassMatrix=reshape(median(All_orderedBootSortedBiomassMatrix),[],36)';
plot_orderedBootSortedBiomassMatrix=All_orderedBootSortedBiomassMatrix;
if ~isempty(NotTypes)
    for i=1:length(NotTypes)
        typeidx = [typeidx; find(sortedBSGroups(:,5)==NotTypes(i))];
    end
    plotBiomassMatrix(typeidx,:)=0;
    plot_orderedBootSortedBiomassMatrix(:,:,typeidx)=0;
end

figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2 scrsz(4)/3.5]);
CIcolor=[140 140 140]./255;

subplot(1,2,1)
colormap(colorcube)
barObj=bar(plotBiomassMatrix',1,'stacked','EdgeColor','none','FaceColor','flat');
for k = 1:size(plotBiomassMatrix',2)
    barObj(k).CData = 37-k;
end
set(gca, 'YScale', 'log')
hold on
plot(prctile(sum(All_orderedBootSortedBiomassMatrix,3),50)','LineWidth',1,'Color',CIcolor)
plot(prctile(sum(plot_orderedBootSortedBiomassMatrix,3),2.5)','k--','LineWidth',2)
plot(prctile(sum(plot_orderedBootSortedBiomassMatrix,3),97.5)','k--','LineWidth',2)
logBiomass=log10(sum(plotBiomassMatrix));
logBiomass(logBiomass==-Inf)=NaN;
RealBiomassPowerLaw=fitlm(BSbins,logBiomass) %with x-value being log body size !!should limit to observed min and max body sizes
biomassPowerLaw=fitlm(1:numBins,logBiomass); %with x-value being bin number
xticks(tickLocations) %40
xticklabels(BSbins(tickLocations))
xlabel 'body size [log_{10}g]'
ylabel 'biomass [Gt]'
ylim ([0.0001 10000])

subplot(1,2,2)
colormap(colorcube)
hold on
barObj=bar(plotBiomassMatrix',1,'stacked','EdgeColor','none','FaceColor','flat');
for k = 1:size(plotBiomassMatrix',2)
    barObj(k).CData = 37-k;
end
xticks(tickLocations) %40
xticklabels(BSbins(tickLocations))
xlabel 'body size [log_{10}g]'
ylabel 'biomass [Gt]'
box on

%Figure S2
figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/5 scrsz(4)/3.5]);
sizeRangePowerLaw=fitlm(log10(BSAllGroups(:,3)),log10(BSAllGroups(:,2))-log10(BSAllGroups(:,1)))
hold on
plot(sizeRangePowerLaw)
scatter(log10(BSAllGroups(:,3)),log10(BSAllGroups(:,2))-log10(BSAllGroups(:,1)),60,'k','filled')
xlabel('mean body size [log_{10}g]','Interpreter','tex')
ylabel('body size range [log_{10}g]','Interpreter','tex')
title('')
legend off
meanSizeRange=mean(log10(BSAllGroups(:,2))-log10(BSAllGroups(:,1)))
sdSizeRange=std(log10(BSAllGroups(:,2))-log10(BSAllGroups(:,1)))
