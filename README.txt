README for The Sizes of Life
Edward W. Tekwa, Katrina A. Catalano, Anna L. Bazzicalupo, Malin L. Pinsky

March 30, 2022

Matlab files:
bodySizeGEVPlot All 2_0cutoff.matbodySizeGEVPlot All no skeleton 2_0cutoff.matbodySizeGEVPlot All ramet 2_0cutoff.matbodySizeGEVPlot_indivGroup.mbodySizeGEVPowerLawsBootstrapped.mfitGEV.mfminsearchbnd.mGEV3pts.mTekwa size biomass data no skeleton.xlsxTekwa size biomass data original.xlsxTekwa size biomass data ramet.xlsx


1.) To fit distributions to each taxonomic group and plot sum size-biomass spectrum, run “bodySizeGEVPlot_indivGroup.m”. This file can be used to generate Figures 1-3 and S2-3. Uncomment line 1, 2, or 3 to load different .mat files that contain the size biomass data with different assumptions. Adjust within biological group minimum and maximum size truncations (lines 11-12). Change cumulative size-biomass spectrum plotting to exclude different biological groups (line 137).

2.) To run bootstrapped linear and Gaussian mixture regressions on data prepared by step above, run "bodySizeGEVPowerLawsBootstrapped.m". This file can be used to generate Figures 4-5 and S4. Change line 7 to specify size-biomass or size-abundance analyses.

3.) The .xlsx files contains the raw data and calculations for group biomass and size estimates.

4.) Other .m files are functions called by the main codes in 1-2.