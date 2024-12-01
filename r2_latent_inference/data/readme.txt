gh146flyaverage.txt & orcoflyaverage.txt are 65-dimensional odor x glomerulus average responses (over lobes & trials) for each fly, obtained by running ORNvsPN_analysis_ALLDATA/plotRandomSubsetPNandORN.m and saving to file

The other files are CSVs filled in by running relevant analysis scripts

OCT-MCH_behavior_persistence.csv:
	- columns obtained by running fly_behavior/PAPERFIGURE_oct_vs_mch_persistence/oct_mch_persistence.m and taking x / y values from the 3 hour behavior persistence figure

ORN_Brp_scores.csv:
	- columns obtained by running IHC/analyze_IHC_brpshort.m and taking x / y values from article Figure 3F, figure #8 in MATLAB script, the measured preference vs. predicted preference using the trained ORN Brp-Short model evaluated on train + test flies

PN_PC2_interp_scores.csv:
	- PN PC2 interpreted (Fig 2C) column: filled in by running PN_analysis/evaluate_allData.m and taking the x values from the measured vs. predicted preference using the interpreted PN PC2 scatterplot (article Figure 2C, figure #18 in MATLAB script)
	- OCT-MCH behavior scores column: y values from above figure
	- PN PC2 scores (Fig 1M) column: filled in by running PN_analysis/evaluateTrainedModelOnTestData.m and taking the x values from the measured vs. predicted preference (article Figure 1M, figure #7 in MATLAB script)

