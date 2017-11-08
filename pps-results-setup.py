## just used to copy values from pps output into single summary files for easier processing

import glob



simulatedFiles = glob.glob("results_GTR/inference_pvalues_effectsizes*.csv")

print(simulatedFiles)

stats = ["mean_rf",
"quantile25",
"quantile50",
"quantile75",
"quantile99",
"quantile999",
"mean_tl",
"var_tl",
"entropy"]

datasetName = "gtrnull"

for stat in stats :


	fileName = "combined_results/" + datasetName + "_" + stat + ".csv"
	outFile = open(fileName, "a")
	outFile.write("Lower 1-tailed,Upper 1-tailed,Two-tailed,Midpoint,Effect Size\n")

	for f in simulatedFiles :

		lines = open(f, 'r').readlines()
		for line in lines :
			rec = line.strip().split(",")
			if rec[0] == stat :
				outFile.write(line.strip().replace(stat+",", "").replace(",", ",")+"\n")
	outFile.close()