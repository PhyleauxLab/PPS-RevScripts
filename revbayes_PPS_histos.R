## used to generate histograms of the pps test statistics

path = "FULL_PATH_HERE/combined_results/"
file.names <- dir(path, pattern ="*.csv")


stats = c("Lower.1.tailed",	"Upper.1.tailed",	"Two.tailed",	"Midpoint",	"Effect.Size")

x <- 1

for(i in 1:length(file.names)){

pp_data <- read.csv(paste(path, file.names[i], sep=""), header = TRUE)
statName <- strsplit(file.names[i], "_", fixed=TRUE)[[1]][-1]
statName <- strsplit(statName, "[.]")[[1]][1]
pdf(paste("plots/",statName, "_plots.pdf", sep=""))
hist(sapply(pp_data["Lower.1.tailed"],as.numeric), main=paste(toupper(statName), "Lower 1 tailed"), xlab="P-Values", col="grey")
hist(sapply(pp_data["Upper.1.tailed"],as.numeric), main=paste(toupper(statName), "Upper 1 tailed"), xlab="P-Values", col="grey")
hist(sapply(pp_data["Two.tailed"],as.numeric), main=paste(toupper(statName), "Two-tailed"), xlab="P-Values", col="grey")
hist(sapply(pp_data["Midpoint"],as.numeric), main=paste(toupper(statName), "Midpoint"), xlab="P-Values", col="grey")
hist(sapply(pp_data["Effect.Size"],as.numeric), main=paste(toupper(statName), "Effect Size"), xlab="Effect Size", col="grey")
dev.off()
} 