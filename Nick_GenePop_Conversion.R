#Nick Example

source("GenoPopConvert.R")

##load in the genepop file with all of the snp data -----------
GenePopData <- read.table("E:/Nick Examples/crab_genepop.txt", header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)

#Grab the outlier list. Note the function requires a vector of characters. 
Outlist=as.vector(read.table("E:/Nick Examples/CrabFinal_outliers.txt",header=T))
Outlist=as.character(Outlist[,1]) #string as a vector


## Create Genepop file with the outliers ------
subset.GenePop(GenePop=GenePopData,
               subs=Outlist,
               dir="Nick_Outliers.txt")

## Create Genepop file with the neutral markers ---------
subset.GenePop(GenePop=GenePopData,
               subs=Outlist,
               dir="Nick_Neutral.txt",keep=FALSE)