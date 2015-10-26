#This function will subset a regular formatted GenePop file. You specify the outliers you would like to 
#subet from the GenePop formated txt file and this function will re-format the text file to  only include the loci 
#required. 

source("GenoPopConvert.R") # source the function (*Note this is relative to the directory)

#GenePopData <- read.table("Test_GenePop_Data.txt", header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
GenePopData <- read.table("Sal.txt", header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)


#Subset out just the outliers
subset.GenePop(GenePop=GenePopData,
               subs=c("190_12","783_49","3075_62","5322_19","7158_20","10139_54","11492_24","13344_39","14239_50"),
               dir="Example_Parsed_outliers_GenePop.txt")

#Subset out just the neutral markers
subset.GenePop(GenePop=GenePopData,
               subs=c("190_12","783_49","3075_62","5322_19","7158_20","10139_54","11492_24","13344_39","14239_50"),
               dir="Example_Parsed_neutral_GenePop.txt",keep=FALSE)

#Subset out populations based on numeric positions
subset.GenePop(GenePop=GenePopData,
               subs=c("190_12","783_49","3075_62","5322_19","7158_20","10139_54","11492_24","13344_39","14239_50"),
               dir="Example_Parsed_outliers_pops_GenePop_Numeric.txt",sPop=c(1,5,6,8))

#Subset out populations based on text based population names
subset.GenePop(GenePop=GenePopData,
               subs=c("190_12","783_49","3075_62","5322_19","7158_20","10139_54","11492_24","13344_39","14239_50"),
               dir="Example_Parsed_outliers_pops_GenePop_Alpha.txt",sPop=c("BDN","DHB","GNR","LHR"))

#Subset out populations based on text based population names but keeping all loci
subset.GenePop(GenePop=GenePopData,
               dir="Example_Parsed_pops_GenePop_Alpha.txt",sPop=c("BDN","DHB","GNR","LHR"))