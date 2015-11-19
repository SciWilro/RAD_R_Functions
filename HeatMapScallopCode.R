## Scallop analysis ###

# Libraries -----------------
library(ggplot2)
library(dplyr)
library(boa)
library(scales)
library(VennDiagram)
library(RCurl)
#library(gstudio)

# Functions -----------------------

#function to grab R code of Github using raw urls
SourceGitFunc <- function(url)
{
  
  require(RCurl)
  script <- getURL(url, ssl.verifypeer = FALSE)
  eval(parse(text = script),envir=.GlobalEnv)
}

#Source the function which will create the plot
SourceGitFunc("https://raw.githubusercontent.com/rystanley/RAD_R_Functions/master/AlleleFreq.R")

# Load data ------------------

GenePopData <- read.table("GenePopFixed.txt",
                          header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)


Outliers <- read.table("AllPopsOutliers.txt",
                       header = TRUE, sep = "\t")

Outliers <- as.character(Outliers$Loci)

Neutral <- read.table("AllPopsNeutral.txt",
                      header=TRUE, sep="\t")

Neutral  <- as.character(Neutral$Loci)

# Get the order of populations --------------
PopOrder <- rev(c("bon","ltb","mgd","nts","psb","bof","ssm","gmi","ssb","gmo","geo","mda")   )


### Create the allele frequency plot --------------
Scallop_Heat_Outlier=AlleleFreqHeatMap(GenePopData,subs = Outliers,keep = TRUE,POP="CHAR",refPop="bon",
                               OrderPops=PopOrder,optimizer = TRUE,standardize = TRUE,plot=TRUE)

Scallop_Heat_Outlier # show the plot output

ggsave("Scallop_Outlier_AlleleFreqHeatMap.png",Scallop_Heat_Outlier,
       width=12,height=8)

## Create the allele frequency plot
Scallop_Heat_Neutral=AlleleFreqHeatMap(GenePopData,subs = Neutral,keep = TRUE,POP="CHAR",refPop="bon",
                                       OrderPops=PopOrder,optimizer = TRUE,standardize = TRUE,plot=TRUE)

Scallop_Heat_Neutral # show the plot output

ggsave("Scallop_Neutral_AlleleFreqHeatMap.png",Scallop_Heat_Neutral,
       width=12,height=8)



### Use the clinal Neutral markers -------------
Cline_Neutral <- read.csv("Cline_Neutral_Loci.csv")
Cline_Neutral <- as.character(Cline_Neutral$Neutral)
Cline_Neutral <- gsub("[a-zA-Z]+","",Cline_Neutral) # remove the allele letters

Scallop_Heat_NeutralCline=AlleleFreqHeatMap(GenePopData,subs = Cline_Neutral,keep = TRUE,POP="CHAR",refPop="bon",
                                            OrderPops=PopOrder,optimizer = TRUE,standardize = TRUE,plot=TRUE)

Scallop_Heat_NeutralCline # show the plot output

ggsave("Scallop_Neutral_Cline_AlleleFreqHeatMap.png",
       Scallop_Heat_NeutralCline,width=12,height=8)

### Use the clinal outlier markers -------------
Cline_Outlier <- read.csv("Cline_Outlier_Loci.csv")
Cline_Outlier <- as.character(Cline_Outlier$Outlier)
Cline_Outlier <- gsub("[a-zA-Z]+","",Cline_Outlier) # remove the allele letters

Scallop_Heat_OutlierCline=AlleleFreqHeatMap(GenePopData,subs = Cline_Outlier,keep = TRUE,POP="CHAR",refPop="bon",
                                            OrderPops=PopOrder,optimizer = TRUE,standardize = TRUE,plot=TRUE)

Scallop_Heat_OutlierCline # show the plot output

ggsave("Scallop_Outlier_Cline_AlleleFreqHeatMap.png",
       Scallop_Heat_OutlierCline,width=12,height=8)

