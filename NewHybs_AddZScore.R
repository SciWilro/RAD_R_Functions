#Add a column to a NH input file to include the known genotype classes (e.g. z1, etc)
setwd("C:/Users/Nick/Desktop")
#Read in the file

NHdata <- read.table("Corr_West.top96_wild.unlinked.S1_R2NH.txt", header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
head(NHdata[1:5,])
#Save the first 5 rows to be added back in at the end.
addbackin<-NHdata[1:5,]
library(tidyr)
#Now read it back in but without the first 5 rows.
NHdata2 <- read.table("Corr_West.top96_wild.unlinked.S1_R2NH.txt", header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE,skip=5)
head(NHdata2)
#Give the column a name for the separate function to work
names(NHdata2)<-"Data2"
Nhdata3<-separate(data = NHdata2,col=Data2,into = c("Individual",rep("SNP",96)),sep = " ")

write.csv(x = Nhdata3,file = "NHdatasplit.csv",row.names = FALSE)

###Next, want to make a csv of your z codes per individual that is separate from the file with all the locus names
#Then merge these together
Nhdata4<-read.csv("NHdatasplit.csv")
Zscorevector<-read.csv("ZvectorforNH.csv")
#Merge
NHfinal<- merge(y=Nhdata4,x=Zscorevector,by="Individual",all=TRUE)
head(NHfinal)

#Now make this an actual file readable by NewHybrids
Loci <- do.call(paste,c(NHfinal[,], sep=" "))
Loci<-gsub("NA"," ",x = Loci)
Loci2<-c(addbackin,Loci)

write.table(Loci2,file = "NHFinalInput.txt",col.names=FALSE,row.names=FALSE,quote=FALSE)
