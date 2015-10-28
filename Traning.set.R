
## This function allows you to separate out a proportion of individuals from each population to create working and training
  ## datasets. (Training sets brah!)

GenePopData <- read.table("Test_GenePop_Data.txt", header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)

GenePop <- GenePopData 

prop.keep = 0.5

training.bra <- function(GenePop, prop.keep, out.name){
  
## GenePop = the genepop file with all loci 
## subs = the loci names of interest
## keep = logical vector which defines whether you want to remove the loci or keep them. 
  ## the default is to keep them (keep=TRUE) assuming you are removing neutral markers and only keeping the subs ("Outliers")
## sPops is the populations of interest. Note these are specified in the order which they appear in the
  ##    original Genepop file. i.e. first pop = 1 second pop = 2  or the text based origin 
  ##    Examples Numeric: sPop=c(1,3,4,7) text: sPop=c("BMR", "GRR","GHR","TRS")

## Function for inserting rows
  insert.vals <- function(Vec,breaks,newVal){
    break.space <- 1:(length(breaks))
    breaks <- breaks+break.space-1 #To space out the insertion points.
    newvec <- rep(NA,length(Vec)+length(breaks)) #Preallocate memory by creating final dataframe.
    for(i in 1:length(breaks)){newvec[breaks[i]]=newVal} #Insert added rows into new dataframe>
    x <- 1:length(newvec)
    x <- x[-(breaks)] #Finding the rows of the new dataframe that will receive old rows
    for(i in 1:length(Vec)){newvec[x[i]]=Vec[i]} 
    return(newvec)}
 
  #Libraries
  
  #Check to make sure the packages required are there
  packages <- c("dplyr", "tidyr", "stringr") ## which packages do we need?
  if (length(setdiff(packages, rownames(installed.packages()))) > 0) { ### checks that the required packages are among those 
    ## installed by comparing the differences in the length of a vector of items included in list (i.e. is package among all installed)
    install.packages(setdiff(packages, rownames(installed.packages())))  ## if a package is not installed, insall it
  } ### this will only work if someone has the CRAN mirror set (I would assume everyone would?)

  #load each library
    require(dplyr)
    require(tidyr)
    require(stringr)
  
## Stacks version information
    stacks.version <- GenePop[1,] # this could be blank or any other source. First row is ignored by GenePop

#Remove first label of the stacks version
    GenePop <- as.vector(GenePop)
    GenePop <- GenePop[-1,]

#Add an index column to Genepop and format as a dataframe
    GenePop <- data.frame(data=GenePop,ind=1:length(GenePop))
    GenePop$data <- as.character(GenePop$data)

#ID the rows which flag the Populations
    Pops  <-  which(GenePop$data == "Pop" | GenePop$data =="pop" | GenePop$data == "POP")
    npops  <-  1:length(Pops)

## Seperate the data into the column headers and the rest
    ColumnData <- GenePop[1:(Pops[1]-1),"data"]
    snpData <- GenePop[Pops[1]:NROW(GenePop),]

#Get a datafile with just the snp data no pops
    tempPops <- which(snpData$data=="Pop")
    snpData <- snpData[-tempPops,]

#Seperate the snpdata
#First we pull out the population data which follows
#"TEXT ,  "
    temp <- separate(snpData,data,into=c("Pops","snps"),sep=",")
    temp$snps <- substring(temp$snps,3) # delete the extra spaces at the beginning
    temp2 <- as.data.frame(do.call(rbind, strsplit(temp$snps," "))) #split characters by spaces
  
    #Contingency to see if R read in the top line as the "stacks version"
    if (length(temp2)!=length(ColumnData)){colnames(temp2) <- c(stacks.version,ColumnData)}
    if (length(temp2)==length(ColumnData)){colnames(temp2) <- ColumnData}
    if (length(temp2)!=length(ColumnData)){stacks.version="No stacks version specified"}
    
## Get the Alpha names from the 
    NamePops=temp[,1] # Names of each
    NameExtract=str_extract(NamePops, "[A-Z]+" ) # extract the text from the individuals names to denote population

## Now add the population tags using npops (number of populations and Pops for the inter differences)
    tPops <- c(Pops,NROW(GenePop))
      PopIDs <- NULL
          for (i in 2:length(tPops)){
            hold <- tPops[i]-tPops[i-1]-1
            if(i==length(tPops)){hold=hold+1}
            pophold <- rep(npops[i-1],hold)
            PopIDs <- c(PopIDs,pophold)
          }
    
    

 temp2$Pop <- PopIDs;rm(hold,pophold,tPops,PopIDs) ### make a vector to denote the populations
  
    
    PopLengths <- table(temp2$Pop)
    
    ## Need to be able to tell what row each individual is in, and what population it is
    ind.vector = c(1:nrow(temp)) ### make a vector that is the number of individuals
    ind.matrix = data.frame(temp2$Pop, ind.vector) ## add populatuions to that
    
   
    
    ### How many populations are there? - need to get a random proportion from each
        ## also will show what numbers can choose between to draw a random sample
        ## i.e. Population 1 is between 1 and 20, because the first individial of Pop2 is at position 21
      
      ## get t
    make.random.groups = (ind.matrix$ind.vector[!duplicated(ind.matrix$temp2.Pop)])
    make.random.groups <- c(make.random.groups, (max(ind.matrix)+1)) ### add the last individual to the list 
      ## so the numbering for the last Population is correct

    ## Will use Sample - which draws between a min and max
    ## will use the numbers in the vector describing the populations to get the first and last indv. in each pop
    min.max.out <- NULL 
    for(i in 1:(length(make.random.groups)-1)){
          a = make.random.groups[i] ## from the first position in population i
          b = make.random.groups[i+1]-1 ## to the first position in population (i+1) - one pisition
          min.max = c(a, b) ## just mash those together
          min.max.out = rbind(min.max.out, min.max) ## let them marinate
      
    } 
    ## Now - have the first and last indiv. in each pop - so chose a random sample of indv. between them
    
    sample.per.pop <- rep(NA, length(PopLengths))
    
    for(i in 1:length(PopLengths)){
      sample.per.pop[i] = round((prop.keep*PopLengths[i]), digits = 0)
      
      
    }
    
   indiv.out <- NULL
      for(i in 1:(length(sample.per.pop))){
   
          from <- min(min.max.out[i,])
          to <- max(min.max.out[i,])
          indv.keep.hold <- sample(to:from, sample.per.pop[i])
          indiv.out <- c(indiv.out, indv.keep.hold)
          }
      
    training.brah.vector <- indiv.out
    working.brah.vector <- setdiff(ind.vector, training.brah.vector)
    
    working.out <- temp2[working.brah.vector,]
    training.out <- temp2[training.brah.vector,]
    
    drop.last.working <- ncol(working.out)
    drop.last.training <- ncol(training.out)
    
    

#Now recompile the GenePop format
    
    #the number of individuals for all popualtions but the last (Pop tagged to the end)
    PopLengths.working <- table(working.out$Pop)
    PopLengths.training <- table(training.out$Pop)
    
   # working.out <- working.out[,-drop.last.working]
  #  training.out <- training.out[,-drop.last.training]
    
    if(length(table(working.out$Pop))==2){PopPosition.working = PopLengths.working+1}
    
    if(length(table(working.out$Pop))>2){ 
          PopPosition.working <- c(PopLengths.working[1]+1,rep(NA,(length(PopLengths.working)-1)))
          for (i in 2:length(PopLengths.working)){
            PopPosition.working[i] <- PopLengths.working[i]+PopPosition.working[i-1]
          }
    }
    
    
    
    if(length(table(training.out$Pop))==2){PopPosition.training = PopLengths.training+1}
    
    if(length(table(training.out$Pop))>2){ 
          PopPosition.training <- c(PopLengths.training[1]+1,rep(NA,(length(PopLengths.training)-1)))
          for (i in 2:length(PopLengths.training)){
            PopPosition.training[i] <- PopLengths.training[i]+PopPosition.training[i-1]
          }
    }
    

    working.out <- working.out[,-ncol(working.out)]
    training.out <- training.out[,-ncol(training.out)]
    
   
    
    # paste together the Loci as one long integer seperated for each loci by a space
    Loci.working <- do.call(paste,c(working.out[,], sep=" "))
    Loci.training <- do.call(paste,c(training.out[,], sep=" "))
    
    #Grab the Population tags that each invididual had following the format ID_,__
    PopVec.working <- paste(gsub(pattern = " ",replacement = "",working.out$Pop)," ,  ",sep="")
    PopVec.training <- paste(gsub(pattern = " ",replacement = "",training.out$Pop)," ,  ",sep="")
    
    #Paste these to the Loci
    Loci.working <- paste(PopVec.working,Loci.working,sep="")
    Loci.training <- paste(PopVec.training,Loci.training,sep="")
    
    where.pops.training <- data.frame(PopLengths.training)[,2]
    where.pops.working <- data.frame(PopLengths.working)[,2]
      
    #Insert the value of "Pop" which partitions the data among populations
    Loci.training <- insert.vals(Vec=Loci.training, breaks = PopPosition.training, newVal = "Pop")
    Loci.working <- insert.vals(Vec=Loci.working,breaks=PopLengths.working,newVal="Pop")
    
    #Add the first "Pop" label
    Loci.training <- c("Pop",Loci.training) 
    Loci.working <- c("Pop", Loci.working)
    ## Add the column labels and the stacks version
    
    ## Add the column labels and the stacks version
    
    #out.name <- deparse(substitute(GenePop))
    out.name.working <- paste(out.name, "working.txt", sep =".")
    out.name.training <- paste(out.name, "training.txt", sep =".")
    
    Output.working <- c(stacks.version,names(working.out)[-length(names(working.out))],Loci.working)
    Output.training <- c(stacks.version,names(training.out)[-length(names(training.out))],Loci.training)
    
    # Save the file - working
    write.table(Output.working,out.name.working,col.names=FALSE,row.names=FALSE)
    write.table(Output.training,out.name.training,col.names=FALSE,row.names=FALSE)
    
    
    
    
} #End function

training.bra(GenePop = GenePopData, prop.keep = 0.5, out.name = "Test")
