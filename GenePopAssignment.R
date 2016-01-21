TrainingGenePop <- function(GenePopData,perc=0.5,Dir,ranSeed=NULL){

  #This function will split a GenePop file into two randomly divided parts based on the percentage.
  #This percentage defines the stratified (population) random sub-sampling to be performed. 
  
  #GenePopData = the genepop file with all loci
  #Read in genepop data
        #EXAMPLE: GenePopData <- read.table("microsat_base",header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
                            
  #perc = the percentage 0-1 (DEFAULT 0.5) of individuals to be included in the 'traning' dataset note that 1-percent will be assigned to the 'assignment' dataset
  #dir = is the directory where the files are to be stored. This must be specified
  #ranSeed = is the randomization seed to be used, this is important if you want to repeat an analysis. Default this is 23. If not specified
  
  # Is the directory specified? If not default to current directory
      if(!exists("Dir",mode = "function")){Dir <- paste0(getwd(),"/")}
  
  #Check to make sure the packages required are there
      packages <- c("dplyr", "tidyr", "stringr") ## which packages do we need?
      if (length(setdiff(packages, rownames(installed.packages()))) > 0) { ### checks that the required packages are among those 
         install.packages(setdiff(packages, rownames(installed.packages())))  ## if a package is not installed, insall it
      } 
  
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
  
    #libraries
      require(dplyr)
      require(tidyr)
      require(stringr)

    #Add an index column to Genepop and format as a dataframe
      GenePop <- as.data.frame(GenePopData)
      GenePop$ind <- rownames(GenePop) #index variable
      colnames(GenePop)[1]="data" # rename column

    #ID the rows which flag the Populations
      Pops  <-  which(GenePop$data == "Pop" | GenePop$data =="pop" | GenePop$data == "POP")
      npops  <-  1:length(Pops)

    ## Seperate the data into the column headers and the rest
      ColumnData <- GenePop[1:(Pops[1]-1),"data"]
      stacks.version <- ColumnData[1]
      ColumnData <- ColumnData[-1] # remove the first stacks.version line if it is there. 

    #Genetic data which is not the column names (note this also includes the population labels)
      snpData <- GenePop[Pops[1]:NROW(GenePop),]
      row.names(snpData) <- NULL # reset the row names

    #Row numbers with the Pop flag
      tempPops <- which(snpData$data=="Pop"| snpData$data =="pop" | snpData$data == "POP") 
      
    #remove the pop labels
      temp <- snpData[-tempPops,] 

    #seperate the data so each population can be teased apart ---------------
      PopData <- separate(temp,data,into=c("Pops","snps"),sep=",") 
      PopData$snps <- substring(PopData$snps,2) # delete the extra spaces at the beginning
      temp2 <- as.data.frame(do.call(rbind, strsplit(PopData$snps," "))) #split characters by spaces

    #Contingency to see if R read in the top line as the "stacks version"
      if (length(temp2)!=length(ColumnData)){ColumnData <- c(stacks.version,ColumnData)} #there was no stacks version so we add the column name back

    #find the number of samples from each strata (population) and then calcualte the number of individuals to be included for 
    #the "training" and "assignment" datasets
      PopNumber <- as.data.frame(table(PopData$Pops))
      PopNumber$training <- floor(PopNumber$Freq*perc)
      PopNumber$assignment <- PopNumber$Freq-PopNumber$training

    ## stratified random sampling by population ---------------
          assignrows <- NULL;trainrows <- NULL
          for(i in unique(PopData$Pops)){
            
            poprows <- which(PopData$Pops==i) # rows corresponding to that population
            if(length(ranSeed)>0){set.seed(ranSeed)} # set to specified seed
            trows <- base::sample(poprows,PopNumber[which(PopNumber$Var1==i),"training"])
            
            if(length(ranSeed)>0){set.seed(ranSeed)} # set to specified seed
            arows <- base::setdiff(poprows,trows)
            assignrows <- c(assignrows,arows)
            trainrows <- c(trainrows,trows)
            rm(poprows)
          }

    #order the row indexes which will populate the new dataframe ---------------
      assignrows <- assignrows[order(assignrows)]
      trainrows <- trainrows[order(trainrows)]


    #Subset the data for each dataset ---------------
      AssignmentData <- temp[assignrows,]
      TrainingData <- temp[trainrows,]


#find out the row positions where "Pop" needs to be added ---------------
    AssignmentLengths <- PopNumber[-nrow(PopNumber),"assignment"]
    AssignmentPosition <- c(AssignmentLengths[1]+1,rep(NA,(length(AssignmentLengths)-1)))
    for (i in 2:length(AssignmentLengths)){
      AssignmentPosition[i] <- AssignmentLengths[i]+AssignmentPosition[i-1]
    }
    
    TrainingLengths <- PopNumber[-nrow(PopNumber),"training"]
    TrainingPosition <- c(TrainingLengths[1]+1,rep(NA,(length(TrainingLengths)-1)))
    for (i in 2:length(TrainingLengths)){
      TrainingPosition[i] <- TrainingLengths[i]+TrainingPosition[i-1]
    }
    
    #Add "Pop" labels
    AData <- insert.vals(Vec = AssignmentData[,1],breaks = AssignmentPosition,newVal = "Pop")
    TData <- insert.vals(Vec = TrainingData[,1],breaks = TrainingPosition,newVal = "Pop")


    ## Save output -----------------
      AssignmentOutput <- c(paste0("Assignment Genepop file ~",(1-perc)*100,"% of data"),ColumnData,"Pop",AData)
      TrainingOutput <- c(paste0("Training Genepop file ~",perc*100,"% of data"),ColumnData,"Pop",TData)

      write.table(AssignmentOutput,paste0(Dir,"AssignmentGenePop_",(1-perc)*100,".txt"),col.names=FALSE,row.names=FALSE,quote=FALSE)
      write.table(TrainingOutput,paste0(Dir,"TrainingGenePop_",perc*100,".txt"),col.names=FALSE,row.names=FALSE,quote=FALSE)
}