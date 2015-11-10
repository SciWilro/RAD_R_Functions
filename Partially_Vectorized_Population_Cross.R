#GenePopData <- read.table("Salmon_South_Coast.txt", header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
GenePopData <- read.table("Example_Parsed_outliers_pops_GenePop_Alpha.txt", header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
GenePop <- GenePopData
ColumnData.Dup = NULL
temp = NULL
temp2 = NULL
a = NULL
b = NULL
Col.name.orig = NULL


## Testing values
#GenePopData <- read.table("Test_GenePop_Data.txt", header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)

GenePopData <- read.table("SC_Pops_East_West.txt", header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)


rm(list=ls())

#source("GenoPopConvert.R")
#GenePopData <- read.table("SC_Pops_WildFarmed_Genepop2.txt", header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)

 #loc <- c("MHC_IA_5500-5800_2", "MHC_IA_5500-5800_3", "MHC_IA_5500-5800_3", "PACAP2GA556_V2", "MHC_1B_SNP8_89170-89270")
# subset.GenePop(GenePop=GenePopData, dir="meta.txt", subs=loc)
 

rm(list=ls())
rm(list=ls())
GenePopData <- read.table("t.sub2.txt", header = FALSE, sep = "\t", qu ote = "", stringsAsFactors = FALSE)



 GenePopData <- read.table("SC_Pops_East_West.txt", header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
 
GenePopData <- read.table("SC_Pops_Wild_Farmed.txt", header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
GenePop <- GenePopData
Nmothers = 20
Nfathers = 20
out.name <- "Farm_Wild"
pop.groups <- c("FRM", "WLD")
family.size <- 50

#do.sex <- function(genotype, Nmothers, Nfathers, Noff, pop, F1, F2 = No, BC = Yes){


 # do.sex <- function(genotype, Nmothers, Nfathers, out.name, family.size, pop.groups=NULL){ 

    #Check to make sure the packages required are there
 # packages <- c("dplyr", "tidyr", "stringr", "plyr") ## which packages do we need?
#  if (length(setdiff(packages, rownames(installed.packages()))) > 0) { ### checks that the required packages are among those 
    ## installed by comparing the differences in the length of a vector of items included in list (i.e. is package among all installed)
 #   install.packages(setdiff(packages, rownames(installed.packages())))  ## if a package is not installed, insall it
#  } ### this will only work if someone has the CRAN mirror set (I would assume everyone would?)

  #load each library
    require(plyr)
    require(dplyr)
    require(tidyr)
    require(stringr)
    
  
 # GenePop <- genotype
  
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
    ColumnData <- GenePop[1:(Pops[1]-1),"data"]  ### SNP Names
    snpData <- GenePop[Pops[1]:NROW(GenePop),]  ### Genotypes

#Get a datafile with just the snp data no pops
    tempPops <- which(snpData$data=="Pop"| snpData$data =="pop" | snpData$data == "POP")
    snpData <- snpData[-tempPops,]

#Seperate the snpdata
#First we pull out the population data which follows
#"TEXT ,  "
    temp <- separate(snpData,data,into=c("Pops","snps"),sep=",")
    temp$snps <- substring(temp$snps,3) # delete the extra spaces at the beginning
    temp2 <- data.frame(do.call(rbind, str_extract_all(temp$snps, "[0-9]{3}"))) 
   
## Going to have to break the two alleles of the SNPS apart - this will thus double the number of columns
    ## SO <- will want to have SNP_A and SNP_A2 
ColumnData2 <- ColumnData ## Duplicatet the SNP names
ColumnData2 <- paste(ColumnData2, "2", sep = ".")    ## add .2 to each duplicated name
     
## can't just append the duplicated names to the end of the original names - have to intersperse them
  ## need a way to do so 
places = rep(1:length(ColumnData)*2) ### creates a list of even numbers 2X as long as the number of columns
  ##i.e. the lenght of the original plus the duplicated names 
  ## - this will also mark the position to insert the duplicates
ColumnData.Dup = rep(NA, times = length(ColumnData)*2) ### make an object to feed names into
for(i in 1:length(ColumnData)){ ### for loop to add original, then duplicated name
  a = places[i]-1 ## Original names go first, and are in the odd positions
  b = places[i] ### Duplicated names go second and are in the even posoitons
  Col.name.orig = ColumnData[i] ## Get the name in the ith position
  Col.name.plus2 = ColumnData2[i] ## get the name in the ith positon
  ColumnData.Dup[a] = Col.name.orig ## add the original name to the new vector
  ColumnData.Dup[b] = Col.name.plus2 ## add the duplicate name to the new vector
}

# temp3 <- data.frame(do.call(cbind, strsplit(temp2," "))) #split characters by spaces
  #View(temp2)
  #View(temp3)
    #Contingency to see if R read in the top line as the "stacks version" -- modified to deal with the duplicated SNP names
    if (length(temp2)!=length(ColumnData.Dup)){colnames(temp2) <- c(stacks.version, paste(stacks.version, "2", sep = "."),ColumnData.Dup)}
    if (length(temp2)==length(ColumnData.Dup)){colnames(temp2) <- ColumnData.Dup}
    #if (length(temp2)/2!=length(ColumnData)){stacks.version="No stacks version specified"}
    



## Get the Alpha names from the 
    NamePops=temp[,1] # Names of each
    #NameExtract=str_extract(NamePops, "[A-Z]+" ) # extract the text from the individuals names to denote population
  
   if(length(pop.groups) == 0){ ### If unique grouping IDs ≠ number of "Pop" user must give vector of groupings 
                                ### equal to number of "Pop" or else the function will fail
    NameExtract=str_extract(NamePops, "[A-z]{3}" ) ### if looking at higher order grouping (i.e. pops in 
        # regions) can have more unique coding than "Pop" - will want to remove original names so can
        ## keep track of which unique groupings cross. i.e. Cross by "Pop", but remember ID of parents
   }
    
   
  # extract the text from the individuals names to denote population
  ## Now add the population tags using npops (number of populations and Pops for the inter differences)
    tPops <- c(Pops,NROW(GenePop))
      PopIDs <- NULL
          for (i in 2:length(tPops)){
            hold <- tPops[i]-tPops[i-1]-1
            if(i==length(tPops)){hold=hold+1}
            pophold <- rep(npops[i-1],hold)
            PopIDs <- c(PopIDs,pophold)
          }
    
    temp2$Pop <- PopIDs;
    #rm(hold,pophold,tPops,PopIDs) ### make a vector to denote the populations

    if(length(pop.groups)!=0){
     hold.names=str_extract(NamePops, "[A-z]{3}" ) ## This may need to be improved in published version
        for(i in 1:length(unique(PopIDs))){
          u.ID.no <- unique(PopIDs)[i]
          to <- min(which(PopIDs==u.ID.no))
          from <- max(which(PopIDs==u.ID.no))
      hold.names[to:from] = paste(pop.groups[i], hold.names[to:from], sep=".")
    }
    NameExtract <- hold.names
    }
    
    #temp2$Pop ## can be used as internal check of code
    
    ## get the nubmer of indivudals within each "Pop" grouping
    PopLengths <- table(temp2$Pop)
    
    ## Need to be able to tell what row each individual is in, and what population it is
    ind.vector = c(1:nrow(temp)) ### make a vector that is the number of individuals
    ind.matrix = data.frame(temp2$Pop, ind.vector) ## add populatuions to that
    
    
    ### make a dataframe to hold the number of mother and father - DF will be compoised of NULLs
    ## need to make the dataframe as long as the longer of the number of mothrs or fathers
    if(Nmothers > Nfathers){
      Nfamilies = matrix(data = rep(NA), nrow = Nmothers, ncol = 2, dimnames = NULL)
      colnames(Nfamilies) <- c("Moth", "Fath")
      Nfamilies = as.data.frame(Nfamilies)
    }
    if(Nmothers < Nfathers){
      Nfamilies = matrix(data = rep(NA), nrow = Nfathers, ncol = 2)
      colnames(Nfamilies) <- c("Moth", "Fath")
      Nfamilies = as.data.frame(Nfamilies)
    }
    if(Nmothers == Nfathers){
      Nfamilies = matrix(data = rep(NA), nrow = Nmothers, ncol = 2)
      colnames(Nfamilies) <- c("Moth", "Fath")
      Nfamilies = as.data.frame(Nfamilies)
    }
    
    ### Fill in the dataframe with the number of mother and father - if one is shorter than the other 
      ## the number will repeat 
    Nfamilies$Moth = rep(1:Nmothers, length = nrow(Nfamilies))
    Nfamilies$Fath = rep(1:Nfathers, length = nrow(Nfamilies))
    Nfamilies$Family = c(1:nrow(Nfamilies)) ## add a column that is family numbers
    
    ### How many populations are there? - need to get a random number of fathers and mothers from each
        ## also will show what numbers can choose between to draw a random sample
        ## i.e. Population 1 is between 1 and 20, because the first individial of Pop2 is at position 21
    make.random.groups = (ind.matrix$ind.vector[!duplicated(ind.matrix$temp2.Pop)])
    make.random.groups <- c(make.random.groups, (max(ind.matrix)+1)) ### add the last individual to the list 
      ## so the numbering for the last Population is correct

    ## Will use Sample - which draws between a min and max 
    ## will use the numbers in the vector describing the populations to get the first and last indv. in each pop
    ## LOOP which gets random rows within each population to use get parents from
    min.max.out <- NULL 
    for(i in 1:(length(make.random.groups)-1)){
          a = make.random.groups[i] ## from the first position in population i
          b = make.random.groups[i+1]-1 ## to the first position in population (i+1) - one pisition
          min.max = c(a, b) ## just mash those together
          min.max.out = rbind(min.max.out, min.max) ## let them marinate
      
    } ## END of random row generation loop
    
    ## Now - have the first and last indiv. in each pop - so chose a random sample of indv. between them
    
    how.many.bang <- Nmothers+Nfathers ## how many parents in each population
    how.many.pops <- length(PopLengths) ## how many populations to create crosses for - recall based on POP
    
    ## LOOP which creates Random mating groups for each POP
    family.list <- NULL # create NULL object to be populated in loop
    for(j in 1:how.many.pops){  ## repeat loop for each population
      
      family.hold <- NULL ## create NULL first time through loop, clears in subsequent looping
  
      from <- min(min.max.out[j,]) ## get the min number on the jth row 
      to <- max(min.max.out[j,]) ## get the max 
      parent.hold <- sample(to:from, how.many.bang) ## randomly sample "how.many.bang" individuals
                ## for each population from the original data between with the start and end of each        
                ## population defined by to and from
    
    ## For simplicty - the mothers are the first Nmother individuals identified at random
    ### get the genotypes of the mothers and fathers from temp2
    mother.get <- temp2[parent.hold[1:Nmothers],]
    father.get <- temp2[parent.hold[((Nmothers+1):length(parent.hold))],]
    
    ### when collected from temp2 - has "Pop" on end - remove
    mother.get <- mother.get[-which(names(mother.get) == "Pop")]
    father.get <- father.get[-which(names(father.get) == "Pop")]
    
    family.hold <- rbind(mother.get, father.get)
    rnames <- paste(NameExtract[as.numeric(rownames(family.hold ))], as.numeric(rownames(family.hold )), sep=".")
    rownames(family.hold) <- rnames
    ## iteratively name each population and asign to the global environment
    assign(x = paste("Family.out", j, sep = "."), value = family.hold, envir = globalenv())
    family.list <- rbind(family.list, paste("Family.out", j, sep = ".")) ## keep list of families to allow
              ## searching and retrieval later
    } # END of loop creating random mating groups


## Do the Pure Crosses    
#uPop.names = unique(NameExtract) ## get the number unique IDs (or pop.groups where given by user)
get.cross.names = NULL
allele.hold <- NULL
pure.prog <- txtProgressBar(min=0, max = length(family.list), style = 3) ## creates a progress bar 
  ## to give the user an indcation of the function status
print("Pure Cross Progress") ## show the user what cross the status bar refers to
for(h in 1:length(family.list)){ ## for each population run the function

family.make = get(family.list[h]) ## get the 'parents' in the hth populatioin
#print(paste0("Pure cross progress = ", ((h/length(family.list))*100)))  

for(i in 1:max(Nfamilies$Family)){ ## get families within population
  
  get.mom.who.bangs <- Nfamilies$Moth[i] ## get the number of the ith mother (number of parents can vary so is not necessarily = i)
  get.dad.who.bangs <- Nfamilies$Fath[i] ## get nubmer of ith father
  
  mom.who.bangs <- family.make[get.mom.who.bangs,] ## get genotype of ith mother
  dad.who.bangs <- family.make[(max(Nfamilies$Moth)+get.dad.who.bangs),] ## dataframe has X mothers 
  ## on top of Y fathers, so ith father is in the max(X) plus number of ith father row
  
  #mom.who.bangs <- cbind(NA, mom.who.bangs)
  #dad.who.bangs <- cbind(NA, dad.who.bangs)
  
  off.name.hold = paste(rownames(family.make[get.mom.who.bangs,]), rownames(family.make[(max(Nfamilies$Moth)+get.dad.who.bangs),]), i, sep = ".") ## this will hold the name of the mother and father so that can check offspring assignment against
  ## Can replace 50 with number of offsrping later
  off.gens.hold = data.frame(matrix(vector(), family.size, (length(ColumnData.Dup)+1))) ### have to create a data frame with the same number of columns as there are alleles, + 1 for the ID, and number of 
                                                            ### rows as there are simulated offpsring 
                                                              ### because it this data frame is specified here, it will thus blank it each time the loop goes through parent pairs
## will make dataframe that is 2 columns (one for each allele) by number of loci rows for mother and father  
mom.matrix <- data.frame(matrix(vector(), length(mom.who.bangs)/2, 2)) ## blank matrix of correct size
mom.matrix$X1 = t(mom.who.bangs[c(T,F)]) ## column of matrix is first allele
mom.matrix$X2 = t(mom.who.bangs[c(F,T)]) ## secomnd column is second allele

## Repeat for father
dad.matrix <- data.frame(matrix(vector(), length(dad.who.bangs)/2, 2))
dad.matrix$X1 <- t(dad.who.bangs[c(T,F)])
dad.matrix$X2 <- t(dad.who.bangs[c(T,F)])
## get 50 offpsring for each parent pair - can change number later if desire
    for(j in 1:family.size){ ### for each parent pair, want to make a number of simulated offpsring. 
    
dad.out <- apply(t(dad.matrix), 2, sample, 1) ## transpose the matrix, then choose one row at random from each column - simulates meiosis with no linkage
mom.out <- apply(t(mom.matrix), 2, sample, 1) ## repeat for father
      
      off.name = paste(off.name.hold, j, sep = "_") ## generates an individual offspring name composed of its parents names (off.name.hold), and the number of the offspring
      #allele.hold = data.frame(allele.hold) ### the allele.hold object needs to be coerced into a data frame format to work with it optimally
      colnames(allele.hold) = NULL ## remove the column names - makes it simpler to work with
    
    allele.hold<- data.frame(c(rbind(mom.out, dad.out))) ### binds alternating 
    off.gens.hold[j,2:length(off.gens.hold)] = t(allele.hold) ## there was an issue with trying to just add the offspring name to the front of the allele.hold object, then adding
                                                              ### this to the off.gens.hold data frame - for some reason, it would only add the level of the 'factor' R was
                                                                ## interpreting the name as - so this starts the genotypes at the second position where they should be, and keeps
                                                                  ## the first blank for now
      off.gens.hold[j,1] = off.name ## this now adds the offspring name to the front row
      colnames(off.gens.hold) = c("ID", ColumnData.Dup)
      allele.hold = NULL ## nulls out the allele.hold object so it can be recycled back through the loop
 
      #print(paste0("J progress = ", ((j/50)*100)))  
      } # J Closure
    gen.dat.file = paste("Pure", off.name.hold,  sep = ".") ## makes an object name consisting of genotypes + the two parent names
    assign(x = gen.dat.file, value = off.gens.hold, envir = globalenv()) ## assigns the object name to the generated simulated offspring genotypes, and then assigns it to the global env.
  get.cross.names = rbind(get.cross.names, gen.dat.file)
} # I Closure
setTxtProgressBar(pure.prog, h)
#print(paste0("I progress = ", ((i/max(Nfamilies$Family))*100)))  
} # H Closure


pop.fam.reps = length(get.cross.names)/how.many.pops
#length(get.cross.names)

#l = 1
#test = get(get.cross.names[1:5])

### since have made N families (greater of number of fathers or mothers entered by user) for each population
  ## will randomly sample 50 offspring from among all families by first cbinding all families into one large dataframe

#l=1
## LOOP that generatest he PURE CROSS RESULTS
get.pure.names <- NULL
for(l in 1:length(pop.groups)){ ## LOOP for the number of populations
  get.string = paste0("Pure.", pop.groups[l]) ## make a string of the name common part of the lth family name
  out.string <- paste0("Pure_", pop.groups[l])
  to.get = which(str_detect(ls(), pattern = get.string) == TRUE) ## which of the names of everything in the global
    ## environment match the names of the lth family? - vector of numbers
  objects = ls() # make a vector of the names of everythign in the Globalenv
  geno.list = list(mget(objects[c(min(to.get):max(to.get))])) ## get a list that is extracted every
    ## Family in the lth population - this is a list of lists
   geno.unlist = unlist(geno.list, recursive = FALSE) ## unlist, so is list of DF
   geno.rbind = do.call(rbind, geno.unlist) ## rbind all dataframes into one
  
   ## ranomly sample 50 individuals amongst all 
   pure.geno.out <- sample_n(tbl = geno.rbind, size = family.size, replace = FALSE)
  ## assign the 50 indv to the global env - this is the OUTPUT/RESULTS
  assign(x = out.string, value = pure.geno.out, envir = globalenv())
  as.character(get.string)
  get.pure.names = c(get.pure.names, out.string) ## get list of naems of Pure Cross results to allow them to be sourced
} ### END PURE CROSS RESULTS LOOP 


### Now make the F1

## make the list of Pure cross names a dataframe
get.pure.names.df <- data.frame(get.pure.names)

# make a NULL object to populate in the loop
F1.fam.name.list <- NULL
for(r in 1:nrow(get.pure.names.df)){ ## START of loop sample random indiv. from each pure cross to use as F1 parents
                                      
  to.get <- as.character(get.pure.names.df[r,]) ## get the name of the rth pure cross family
  F1.offspring <- get(to.get) ## hold the dataframe for use
  
  F1.parent.hold <- sample_n(tbl = F1.offspring, size = how.many.bang, replace = FALSE) ## randomly sample 
            ## X fathers + Y mothers from among the pure cross genotypes
  
  F1.fam.name <- paste(to.get, "parents", sep = "_") ## make a name for held parents
  assign(x = F1.fam.name, value = F1.parent.hold, envir = globalenv()) ## assign the held parents to the global env
  F1.fam.name.list <- c(F1.fam.name.list, F1.fam.name) ## keep a list of names to allow retrieval
  
}## END of loop to get F1 parents

## Now must get all unique combinations of Pure cross families - i.e. create unidirectional crosses amongst
## all populations
F1.par.combs <- combn(F1.fam.name.list, 2, paste, collapse = '+') ## get the combinations, combine with a "+" so can easily subset later
F1.par.combs <- data.frame(do.call(rbind, str_split(F1.par.combs, "\\+"))) ## use the "+" to split into a datframe



#if(nrow(F1.par.combs)<2){
#  F1.m.f.rev <- data.frame(F1.par.combs[,2], F1.par.combs[,1]) ## reverse the order of the PCs
#  names(F1.m.f.rev) <- names(F1.par.combs) ## give the reversed pops teh same name to allow rbind
# F1.par.combs <-  rbind(F1.par.combs, F1.m.f.rev) ## rbind orignal and reversed
  
#}# END of if statement

## Randomly sample "parents" from the pure crosses to be used in the F1
## loop is essentially the same as the one used to get the parents for the PC

j=1 ## used for internal code checking

 F1.family.list <- NULL 
  
    for(j in 1:nrow(F1.par.combs)){  ## Start of loop to get F1 parents
      family.hold <- NULL
    
    mother.get.word.F1 <- as.character(F1.par.combs[j,1]) # mother comes from the pure cross indv in the data framed named in jrth row of the first column
    father.get.word.F1 <- as.character(F1.par.combs[j,2]) # father comes from the pure ccross indv in the data framed named in jrth row of the second column
    
   # mother.pop.name.F1 <- str_extract(string = mother.get.word.F1, pattern = "_..._") ## used in alternate naming version
   # father.pop.name.F1 <- str_extract(string = father.get.word.F1, pattern = "_..._")
    mother.pop.name.F1 <- str_extract(string = mother.get.word.F1, pattern = "Pure_[A-z]{3}")
    father.pop.name.F1 <- str_extract(string = father.get.word.F1, pattern = "Pure_[A-z]{3}")
    mother.pop.name.F1 <- paste("F1par", mother.pop.name.F1, sep="_")
    
    mother.pop.F1 <- get(mother.get.word.F1) ## retrieve the population from the globalEnv
    father.pop.F1 <- get(father.get.word.F1)
    
    mother.get.F1 <- sample_n(tbl = mother.pop.F1, size = Nmothers) ## randomly select Nmothers from the mothers populaiton
    father.get.F1 <- sample_n(tbl = father.pop.F1, size = Nfathers) ## do the same for the fathers
    
    
    family.hold.F1 <- rbind(mother.get.F1, father.get.F1) ## bind the mothers above the fathers
    
    Fam.name.F1 <- paste(mother.pop.name.F1, father.pop.name.F1,  sep = ".") ## add the names of the mother and father pops so can track the genetics
    
    assign(x = Fam.name.F1, value = family.hold.F1, envir = globalenv())
    F1.family.list <- c(F1.family.list, Fam.name.F1)
      
    }
 

## MAKE THE F1 - the process is the same as for the Pure Crosses 
 F1.get.cross.names <- NULL
 F1.prog <- txtProgressBar(min=0, max = length(F1.family.list), style = 3)
 print("F1 Cross Progress")
 for(h in 1:length(F1.family.list)){

F1.family.make = get(F1.family.list[h])
#print(paste0("F1 cross progress = ", ((h/length(F1.family.list))*100)))  
setTxtProgressBar(F1.prog, h)
for(i in 1:max(Nfamilies$Family)){ ## get families within population
  
  get.mom.who.bangs <- Nfamilies$Moth[i]
  get.dad.who.bangs <- Nfamilies$Fath[i]
  
  mom.who.bangs <- F1.family.make[get.mom.who.bangs,-1] ## have to remove the first column because it is in as the population not a genotype
  dad.who.bangs <- F1.family.make[(max(Nfamilies$Moth)+get.dad.who.bangs),-1]
  
  #mom.who.bangs <- cbind(NA, mom.who.bangs)
  #dad.who.bangs <- cbind(NA, dad.who.bangs)
  
  off.name.hold = paste(rownames(F1.family.make[get.mom.who.bangs,]), ".zF1z.", gsub(x = rownames(F1.family.make[(max(Nfamilies$Moth)+get.dad.who.bangs),]), pattern = "Pure_", replacement = ""), i, sep = "X") ## this will hold the name of the mother and father so that can check offspring assignment against
  ## Can replace 50 with number of offsrping later
  off.gens.hold = data.frame(matrix(vector(), family.size, (length(ColumnData.Dup)+1))) ### have to create a data frame with the same number of columns as there are alleles, + 1 for the ID, and number of 
                                                            ### rows as there are simulated offpsring 
                                                                ### because it this data frame is specified here, it will thus blank it each time the loop goes through parent pairs
  
    mom.matrix <- data.frame(matrix(vector(), length(mom.who.bangs)/2, 2))
mom.matrix$X1 = t(mom.who.bangs[c(T,F)])
mom.matrix$X2 = t(mom.who.bangs[c(F,T)])


dad.matrix <- data.frame(matrix(vector(), length(dad.who.bangs)/2, 2))
dad.matrix$X1 <- t(dad.who.bangs[c(T,F)])
dad.matrix$X2 <- t(dad.who.bangs[c(T,F)])
    for(j in 1:family.size){ ### for each parent pair, want to make a number of simulated offpsring. 
    
dad.out <- apply(t(dad.matrix), 2, sample, 1)
mom.out <- apply(t(mom.matrix), 2, sample, 1)
      
      off.name = paste(off.name.hold, j, sep = "_") ## generates an individual offspring name composed of its parents names (off.name.hold), and the number of the offspring
      #allele.hold = data.frame(allele.hold) ### the allele.hold object needs to be coerced into a data frame format to work with it optimally
      colnames(allele.hold) = NULL ## remove the column names - makes it simpler to work with
    
    allele.hold<- data.frame(c(rbind(mom.out, dad.out)))
    off.gens.hold[j,2:length(off.gens.hold)] = t(allele.hold) ## there was an issue with trying to just add the offspring name to the front of the allele.hold object, then adding
                                                              ### this to the off.gens.hold data frame - for some reason, it would only add the level of the 'factor' R was
                                                                ## interpreting the name as - so this starts the genotypes at the second position where they should be, and keeps
                                                                  ## the first blank for now
      off.gens.hold[j,1] = off.name ## this now adds the offspring name to the front row
      colnames(off.gens.hold) = c("ID", ColumnData.Dup)
      allele.hold = NULL ## nulls out the allele.hold object so it can be recycled back through the loop
 
      } # J Closure
  gen.fam <- paste0("A1_fam_", h,  "_")
  gen.dat.file <- gsub(x =  off.name.hold, pattern = "Pure.", replacement = gen.fam)
    
    assign(x = gen.dat.file, value = off.gens.hold, envir = globalenv()) ## assigns the object name to the generated simulated offspring genotypes, and then assigns it to the global env.
  F1.get.cross.names = rbind(F1.get.cross.names, gen.dat.file)
} # I Closure

#print(paste0("I progress = ", ((i/max(Nfamilies$Family))*100)))  
} # H Closure
 

# length(gen.dat.file[1]) ### used for internal error checking
 
 
 
l=1 ## used for internal error checking
 get.F1.names = NULL
for(l in 1:length(F1.family.list)){
  get.string.F1 = paste0("A1_fam_", l,  "_") ## Cerate a string to search for based on the known beginning of the F1 cross generating function
  to.get.F1 = which(str_detect(ls(), pattern = get.string.F1) == TRUE) ## make a vector of numbers which are the positions in the globalEnv where
                                                                        ## the dataframes which match the string are
  objects.F1 = ls() ## make a vector that is the names in the globalEmv
  geno.list.F1 = list(mget(objects.F1[c(min(to.get.F1):max(to.get.F1))])) ## get all those datafreames from the globalEnv that are in the positioins 
   geno.unlist.F1 = unlist(geno.list.F1, recursive = FALSE)
   geno.rbind.F1 = do.call(rbind, geno.unlist.F1)
  
   geno.out.F1 <- sample_n(tbl = geno.rbind.F1, size = family.size, replace = FALSE)
  
   ## remove the extraneous from the name vector
   nam.geno.hold <- names(geno.list.F1[[1]])[1]
  #nam.geno.hold <- str_split(string = nam.geno.hold, pattern = get.string.F1)[[1]][2]
  nam.geno.hold<- gsub(pattern = "[0-9]", replacement = "", x = nam.geno.hold)
  nam.geno.hold<- gsub(pattern = "X", replacement = "", x = nam.geno.hold)
  nam.geno.hold <- gsub(pattern = "\\.{2,}", replacement = ".", nam.geno.hold)
  nam.geno.hold <- gsub(pattern = "fam__", replacement = "", nam.geno.hold)
  nam.geno.hold <- gsub(pattern = "\\.A_", replacement = ".", nam.geno.hold)
  nam.geno.hold <- gsub(pattern = "A_", replacement = "F1_", nam.geno.hold)
  
  assign(x = nam.geno.hold, value = geno.out.F1, envir = globalenv())
  get.F1.names = c(get.F1.names, nam.geno.hold)
} 
 
 
 
 #rm(F1.get.cross.names)

F2.par.combs <- combn(get.F1.names, 2, paste, collapse = '+') ### get all combinations of the families created in teh F1
F2.par.combs <- data.frame(do.call(rbind, str_split(F2.par.combs, "\\+"))) 

if(length(get.F1.names)==1){
  F2.par.combs <- data.frame(cbind(get.F1.names, get.F1.names))
}

q=1 ## used for internal code checking
F2.family.list <- NULL

for(q in 1:nrow(F2.par.combs)){  
      family.hold <- NULL
    
    mother.get.word.F2 <- as.character(F2.par.combs[q,1])
    father.get.word.F2 <- as.character(F2.par.combs[q,2])
    
    mother.pop.name.F2 <- unlist(str_split(string = mother.get.word.F2, pattern = "F1_"))[2]
    father.pop.name.F2 <-  unlist(str_split(string = father.get.word.F2, pattern = "F1_"))[2]
    
    mother.pop.F2 <- get(mother.get.word.F2) 
    father.pop.F2 <- get(father.get.word.F2)
    
    mother.get.F2 <- sample_n(tbl = mother.pop.F2, size = Nmothers)
    father.get.F2 <- sample_n(tbl = father.pop.F2, size = Nfathers)
    
    mother.pop.name.F2 <- gsub(x = mother.pop.name.F2, pattern = "Pure_", replacement = "")
    father.pop.name.F2 <- gsub(x = father.pop.name.F2, pattern = "Pure_", replacement = "")
    
    
    family.hold.F2 <- rbind(mother.get.F2, father.get.F2)
    
    Fam.name.F2 <- paste("F2.", mother.pop.name.F2, "x.", father.pop.name.F2,  sep = "")
    
    assign(x = Fam.name.F2, value = family.hold.F2, envir = globalenv())
    F2.family.list <- rbind(F2.family.list, Fam.name.F2)
      
}

#rownames(F2_Aqu.WLDxWLD.Aqu)

### Make the F2 Families

F2.get.cross.names <- NULL
F2.prog <- txtProgressBar(min = 0, max = length(F2.family.list), style = 3)
 print("F2 Cross Progress")
 for(h in 1:length(F2.family.list)){

F2.family.make = get(F2.family.list[h])
#print(paste0("F2 cross progress = ", ((h/length(F2.family.list))*100)))  
setTxtProgressBar(F2.prog, h)
for(i in 1:max(Nfamilies$Family)){ ## get families within population
  
  get.mom.who.bangs <- Nfamilies$Moth[i]
  get.dad.who.bangs <- Nfamilies$Fath[i]
  
  mom.who.bangs <- F2.family.make[get.mom.who.bangs,-1]
  dad.who.bangs <- F2.family.make[(max(Nfamilies$Moth)+get.dad.who.bangs),-1]
  
  #mom.who.bangs <- cbind(NA, mom.who.bangs)
  #dad.who.bangs <- cbind(NA, dad.who.bangs)
  
  off.name.hold = paste(rownames(F2.family.make[get.mom.who.bangs,]), ".bF2b.", gsub(x = rownames(F2.family.make[(max(Nfamilies$Moth)+get.dad.who.bangs),]), pattern = "F1.", replacement = ""), i, sep = "X")
  ## Can replace 50 with number of offsrping later
  off.gens.hold = data.frame(matrix(vector(), family.size, (length(ColumnData.Dup)+1))) ### have to create a data frame with the same number of columns as there are alleles, + 1 for the ID, and number of 
                                                             ### because it this data frame is specified here, it will thus blank it each time the loop goes through parent pairs
  
    mom.matrix <- data.frame(matrix(vector(), length(mom.who.bangs)/2, 2))
mom.matrix$X1 = t(mom.who.bangs[c(T,F)])
mom.matrix$X2 = t(mom.who.bangs[c(F,T)])


dad.matrix <- data.frame(matrix(vector(), length(dad.who.bangs)/2, 2))
dad.matrix$X1 <- t(dad.who.bangs[c(T,F)])
dad.matrix$X2 <- t(dad.who.bangs[c(T,F)])
    for(j in 1:family.size){ ### for each parent pair, want to make a number of simulated offpsring. 
    
dad.out <- apply(t(dad.matrix), 2, sample, 1)
mom.out <- apply(t(mom.matrix), 2, sample, 1)
      
      off.name = paste(off.name.hold, j, sep = "_") ## generates an individual offspring name composed of its parents names (off.name.hold), and the number of the offspring
      #allele.hold = data.frame(allele.hold) ### the allele.hold object needs to be coerced into a data frame format to work with it optimally
      colnames(allele.hold) = NULL ## remove the column names - makes it simpler to work with
    
    allele.hold<- data.frame(c(rbind(mom.out, dad.out)))
    off.gens.hold[j,2:length(off.gens.hold)] = t(allele.hold) ## there was an issue with trying to just add the offspring name to the front of the allele.hold object, then adding
                                                              ### this to the off.gens.hold data frame - for some reason, it would only add the level of the 'factor' R was
                                                                ## interpreting the name as - so this starts the genotypes at the second position where they should be, and keeps
                                                                  ## the first blank for now
      off.gens.hold[j,1] = off.name ## this now adds the offspring name to the front row
      colnames(off.gens.hold) = c("ID", ColumnData.Dup)
      allele.hold = NULL ## nulls out the allele.hold object so it can be recycled back through the loop
 
      } # J Closure
    gen.fam <- paste0("B2_fam_", h,  "_")
  gen.dat.file <- gsub(x =  off.name.hold, pattern = "A1_", replacement = gen.fam)
   
 
    assign(x = gen.dat.file, value = off.gens.hold, envir = globalenv()) ## assigns the object name to the generated simulated offspring genotypes, and then assigns it to the global env.
  F2.get.cross.names = rbind(F2.get.cross.names, gen.dat.file)
} # I Closure

#print(paste0("I progress = ", ((i/max(Nfamilies$Family))*100)))  
} # H Closure
 
 
 
 
l=1 ## used for internal error checking
 ## loop is essentially the same as the the one for the F1 families - just modified to fit the name structure of the F2
 get.F2.names = NULL
for(l in 1:length(F2.family.list)){
  get.string.F2 = paste0("B2_fam_", l, "_", "fam_[0-9]{1,}_")
  to.get.F2 = which(str_detect(ls(), pattern = get.string.F2) == TRUE)
  objects.F2 = ls()
  geno.list.F2 = list(mget(objects.F2[c(min(to.get.F2):max(to.get.F2))]))
   geno.unlist.F2 = unlist(geno.list.F2, recursive = FALSE)
   geno.rbind.F2 = do.call(rbind, geno.unlist.F2)
  
   
   geno.out.F2 <- sample_n(tbl = geno.rbind.F2, size = family.size, replace = FALSE)
   
   nam.geno.hold <- names(geno.list.F2[[1]])[1]
  #nam.geno.hold <- str_split(string = nam.geno.hold, pattern = get.string.F2)[[1]][2]
  nam.geno.hold<- gsub(pattern = "[0-9]", replacement = "", x = nam.geno.hold)
  nam.geno.hold<- gsub(pattern = "X", replacement = "", x = nam.geno.hold)
  nam.geno.hold<- gsub(pattern = "fam", replacement = "", x = nam.geno.hold)
  nam.geno.hold <- gsub(pattern = "\\_{,}", replacement = ".", nam.geno.hold)
  nam.geno.hold <- gsub(pattern = "\\.{2,}", replacement = ".", nam.geno.hold)
  nam.geno.hold <- gsub(pattern = "B.", replacement = "", nam.geno.hold)
  nam.geno.hold <- paste("F2", nam.geno.hold, sep="_")

  assign(x = nam.geno.hold, value = geno.out.F2, envir = globalenv())
  get.F2.names = rbind(get.F2.names, nam.geno.hold)
} 
  
 
 BC.names.list <- c(get.pure.names, get.F1.names)
 
BC.par.combs <- combn(BC.names.list, 2, paste, collapse = '+')
BC.par.combs <- data.frame(do.call(rbind, str_split(BC.par.combs, "\\+"))) 

## when creating all combinations, also creates one which are actually pure crosses and F2 crosses - remove them

if(any(str_detect(string = BC.par.combs$X1, pattern = "F1_")==TRUE)){
BC.par.combs <- BC.par.combs[-which(str_detect(string = BC.par.combs$X1, pattern = "F1_")==TRUE),] ## This causes it to fail when there are none to remove - why?????
}
if(any(str_detect(string = BC.par.combs$X2, pattern = "Pure_")==TRUE)){
BC.par.combs <- BC.par.combs[-which(str_locate(string = BC.par.combs$X2, pattern = "Pure_")==TRUE),]
}

BC.par.combs <- droplevels(BC.par.combs)

BC.family.list <- NULL

for(q in 1:nrow(BC.par.combs)){  
      family.hold <- NULL
    
    mother.get.word.BC <- as.character(BC.par.combs[q,1])
    father.get.word.BC <- as.character(BC.par.combs[q,2])
    
    mother.pop.name.BC <- unlist(str_split(string = mother.get.word.BC, pattern = "Pure_"))[2] ## 
    father.pop.name.BC <-  unlist(str_split(string = father.get.word.BC, pattern = "F1_"))[2]
    
    mother.pop.BC <- get(mother.get.word.BC) 
    father.pop.BC <- get(father.get.word.BC)
    
    mother.get.BC <- sample_n(tbl = mother.pop.BC, size = Nmothers)
    father.get.BC <- sample_n(tbl = father.pop.BC, size = Nfathers)
    
    
    family.hold.BC <- rbind(mother.get.BC, father.get.BC)
    
    Fam.name.BC <- paste("BC.pure", mother.pop.name.BC, sep = ".")
    Fam.name.BC <- paste(Fam.name.BC, father.pop.name.BC,  sep = "x")
    
    
    assign(x = Fam.name.BC, value = family.hold.BC, envir = globalenv())
    BC.family.list <- c(BC.family.list, Fam.name.BC)
      
    } 
 
i=1

BC.get.cross.names <- NULL
BC.prog <- txtProgressBar(min = 0, max = length(BC.family.list), style = 3)
 print("Back Cross Progress")
 for(h in 1:length(BC.family.list)){

BC.family.make = get(BC.family.list[h])
#print(paste0("BC cross progress = ", ((h/length(BC.family.list))*100)))  
setTxtProgressBar(BC.prog, h)
for(i in 1:max(Nfamilies$Family)){ ## get families within population
  
  get.mom.who.bangs <- Nfamilies$Moth[i]
  get.dad.who.bangs <- Nfamilies$Fath[i]
  
  mom.who.bangs <- BC.family.make[get.mom.who.bangs,-1]
  dad.who.bangs <- BC.family.make[(max(Nfamilies$Moth)+get.dad.who.bangs),-1]
  
  #mom.who.bangs <- cbind(NA, mom.who.bangs)
  #dad.who.bangs <- cbind(NA, dad.who.bangs)
  
  off.name.hold = paste(rownames(BC.family.make[get.mom.who.bangs,]), "cBCc", rownames(BC.family.make[(max(Nfamilies$Moth)+get.dad.who.bangs),]), i, sep = "X")
  ## Can replace 50 with number of offsrping later
  off.gens.hold = data.frame(matrix(vector(), family.size, (length(ColumnData.Dup)+1))) ### have to create a data frame with the same number of columns as there are alleles, + 1 for the ID, and number of 
                                                                ### because it this data frame is specified here, it will thus blank it each time the loop goes through parent pairs
  
    mom.matrix <- data.frame(matrix(vector(), length(mom.who.bangs)/2, 2))
mom.matrix$X1 = t(mom.who.bangs[c(T,F)])
mom.matrix$X2 = t(mom.who.bangs[c(F,T)])


dad.matrix <- data.frame(matrix(vector(), length(dad.who.bangs)/2, 2))
dad.matrix$X1 <- t(dad.who.bangs[c(T,F)])
dad.matrix$X2 <- t(dad.who.bangs[c(T,F)])
    for(j in 1:family.size){ ### for each parent pair, want to make a number of simulated offpsring. 
    
dad.out <- apply(t(dad.matrix), 2, sample, 1)
mom.out <- apply(t(mom.matrix), 2, sample, 1)
      
      off.name = paste(off.name.hold, j, sep = "_") ## generates an individual offspring name composed of its parents names (off.name.hold), and the number of the offspring
      #allele.hold = data.frame(allele.hold) ### the allele.hold object needs to be coerced into a data frame format to work with it optimally
      colnames(allele.hold) = NULL ## remove the column names - makes it simpler to work with
    
    allele.hold<- data.frame(c(rbind(mom.out, dad.out)))
    off.gens.hold[j,2:length(off.gens.hold)] = t(allele.hold) ## there was an issue with trying to just add the offspring name to the front of the allele.hold object, then adding
                                                              ### this to the off.gens.hold data frame - for some reason, it would only add the level of the 'factor' R was
                                                                ## interpreting the name as - so this starts the genotypes at the second position where they should be, and keeps
                                                                  ## the first blank for now
      off.gens.hold[j,1] = off.name ## this now adds the offspring name to the front row
      colnames(off.gens.hold) = c("ID", ColumnData.Dup)
      allele.hold = NULL ## nulls out the allele.hold object so it can be recycled back through the loop
 
      } # J Closure
  
  #gen.dat.file <- gsub(pattern = "")
    gen.fam <- paste0("CC_fam_", h,  "_")
  gen.dat.file <- gsub(x =  off.name.hold, pattern = "Pure", replacement = gen.fam)
 
    assign(x = gen.dat.file, value = off.gens.hold, envir = globalenv()) ## assigns the object name to the generated simulated offspring genotypes, and then assigns it to the global env.
  BC.get.cross.names = rbind(BC.get.cross.names, gen.dat.file)
} # I Closure

#print(paste0("I progress = ", ((i/max(Nfamilies$Family))*100)))  
} # H Closure
 

l=1
get.BC.names = NULL
for(l in 1:length(BC.family.list)){
  get.string.BC = (BC.family.list[l])
  to.get.BC = which(str_detect(ls(), pattern = paste0("CC_fam_",l,"_")) == TRUE)
  objects.BC = ls()
  geno.list.BC = list(mget(objects.BC[c(min(to.get.BC):max(to.get.BC))]))
   geno.unlist.BC = unlist(geno.list.BC, recursive = FALSE)
   geno.rbind.BC = do.call(rbind, geno.unlist.BC)
  out.string.BC <- gsub(x = get.string.BC, pattern = "BC.pure.", replacement = "BC_")
   
   geno.out.BC <- sample_n(tbl = geno.rbind.BC, size = family.size, replace = FALSE)

  assign(x = out.string.BC, value = geno.out.BC, envir = globalenv())
  get.BC.names = rbind(get.BC.names, out.string.BC)
  #get.BC.names <- gsub(pattern = "BC.", replacement = "BC_", x = get.BC.names)
} 
  

get.pure.names

get.F1.names

get.F2.names

get.BC.names


### Edit the individual names

b=1
for(b in 1:length(get.pure.names)){
  
  fam.to.convert.indv.names <- get.pure.names[b]
  
 fam.to.cv.names <- data.frame(get(fam.to.convert.indv.names))
  cv.nm <- fam.to.cv.names$ID
  
    cv.nm <- gsub(pattern = "[0,2-9]", x = cv.nm, replacement = ".")
    cv.nm <- gsub(pattern = "X", x = cv.nm, replacement = ".")
    cv.nm <- gsub(pattern = "\\.{2,}", x = cv.nm, replacement = ".")
    cv.nm <- gsub(pattern = "\\.[0-9]", x = cv.nm, replacement = "")
    cv.nm <- gsub(pattern = "z", x = cv.nm, replacement = ".")
    cv.nm <- gsub(pattern = "\\.1", x = cv.nm, replacement = ".")
    cv.nm <- gsub(pattern = "\\_1", x = cv.nm, replacement = ".")
    cv.nm <- gsub(pattern = "\\.\\_", x = cv.nm, replacement = "")
    cv.nm <- gsub(pattern = "\\_\\.", x = cv.nm, replacement = ".")
    cv.nm <- gsub(pattern = "\\.{2,}", x = cv.nm, replacement = ".")
    cv.nm <- gsub(pattern = "\\.$", x = cv.nm, replacement = "")
    cv.nm <- gsub(pattern = "[0-9]$", x = cv.nm, replacement = "")

    fam.to.cv.names$ID <- cv.nm
  
  
  assign(x = fam.to.convert.indv.names, value = fam.to.cv.names, envir = globalenv())
}

 
b=1
for(b in 1:length(get.F1.names)){
  
  fam.to.convert.indv.names <- get.F1.names[b]
  
 fam.to.cv.names <- data.frame(get(fam.to.convert.indv.names))
  cv.nm <- fam.to.cv.names$ID
  
    cv.nm <- gsub(pattern = "[0-9]", x = cv.nm, replacement = ".")
    cv.nm <- gsub(pattern = "X", x = cv.nm, replacement = ".")
    cv.nm <- gsub(pattern = "\\_", x = cv.nm, replacement = "")
    cv.nm <- gsub(pattern = "\\.{2,}", x = cv.nm, replacement = ".")
    cv.nm <- gsub(pattern = "Pure.", x = cv.nm, replacement = "")
    cv.nm <- gsub(pattern = "zF.z.", x = cv.nm, replacement = "F1.")
    cv.nm <- gsub(pattern = "\\.$", x = cv.nm, replacement = "")
    
    cv.nm <- paste("F1", cv.nm, sep="_")

    fam.to.cv.names$ID <- cv.nm
  
  
  assign(x = fam.to.convert.indv.names, value = fam.to.cv.names, envir = globalenv())
}


for(b in 1:length(get.F2.names)){
  
  fam.to.convert.indv.names <- get.F2.names[b]
  
 fam.to.cv.names <- data.frame(get(fam.to.convert.indv.names))
  cv.nm <- fam.to.cv.names$ID
    
    cv.nm <- gsub(pattern = "[0-9]", x = cv.nm, replacement = ".")
    cv.nm <- gsub(pattern = "X", x = cv.nm, replacement = ".")
    cv.nm <- gsub(pattern = "\\_", x = cv.nm, replacement = ".")
    cv.nm <- gsub(pattern = "A\\.", x = cv.nm, replacement = ".")
    cv.nm <- gsub(pattern = "\\.{2,}", x = cv.nm, replacement = ".")
    cv.nm <- gsub(pattern = "\\.$", x = cv.nm, replacement = "")
    cv.nm <- gsub(pattern = "fam", x = cv.nm, replacement = ".")
    cv.nm <- gsub(pattern = "^\\.", x = cv.nm, replacement = "")
    cv.nm <- gsub(pattern = "zF.z", x = cv.nm, replacement = "F1")
    cv.nm <- gsub(pattern = "z", x = cv.nm, replacement = "F1")
    cv.nm <- gsub(pattern = "bF.b", x = cv.nm, replacement = "F2")
    cv.nm <- gsub(pattern = "\\.{2,}", x = cv.nm, replacement = ".")
    cv.nm <- gsub(pattern = "^\\.", x = cv.nm, replacement = "F2_")

    fam.to.cv.names$ID <- cv.nm
  
  
  assign(x = fam.to.convert.indv.names, value = fam.to.cv.names, envir = globalenv())
}

for(b in 1:length(get.BC.names)){
  
  fam.to.convert.indv.names <- get.BC.names[b]
  
  fam.to.cv.names <- get(fam.to.convert.indv.names)
  cv.nm <- fam.to.cv.names$ID
  
    cv.nm <- gsub(pattern = "[0-9]", x = cv.nm, replacement = ".")
    cv.nm <- gsub(pattern = "X", x = cv.nm, replacement = ".")
    cv.nm <- gsub(pattern = "\\_", x = cv.nm, replacement = ".")
    cv.nm <- gsub(pattern = "A\\.", x = cv.nm, replacement = ".")
    cv.nm <- gsub(pattern = "\\.{2,}", x = cv.nm, replacement = ".")
    cv.nm <- gsub(pattern = "\\.$", x = cv.nm, replacement = "")
    cv.nm <- gsub(pattern = "fam", x = cv.nm, replacement = ".")
    cv.nm <- gsub(pattern = "^\\.", x = cv.nm, replacement = "")
    cv.nm <- gsub(pattern = "zF.z", x = cv.nm, replacement = "F1")
    cv.nm <- gsub(pattern = "z", x = cv.nm, replacement = "F1")
    cv.nm <- gsub(pattern = "cBCc", x = cv.nm, replacement = "BC")
    cv.nm <- gsub(pattern = "\\.{2,}", x = cv.nm, replacement = ".")
    cv.nm <- gsub(pattern = "Pure.", x = cv.nm, replacement = "Pure_")
    fam.to.cv.names$ID <- cv.nm
  
  
  assign(x = fam.to.convert.indv.names, value = fam.to.cv.names, envir = globalenv())
}



### Edit the individual names




for(b in 1:length(get.pure.names)){
  #test.fam <- Pure_BDN
  fam.to.bind.name <- get.pure.names[b]
  #fam.to.bind.name <- "test.fam"
  fam.to.bind <- get(fam.to.bind.name)
  indiv.hold <- fam.to.bind[,1]
  loci.bind <- which(str_detect(string = colnames(fam.to.bind), pattern = "\\.2")==TRUE)
    col.out <- NULL
    for(s in 1:length(loci.bind)){
      place.1 <- (loci.bind[s]-1)
      place.2 <- loci.bind[s]
      hold.col <- paste0(fam.to.bind[,place.1], fam.to.bind[,place.2])
      col.out <- cbind(col.out, hold.col)

        }
  
  fam.reord <- cbind(indiv.hold,col.out)
  colnames(fam.reord) <- c(colnames(fam.to.bind[1]), colnames(fam.to.bind[c((loci.bind-1))]))
  assign(x = fam.to.bind.name, fam.reord, envir = globalenv())
  
}

for(b in 1:length(get.pure.names)){
  
  fam.to.remove.untyped.name <- get.pure.names[b]
  
  fam.to.remove.untyped <- get(fam.to.remove.untyped.name)
  fam.to.remove.untyped[which(str_detect(string = fam.to.remove.untyped, pattern = "000")==TRUE)] = "000000"
  assign(x = fam.to.remove.untyped.name, value = fam.to.remove.untyped, envir = globalenv())
}
  


for(b in 1:length(get.F1.names)){
  #test.fam <- Pure_BDN
  fam.to.bind.name <- get.F1.names[b]
  #fam.to.bind.name <- "test.fam"
  fam.to.bind <- get(fam.to.bind.name)
  indiv.hold <- fam.to.bind[,1]
  loci.bind <- which(str_detect(string = colnames(fam.to.bind), pattern = "\\.2")==TRUE)
    col.out <- NULL
    for(s in 1:length(loci.bind)){
      place.1 <- (loci.bind[s]-1)
      place.2 <- loci.bind[s]
      hold.col <- paste0(fam.to.bind[,place.1], fam.to.bind[,place.2])
      col.out <- cbind(col.out, hold.col)

        }
  
  fam.reord <- cbind(indiv.hold,col.out)
  colnames(fam.reord) <- c(colnames(fam.to.bind[1]), colnames(fam.to.bind[c((loci.bind-1))]))
  assign(x = fam.to.bind.name, fam.reord, envir = globalenv())
  
}

for(b in 1:length(get.F1.names)){
  
  fam.to.remove.untyped.name <- get.F1.names[b]
  
  fam.to.remove.untyped <- get(fam.to.remove.untyped.name)
  fam.to.remove.untyped[which(str_detect(string = fam.to.remove.untyped, pattern = "000")==TRUE)] = "000000"
  assign(x = fam.to.remove.untyped.name, value = fam.to.remove.untyped, envir = globalenv())
}


for(b in 1:length(get.F2.names)){
  #test.fam <- Pure_BDN
  fam.to.bind.name <- get.F2.names[b]
  #fam.to.bind.name <- "test.fam"
  fam.to.bind <- get(fam.to.bind.name)
  indiv.hold <- fam.to.bind[,1]
  loci.bind <- which(str_detect(string = colnames(fam.to.bind), pattern = "\\.2")==TRUE)
    col.out <- NULL
    for(s in 1:length(loci.bind)){
      place.1 <- (loci.bind[s]-1)
      place.2 <- loci.bind[s]
      hold.col <- paste0(fam.to.bind[,place.1], fam.to.bind[,place.2])
      col.out <- cbind(col.out, hold.col)

        }
  
  fam.reord <- cbind(indiv.hold,col.out)
  colnames(fam.reord) <- c(colnames(fam.to.bind[1]), colnames(fam.to.bind[c((loci.bind-1))]))
  assign(x = fam.to.bind.name, fam.reord, envir = globalenv())
  
}

for(b in 1:length(get.F2.names)){
  
  fam.to.remove.untyped.name <- get.F2.names[b]
  
  fam.to.remove.untyped <- get(fam.to.remove.untyped.name)
  fam.to.remove.untyped[which(str_detect(string = fam.to.remove.untyped, pattern = "000")==TRUE)] = "000000"
  assign(x = fam.to.remove.untyped.name, value = fam.to.remove.untyped, envir = globalenv())
}

b=1
for(b in 1:length(get.BC.names)){
  #test.fam <- Pure_BDN
  fam.to.bind.name <- get.BC.names[b]
  #fam.to.bind.name <- "test.fam"
  fam.to.bind <- get(fam.to.bind.name)
  indiv.hold <- fam.to.bind[,1]
  loci.bind <- which(str_detect(string = colnames(fam.to.bind), pattern = "\\.2")==TRUE)
    col.out <- NULL
    for(s in 1:length(loci.bind)){
      place.1 <- (loci.bind[s]-1)
      place.2 <- loci.bind[s]
      hold.col <- paste0(fam.to.bind[,place.1], fam.to.bind[,place.2])
      col.out <- cbind(col.out, hold.col)

        }
  
  fam.reord <- cbind(indiv.hold,col.out)
  colnames(fam.reord) <- c(colnames(fam.to.bind[1]), colnames(fam.to.bind[c((loci.bind-1))]))
  assign(x = fam.to.bind.name, fam.reord, envir = globalenv())
  
}

for(b in 1:length(get.BC.names)){
  
  fam.to.remove.untyped.name <- get.BC.names[b]
  
  fam.to.remove.untyped <- get(fam.to.remove.untyped.name)
  fam.to.remove.untyped[which(str_detect(string = fam.to.remove.untyped, pattern = "000")==TRUE)] = "000000"
  assign(x = fam.to.remove.untyped.name, value = fam.to.remove.untyped, envir = globalenv())
}




#Now recompile the GenePop format


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
    
    #the number of individuals for all popualtions but the last (Pop tagged to the end)
    
p.names.get.nrow <- get(get.pure.names[1])
F1.names.get.nrow <- get(get.F1.names[1])
F2.names.get.nrow <- get(get.F2.names[1])
BC.names.get.nrow <- get(get.BC.names[1])

## n.pops for pure already defined
n.pop.F1 <- nrow(data.frame(get.F1.names))
n.pop.F2 <- nrow(data.frame(get.F2.names))
n.pop.BC <- nrow(data.frame(get.BC.names))

    PopLengths.Pure <- nrow(p.names.get.nrow)
    PopLengths.F1 <- nrow(F1.names.get.nrow)
    PopLengths.F2 <- nrow(F2.names.get.nrow)
    PopLengths.BC <- nrow(BC.names.get.nrow)
    
    
    PopPosition.Pure <- NULL
    PopPosition.Pure <- c(PopLengths.Pure[1]+1,rep(NA,length(npops)-1))
      for (i in 2:length(PopPosition.Pure)){
          PopPosition.Pure[i] <- PopLengths.Pure+PopPosition.Pure[i-1]
      }
  
      
      PopPosition.F1 <- NULL    
       PopPosition.F1 <- c(PopLengths.F1[1]+1,rep(NA,(n.pop.F1)-1))
      for (i in 2:length(PopPosition.F1)){
          PopPosition.F1[i] <- PopLengths.F1+PopPosition.F1[i-1]
      }
  
     
       PopPosition.F2 <- NULL
       PopPosition.F2 <- c(PopLengths.F2[1]+1,rep(NA,(n.pop.F2)-1))
      for (i in 2:length(PopPosition.F2)){
          PopPosition.F2[i] <- PopLengths.F2+PopPosition.F2[i-1]
      }
       
       
       PopPosition.BC <- NULL
       PopPosition.BC <- c(PopLengths.BC[1]+1,rep(NA,(n.pop.BC)-1))
      for (i in 2:length(PopPosition.BC)){
          PopPosition.BC[i] <- PopLengths.BC+PopPosition.BC[i-1]
      }
 
       PopPosition.Pure <- PopPosition.Pure[-length(PopPosition.Pure)]
       PopPosition.F1 <- PopPosition.F1[-length(PopPosition.F1)]
       PopPosition.F2 <- PopPosition.F2[-length(PopPosition.F2)]
       PopPosition.BC <- PopPosition.BC[-length(PopPosition.BC)]

  Pure.out <- NULL
       for(v in 1:length(npops)){
         fam.to.bind.name <- get.pure.names[v]
          fam.to.bind <- get(fam.to.bind.name)
          Pure.out <- rbind(Pure.out, fam.to.bind)
       }
  
  
    F1.out <- NULL
       for(v in 1:n.pop.F1){
         fam.to.bind.name <- get.F1.names[v]
          fam.to.bind <- get(fam.to.bind.name)
          F1.out <- rbind(F1.out, fam.to.bind)
       }
    
     F2.out <- NULL
       for(v in 1:n.pop.F2){
         fam.to.bind.name <- get.F2.names[v]
          fam.to.bind <- get(fam.to.bind.name)
          F2.out <- rbind(F2.out, fam.to.bind)
       }

        BC.out <- NULL
       for(v in 1:n.pop.BC){
         fam.to.bind.name <- get.BC.names[v]
          fam.to.bind <- get(fam.to.bind.name)
          BC.out <- rbind(BC.out, fam.to.bind)
       }
   
        All.out <- rbind(Pure.out, F1.out, F2.out, BC.out)
    
    # paste together the Loci as one long integer seperated for each loci by a space
    Loci.pure <- do.call(paste,c(data.frame(Pure.out[,]), sep=" "))
    Loci.F1 <- do.call(paste,c(data.frame(F1.out[,]), sep=" "))
    Loci.F2 <- do.call(paste,c(data.frame(F2.out[,]), sep=" "))
    Loci.BC <- do.call(paste,c(data.frame(BC.out[,]), sep=" "))
    Loci.all <- do.call(paste,c(data.frame(All.out[,]), sep=" "))
    
        PopPosition.ALL <- NULL
        n.pop.ALL <- (max(npops) + n.pop.F1 + n.pop.F2 + n.pop.BC)
       PopPosition.ALL <- c(PopLengths.BC[1]+1,rep(NA,(n.pop.ALL)-1))
      for (i in 2:length(PopPosition.ALL)){
          PopPosition.ALL[i] <- PopLengths.BC+PopPosition.ALL[i-1]
      }
    
      
    #Insert the value of "Pop" which partitions the data among populations
    Loci.pure <- insert.vals(Vec=Loci.pure, breaks = PopPosition.Pure, newVal = "Pop")
    Loci.F1 <- insert.vals(Vec=Loci.F1, breaks = PopPosition.F1, newVal = "Pop")
    Loci.F2 <- insert.vals(Vec=Loci.F2, breaks = PopPosition.F2, newVal = "Pop")
    Loci.BC <- insert.vals(Vec=Loci.BC, breaks = PopPosition.BC, newVal = "Pop")
    Loci.all <- insert.vals(Vec=Loci.all, breaks = PopPosition.ALL, newVal = "Pop")
    
    
    #Add the first "Pop" label
    Loci.pure <- c("Pop",Loci.pure) 
    Loci.F1 <- c("Pop", Loci.F1)
    Loci.F2 <- c("Pop", Loci.F2)
    Loci.BC <- c("Pop", Loci.BC)
    Loci.all <- c("Pop", Loci.all)
    
    
    
    ## Add the column labels and the stacks version
    names.to.give <- colnames(p.names.get.nrow)
    names.to.give <- names.to.give[-1]
    names.to.give.Pure <- c("Simulated Pure Crosses", names.to.give)
    names.to.give.F1 <- c("Simulated F1 Crosses", names.to.give)
    names.to.give.F2 <- c("Simulated F2 Crosses", names.to.give)
    names.to.give.BC <- c("Simulated Back-Crosses", names.to.give)
    names.to.give.all <- c("Simulated Pure Crosses, F1, F2 and Back-Crosses", names.to.give)
    
    Loci.pure <- c(names.to.give.Pure,Loci.pure) 
    Loci.F1 <- c(names.to.give.F1, Loci.F1)
    Loci.F2 <- c(names.to.give.F2, Loci.F2)
    Loci.BC <- c(names.to.give.BC, Loci.BC)
    Loci.all <- c(names.to.give.all, Loci.all)
    
    
  
    
    #out.name <- deparse(substitute(GenePop))
    out.name.pure<- paste(out.name, "pure.crosses.txt", sep =".")
    out.name.F1 <- paste(out.name, "F1.crosses.txt", sep =".")
    out.name.F2 <- paste(out.name, "F2.crosses.txt", sep =".")
    out.name.BC <- paste(out.name, "BC.crosses.txt", sep =".")
    out.name.all <- paste(out.name, "All.crosses.txt", sep =".")
    
    

    
    # Save the file - working
    write.table(x = Loci.pure,file = out.name.pure, col.names=FALSE,row.names=FALSE)
    write.table(x = Loci.F1,file = out.name.F1, col.names=FALSE,row.names=FALSE)
    write.table(x = Loci.F2,file = out.name.F2, col.names=FALSE,row.names=FALSE)
    write.table(x = Loci.BC,file = out.name.BC, col.names=FALSE,row.names=FALSE)
    write.table(x = Loci.all,file = out.name.all, col.names=FALSE,row.names=FALSE)
    
    



rm(list=ls())














rm(list = family.list)
rm(list = get.cross.names)
rm(list = F1.fam.name.list)
rm(list = F1.get.cross.names)
rm(list = F2.family.list)
rm(list = F2.get.cross.names)
rm(list = BC.family.list)
rm(list = BC.get.cross.names)
rm(mom.who.bangs)
rm(dad.who.bangs)
rm(ind.matrix)
rm(min.max.out)
rm(allele.dad)
rm(allele.mom)
rm(allele.hold)
rm(family.list)
rm(get.cross.names)
rm(get.cross.names)
rm(get.pure.names)
rm(get.pure.names.df)
rm(BC.get.cross.names)
rm(BC.family.list)
rm(BC.par.combs)
rm(BC.names.list)
rm(get.BC.names)
rm(F1.get.cross.names)
rm(F1.family.list)
rm(F1.par.combs)
rm(F1.fam.name.list)
rm(get.F1.names)
rm(F2.get.cross.names)
rm(F2.family.list)
rm(F2.par.combs)
rm(get.F2.names)

}



do.sex(genotype = GenePopData, Nmothers = 5, Nfathers = 5)
rm(list=ls())