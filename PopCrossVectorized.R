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
rm(list=ls())

GenePopData <- read.table("SC_Pops_WildFarmed_Genepop2.txt", header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
GenePopData <- read.table("Crab_Outliers_Genepop.txt", header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)

GenePopData <- read.table("Salmon_South_Coast_noOutlier.txt", header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
GenePop <- GenePopData
Nmothers = 5
Nfathers = 5
out.name <- "SexFnTest"


#do.sex <- function(genotype, Nmothers, Nfathers, Noff, pop, F1, F2 = No, BC = Yes){


 # do.sex <- function(genotype, Nmothers, Nfathers, out.name, populations=NULL){ 
    #Check to make sure the packages required are there
 # packages <- c("dplyr", "tidyr", "stringr", "plyr") ## which packages do we need?
  if (length(setdiff(packages, rownames(installed.packages()))) > 0) { ### checks that the required packages are among those 
    ## installed by comparing the differences in the length of a vector of items included in list (i.e. is package among all installed)
    install.packages(setdiff(packages, rownames(installed.packages())))  ## if a package is not installed, insall it
  } ### this will only work if someone has the CRAN mirror set (I would assume everyone would?)

  #load each library
    require(dplyr)
    require(tidyr)
    require(stringr)
    require(plyr)
  
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
  NameExtract=str_extract(NamePops, "[A-Z]{3}" ) # extract the text from the individuals names to denote population
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
  
    temp2$Pop
    
    ### FOR SOME REASON THIS LOOP IS NOT ADDING THE FINAL POP!!!
    
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
    
    how.many.bang <- Nmothers+Nfathers
    how.many.pops <- length(PopLengths)
    family.list <- NULL
  
    for(j in 1:how.many.pops){  
      family.hold <- NULL
    
   
      from <- min(min.max.out[j,])
      to <- max(min.max.out[j,])
      
      if(how.many.bang<=length(c(to:from))){
      parent.hold <- sample(to:from, how.many.bang)
      }
      if(how.many.bang>length(c(to:from))){
      parent.hold <- sample(to:from, how.many.bang, replace = TRUE)
      }
  
   
  
    
    ## For simplicty - the mothers are the first Nmother individuals identified at random
    ### get the genotypes of the mothers and fathers from temp2
    mother.get <- temp2[parent.hold[1:Nmothers],]
    father.get <- temp2[parent.hold[((Nmothers+1):length(parent.hold))],]
    
    ### when collected from temp2 - has "Pop" on end - remove
    mother.get <- mother.get[-which(names(mother.get) == "Pop")]
    father.get <- father.get[-which(names(father.get) == "Pop")]
    
    family.hold <- rbind(mother.get, father.get)
    
    assign(x = paste("Family.out", j, sep = "."), value = family.hold, envir = globalenv())
    family.list <- rbind(family.list, paste("Family.out", j, sep = "."))
      
    }
## MAKE PURE CROSSES
uPop.names = unique(NameExtract)
get.cross.names = NULL
allele.hold <- NULL
pure.prog <- txtProgressBar(min=0, max = length(family.list), style = 3)
print("Pure Cross Progress")
for(h in 1:length(family.list)){

family.make = get(family.list[h])
#print(paste0("Pure cross progress = ", ((h/length(family.list))*100)))  
setTxtProgressBar(pure.prog, h)
for(i in 1:max(Nfamilies$Family)){ ## get families within population
  
  get.mom.who.bangs <- Nfamilies$Moth[i]
  get.dad.who.bangs <- Nfamilies$Fath[i]
  
  mom.who.bangs <- family.make[get.mom.who.bangs,]
  dad.who.bangs <- family.make[(max(Nfamilies$Moth)+get.dad.who.bangs),]
  
  #mom.who.bangs <- cbind(NA, mom.who.bangs)
  #dad.who.bangs <- cbind(NA, dad.who.bangs)
  
  off.name.hold = paste(uPop.names[h], i, sep = ".") ## this will hold the name of the mother and father so that can check offspring assignment against
  ## Can replace 50 with number of offsrping later
  off.gens.hold = data.frame(matrix(vector(), 50, (length(ColumnData.Dup)+1))) ### have to create a data frame with the same number of columns as there are alleles, + 1 for the ID, and number of 
                                                            ### rows as there are simulated offpsring 
                                                              ### because it this data frame is specified here, it will thus blank it each time the loop goes through parent pairs
  
    mom.matrix <- data.frame(matrix(vector(), length(mom.who.bangs)/2, 2))
mom.matrix$X1 = t(mom.who.bangs[c(T,F)])
mom.matrix$X2 = t(mom.who.bangs[c(F,T)])


dad.matrix <- data.frame(matrix(vector(), length(dad.who.bangs)/2, 2))
dad.matrix$X1 <- t(dad.who.bangs[c(T,F)])
dad.matrix$X2 <- t(dad.who.bangs[c(T,F)])
    for(j in 1:50){ ### for each parent pair, want to make a number of simulated offpsring. 
    
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
      #print(paste0("J progress = ", ((j/50)*100)))  
      } # J Closure
    gen.dat.file = paste("Pure", off.name.hold,  sep = "_") ## makes an object name consisting of genotypes + the two parent names
    assign(x = gen.dat.file, value = off.gens.hold, envir = globalenv()) ## assigns the object name to the generated simulated offspring genotypes, and then assigns it to the global env.
  get.cross.names = rbind(get.cross.names, gen.dat.file)
} # I Closure
} ## H Closure



pop.fam.reps = length(get.cross.names)/how.many.pops
#length(get.cross.names)

#l = 1
#test = get(get.cross.names[1:5])
get.pure.names = NULL
to.get <- NULL
get.string <- NULL
objects <- NULL
geno.list <- NULL
geno.unlist <- NULL

for(l in 1:length(uPop.names)){
  get.string = paste0("Pure_", uPop.names[l])
  to.get = which(str_detect(ls(), pattern = get.string) == TRUE)
  objects = ls()
  geno.list = list(mget(objects[c(min(to.get):max(to.get))]))
   geno.unlist = unlist(geno.list, recursive = FALSE)
   geno.rbind = do.call(rbind, geno.unlist)
  
   
   pure.geno.out <- sample_n(tbl = geno.rbind, size = 50, replace = FALSE)

  assign(x = get.string, value = pure.geno.out, envir = globalenv())
  get.pure.names = rbind(get.pure.names, get.string)
} 



### Now make the F1

get.pure.names.df <- data.frame(get.pure.names)
F1.fam.name.list <- NULL
for(r in 1:nrow(get.pure.names.df)){
  
  to.get <- as.character(get.pure.names.df[r,])
  F1.offspring <- get(to.get)
  
  F1.parent.hold <- sample_n(tbl = F1.offspring, size = how.many.bang, replace = FALSE)
  
  F1.fam.name <- paste(to.get, "parents", sep = "_")
  assign(x = F1.fam.name, value = F1.parent.hold, envir = globalenv())
  F1.fam.name.list <- c(F1.fam.name.list, F1.fam.name)
  
}

F1.par.combs <- combn(F1.fam.name.list, 2, paste, collapse = '+')
F1.par.combs <- data.frame(do.call(rbind, str_split(F1.par.combs, "\\+"))) 


 F1.family.list <- NULL
  
    for(j in 1:nrow(F1.par.combs)){  
      family.hold <- NULL
    
    mother.get.word.F1 <- as.character(F1.par.combs[j,1])
    father.get.word.F1 <- as.character(F1.par.combs[j,2])
    
   # mother.pop.name.F1 <- str_extract(string = mother.get.word.F1, pattern = "_..._")
   # father.pop.name.F1 <- str_extract(string = father.get.word.F1, pattern = "_..._")
    mother.pop.name.F1 <- str_extract(string = mother.get.word.F1, pattern = "[A-Z]{3}")
    father.pop.name.F1 <- str_extract(string = father.get.word.F1, pattern = "[A-Z]{3}")
    
    mother.pop.F1 <- get(mother.get.word.F1) 
    father.pop.F1 <- get(father.get.word.F1)
    
    mother.get.F1 <- sample_n(tbl = mother.pop.F1, size = Nmothers)
    father.get.F1 <- sample_n(tbl = father.pop.F1, size = Nfathers)
    
    
    family.hold.F1 <- rbind(mother.get.F1, father.get.F1)
    
    Fam.name.F1 <- paste("F1.par", mother.pop.name.F1, father.pop.name.F1,  sep = ".")
    
    assign(x = Fam.name.F1, value = family.hold.F1, envir = globalenv())
    F1.family.list <- rbind(F1.family.list, Fam.name.F1)
      
    }
 
 ## MAKE THE F1
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
  
  mom.who.bangs <- F1.family.make[get.mom.who.bangs,-1]
  dad.who.bangs <- F1.family.make[(max(Nfamilies$Moth)+get.dad.who.bangs),-1]
  
  #mom.who.bangs <- cbind(NA, mom.who.bangs)
  #dad.who.bangs <- cbind(NA, dad.who.bangs)
  
  off.name.hold = paste(F1.family.list[h], i, sep = ".X.") ## this will hold the name of the mother and father so that can check offspring assignment against
  ## Can replace 50 with number of offsrping later
  off.gens.hold = data.frame(matrix(vector(), 50, (length(ColumnData.Dup)+1))) ### have to create a data frame with the same number of columns as there are alleles, + 1 for the ID, and number of 
                                                            ### rows as there are simulated offpsring 
                                                              ### because it this data frame is specified here, it will thus blank it each time the loop goes through parent pairs
    
        mom.matrix <- data.frame(matrix(vector(), length(mom.who.bangs)/2, 2))
mom.matrix$X1 = t(mom.who.bangs[c(T,F)])
mom.matrix$X2 = t(mom.who.bangs[c(F,T)])


dad.matrix <- data.frame(matrix(vector(), length(dad.who.bangs)/2, 2))
dad.matrix$X1 <- t(dad.who.bangs[c(T,F)])
dad.matrix$X2 <- t(dad.who.bangs[c(T,F)])
    for(j in 1:50){ ### for each parent pair, want to make a number of simulated offpsring. 
    
dad.out <- apply(t(dad.matrix), 2, sample, 1)
mom.out <- apply(t(mom.matrix), 2, sample, 1)
      
      off.name = paste(off.name.hold, j, sep = "_") ## generates an individual offspring name composed of its parents names (off.name.hold), and the number of the offspring
      #allele.hold = data.frame(allele.hold) ### the allele.hold object needs to be coerced into a data frame format to work with it optimally
      colnames(allele.hold) = NULL ## remove the column names - makes it simpler to work with
    
    allele.hold<- data.frame(c(rbind(mom.out, dad.out)))
    off.gens.hold[j,2:length(off.gens.hold)] = t(allele.hold) 
      off.gens.hold[j,1] = off.name ## this now adds the offspring name to the front row
      colnames(off.gens.hold) = c("ID", ColumnData.Dup)
      allele.hold = NULL ## nulls out the allele.hold object so it can be recycled back through the loop
      #print(paste0("J progress = ", ((j/50)*100)))  
      } # J Closure
  
  gen.dat.file <- unlist(str_split(string = off.name.hold, pattern = "F1.par."))[2]
    
  gen.dat.file = paste("F1", gen.dat.file,  sep = "_") ## makes an object name consisting of genotypes + the two parent names
    assign(x = gen.dat.file, value = off.gens.hold, envir = globalenv()) ## assigns the object name to the generated simulated offspring genotypes, and then assigns it to the global env.
  F1.get.cross.names = rbind(F1.get.cross.names, gen.dat.file)
} # I Closure

#print(paste0("I progress = ", ((i/max(Nfamilies$Family))*100)))  
} # H Closure
 
 #l=1
 get.F1.names = NULL
for(l in 1:length(F1.family.list)){
  get.string.F1 = (F1.family.list[l])
  to.get.F1 = which(str_detect(ls(), pattern = get.string) == TRUE)
  objects.F1 = ls()
  geno.list.F1 = list(mget(objects.F1[c(min(to.get.F1):max(to.get.F1))]))
   geno.unlist.F1 = unlist(geno.list.F1, recursive = FALSE)
   geno.rbind.F1 = do.call(rbind, geno.unlist.F1)
  
   
   geno.out.F1 <- sample_n(tbl = geno.rbind.F1, size = 50, replace = FALSE)

  assign(x = get.string.F1, value = geno.out.F1, envir = globalenv())
  get.F1.names = rbind(get.F1.names, get.string.F1)
} 
 
 #rm(F1.get.cross.names)
F2.par.combs <- combn(get.F1.names, 2, paste, collapse = '+')
F2.par.combs <- data.frame(do.call(rbind, str_split(F2.par.combs, "\\+"))) 

F2.family.list <- NULL

for(q in 1:nrow(F2.par.combs)){  
      family.hold <- NULL
    
    mother.get.word.F2 <- as.character(F2.par.combs[q,1])
    father.get.word.F2 <- as.character(F2.par.combs[q,2])
    
    mother.pop.name.F2 <- unlist(str_split(string = mother.get.word.F2, pattern = "F1.par."))[2]
    father.pop.name.F2 <-  unlist(str_split(string = father.get.word.F2, pattern = "F1.par."))[2]
    
    mother.pop.F2 <- get(mother.get.word.F2) 
    father.pop.F2 <- get(father.get.word.F2)
    
    mother.get.F2 <- sample_n(tbl = mother.pop.F2, size = Nmothers)
    father.get.F2 <- sample_n(tbl = father.pop.F2, size = Nfathers)
    
    
    family.hold.F2 <- rbind(mother.get.F2, father.get.F2)
    
    Fam.name.F2 <- paste("F2.par.", mother.pop.name.F2, "x", father.pop.name.F2,  sep = "")
    
    assign(x = Fam.name.F2, value = family.hold.F2, envir = globalenv())
    F2.family.list <- rbind(F2.family.list, Fam.name.F2)
      
    }


## MAKE THE F2
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
  
  off.name.hold = paste(F2.family.list[h], i, sep = ".X.") ## this will hold the name of the mother and father so that can check offspring assignment against
  ## Can replace 50 with number of offsrping later
  off.gens.hold = data.frame(matrix(vector(), 50, (length(ColumnData.Dup)+1))) ### have to create a data frame with the same number of columns as there are alleles, + 1 for the ID, and number of 
                                                            ### rows as there are simulated offpsring 
                                                              ### because it this data frame is specified here, it will thus blank it each time the loop goes through parent pair
    
        mom.matrix <- data.frame(matrix(vector(), length(mom.who.bangs)/2, 2))
mom.matrix$X1 = t(mom.who.bangs[c(T,F)])
mom.matrix$X2 = t(mom.who.bangs[c(F,T)])


dad.matrix <- data.frame(matrix(vector(), length(dad.who.bangs)/2, 2))
dad.matrix$X1 <- t(dad.who.bangs[c(T,F)])
dad.matrix$X2 <- t(dad.who.bangs[c(T,F)])
    for(j in 1:50){ ### for each parent pair, want to make a number of simulated offpsring. 
    
dad.out <- apply(t(dad.matrix), 2, sample, 1)
mom.out <- apply(t(mom.matrix), 2, sample, 1)
      
      off.name = paste(off.name.hold, j, sep = "_") ## generates an individual offspring name composed of its parents names (off.name.hold), and the number of the offspring
      #allele.hold = data.frame(allele.hold) ### the allele.hold object needs to be coerced into a data frame format to work with it optimally
      colnames(allele.hold) = NULL ## remove the column names - makes it simpler to work with
    
    allele.hold<- data.frame(c(rbind(mom.out, dad.out)))
    off.gens.hold[j,2:length(off.gens.hold)] = t(allele.hold) 
      
      
      off.gens.hold[j,1] = off.name ## this now adds the offspring name to the front row
      colnames(off.gens.hold) = c("ID", ColumnData.Dup)
      allele.hold = NULL ## nulls out the allele.hold object so it can be recycled back through the loop
      #print(paste0("J progress = ", ((j/50)*100)))  
      } # J Closure
  
  gen.dat.file <- unlist(str_split(string = off.name.hold, pattern = "F2.par."))[2]
    
  gen.dat.file = paste("F2", gen.dat.file,  sep = "_") ## makes an object name consisting of genotypes + the two parent names
    assign(x = gen.dat.file, value = off.gens.hold, envir = globalenv()) ## assigns the object name to the generated simulated offspring genotypes, and then assigns it to the global env.
  F2.get.cross.names = rbind(F2.get.cross.names, gen.dat.file)
} # I Closure

#print(paste0("I progress = ", ((i/max(Nfamilies$Family))*100)))  
} # H Closure
 
# l=1
 get.F2.names = NULL
for(l in 1:length(F2.family.list)){
  get.string.F2 = (F2.family.list[l])
  to.get.F2 = which(str_detect(ls(), pattern = get.string) == TRUE)
  objects.F2 = ls()
  geno.list.F2 = list(mget(objects.F2[c(min(to.get.F2):max(to.get.F2))]))
   geno.unlist.F2 = unlist(geno.list.F2, recursive = FALSE)
   geno.rbind.F2 = do.call(rbind, geno.unlist.F2)
  
   
   geno.out.F2 <- sample_n(tbl = geno.rbind.F2, size = 50, replace = FALSE)

  assign(x = get.string.F2, value = geno.out.F2, envir = globalenv())
  get.F2.names = rbind(get.F2.names, get.string.F2)
} 
  

 BC.names.list <- rbind(get.pure.names, get.F1.names)
 
BC.par.combs <- combn(BC.names.list, 2, paste, collapse = '+')
BC.par.combs <- data.frame(do.call(rbind, str_split(BC.par.combs, "\\+"))) 


BC.par.combs <- BC.par.combs[-which(str_detect(string = BC.par.combs$X1, pattern = "F1.par.")==TRUE),]
BC.par.combs <- BC.par.combs[-which(str_locate(string = BC.par.combs$X2, pattern = "Pure_")==TRUE),]


BC.family.list <- NULL

for(q in 1:nrow(BC.par.combs)){  
      family.hold <- NULL
    
    mother.get.word.BC <- as.character(BC.par.combs[q,1])
    father.get.word.BC <- as.character(BC.par.combs[q,2])
    
    mother.pop.name.BC <- unlist(str_split(string = mother.get.word.BC, pattern = "Pure_"))[2]
    father.pop.name.BC <-  unlist(str_split(string = father.get.word.BC, pattern = "F1.par."))[2]
    
    mother.pop.BC <- get(mother.get.word.BC) 
    father.pop.BC <- get(father.get.word.BC)
    
    mother.get.BC <- sample_n(tbl = mother.pop.BC, size = Nmothers)
    father.get.BC <- sample_n(tbl = father.pop.BC, size = Nfathers)
    
    
    family.hold.BC <- rbind(mother.get.BC, father.get.BC)
    
    Fam.name.BC <- paste("BC.pure", mother.pop.name.BC, sep = ".")
    Fam.name.BC <- paste(Fam.name.BC, father.pop.name.BC,  sep = "x")
    
    
    assign(x = Fam.name.BC, value = family.hold.BC, envir = globalenv())
    BC.family.list <- rbind(BC.family.list, Fam.name.BC)
      
    } 
 

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
  
  off.name.hold = paste(BC.family.list[h], i, sep = ".X.") ## this will hold the name of the mother and father so that can check offspring assignment against
  ## Can replace 50 with number of offsrping later
  off.gens.hold = data.frame(matrix(vector(), 50, (length(ColumnData.Dup)+1))) ### have to create a data frame with the same number of columns as there are alleles, + 1 for the ID, and number of 
                                                            ### rows as there are simulated offpsring 
                                                              ### because it this data frame is specified here, it will thus blank it each time the loop goes through parent pairs
   
    
        mom.matrix <- data.frame(matrix(vector(), length(mom.who.bangs)/2, 2))
mom.matrix$X1 = t(mom.who.bangs[c(T,F)])
mom.matrix$X2 = t(mom.who.bangs[c(F,T)])


dad.matrix <- data.frame(matrix(vector(), length(dad.who.bangs)/2, 2))
dad.matrix$X1 <- t(dad.who.bangs[c(T,F)])
dad.matrix$X2 <- t(dad.who.bangs[c(T,F)])
    for(j in 1:50){ ### for each parent pair, want to make a number of simulated offpsring. 
    
dad.out <- apply(t(dad.matrix), 2, sample, 1)
mom.out <- apply(t(mom.matrix), 2, sample, 1)
      
      off.name = paste(off.name.hold, j, sep = "_") ## generates an individual offspring name composed of its parents names (off.name.hold), and the number of the offspring
      #allele.hold = data.frame(allele.hold) ### the allele.hold object needs to be coerced into a data frame format to work with it optimally
      colnames(allele.hold) = NULL ## remove the column names - makes it simpler to work with
    
    allele.hold<- data.frame(c(rbind(mom.out, dad.out)))
    off.gens.hold[j,2:length(off.gens.hold)] = t(allele.hold) 
      
      
      off.gens.hold[j,1] = off.name ## this now adds the offspring name to the front row
      colnames(off.gens.hold) = c("ID", ColumnData.Dup)
      allele.hold = NULL ## nulls out the allele.hold object so it can be recycled back through the loop
      #print(paste0("J progress = ", ((j/50)*100)))  
      } # J Closure
  
  gen.dat.file <- unlist(str_split(string = off.name.hold, pattern = "BC.par."))[2]
    
  gen.dat.file = paste("BC", gen.dat.file,  sep = "_") ## makes an object name consisting of genotypes + the two parent names
    assign(x = gen.dat.file, value = off.gens.hold, envir = globalenv()) ## assigns the object name to the generated simulated offspring genotypes, and then assigns it to the global env.
  BC.get.cross.names = rbind(BC.get.cross.names, gen.dat.file)
} # I Closure

#print(paste0("I progress = ", ((i/max(Nfamilies$Family))*100)))  
} # H Closure
 


get.BC.names = NULL
for(l in 1:length(BC.family.list)){
  get.string.BC = (BC.family.list[l])
  to.get.BC = which(str_detect(ls(), pattern = get.string) == TRUE)
  objects.BC = ls()
  geno.list.BC = list(mget(objects.BC[c(min(to.get.BC):max(to.get.BC))]))
   geno.unlist.BC = unlist(geno.list.BC, recursive = FALSE)
   geno.rbind.BC = do.call(rbind, geno.unlist.BC)
  
   
   geno.out.BC <- sample_n(tbl = geno.rbind.BC, size = 50, replace = FALSE)

  assign(x = get.string.BC, value = geno.out.BC, envir = globalenv())
  get.BC.names = rbind(get.BC.names, get.string.BC)
} 
  

get.pure.names

get.F1.names

get.F2.names

get.BC.names

nrow(BC.parXBDNXBDN.DHB)

test.names <- get.pure.names

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
n.pop.F1 <- nrow(get.F1.names)
n.pop.F2 <- nrow(get.F2.names)
n.pop.BC <- nrow(get.BC.names)

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
       for(v in 1:n.pop.F1){
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
