#GenePopData <- read.table("Salmon_South_Coast.txt", header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
GenePopData <- read.table("Test_GenePop_Data.txt", header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
GenePop <- GenePopData
ColumnData.Dup = NULL
temp = NULL
temp2 = NULL
a = NULL
b = NULL
Col.name.orig = NULL


do.sex <- function(genotype, Nmothers, Nfathers, Noff, pop, F1, F2 = No, BC = Yes){
  
  GenePop <- genotype
  
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
    ColumnData <- GenePop[1:(Pops[1]-1),"data"]  ### SNP Names
    snpData <- GenePop[Pops[1]:NROW(GenePop),]  ### Genotypes

#Get a datafile with just the snp data no pops
    tempPops <- which(snpData$data=="Pop")
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
  
    temp2$Pop
    
    PopLengths <- table(temp2$Pop)[-length(table(temp2$Pop))]
    
    ## Need to be able to tell what row each individual is in, and what population it is
    ind.vector = c(1:nrow(temp)) ### make a vector that is the number of individuals
    ind.matrix = data.frame(temp2$Pop, ind.vector) ## add populatuions to that
    
    Nmothers = 5
    Nfathers = 4
    
    ### make a dataframe to hold the number of mother and father - DF will be compoised of NULLs
    ## need to make the dataframe as long as the longer of the number of mothrs or fathers
    if(Nmothers > Nfathers){
      Nfamilies = matrix(data = rep(NA), nrow = Nmothers, ncol = 2, dimnames = NULL)
      colnames(Nfamilies) <- c("Moth", "Fath")
      Nfamilies = as.data.frame(Nfamilies)
    }
    if(Nmothers < Nfathers){
      Nfamilies = matrix(data = rep(NA), nrow = Nfathers, ncol = 2)
    }
    if(Nmothers == Nfathers){
      Nfamilies = matrix(data = rep(NA), nrow = Nmothers, ncol = 2)
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
    for(i in 1:how.many.bang){
   
      from <- min(min.max.out[i,])
      to <- max(min.max.out[i,])
      parent.hold <- sample(to:from, how.many.bang)
   
        }
    
    ## For simplicty - the mothers are the first Nmother individuals identified at random
    ### get the genotypes of the mothers and fathers from temp2
    mother.get <- temp2[parent.hold[1:Nmothers],]
    father.get <- temp2[parent.hold[((Nmothers+1):length(parent.hold))],]
    
    ### when collected from temp2 - has "Pop" on end - remove
    mother.get <- mother.get[-which(names(mother.get) == "Pop")]
    father.get <- father.get[-which(names(father.get) == "Pop")]
    
    family.hold <- rbind(mother.get, father.get)
    
    Family.out = family.hold
    
    assign(x = paste("Family.out", j, sep = "."), value = family.hold, envir = globalenv())
    family.list <- rbind(family.list, paste("Family.out", j, sep = "."))
      
    }
    

}

h = 1
i = 1
k = 1
uPop.names = unique(NameExtract)
get.cross.names = NULL

for(h in 1:length(family.list)){

family.make = get(family.list[h])
print(paste0("H progress = ", ((h/length(family.list))*100)))  

for(i in 1:max(Nfamilies$Family)){ ## get families within population
  
  get.mom.who.bangs <- Nfamilies$Moth[i]
  get.dad.who.bangs <- Nfamilies$Fath[i]
  
  mom.who.bangs <- Family.out[get.mom.who.bangs,]
  dad.who.bangs <- Family.out[(max(Nfamilies$Moth)+get.dad.who.bangs),]
  
  #mom.who.bangs <- cbind(NA, mom.who.bangs)
  #dad.who.bangs <- cbind(NA, dad.who.bangs)
  
  off.name.hold = paste(uPop.names[h], uPop.names[h], i, sep = ".X.") ## this will hold the name of the mother and father so that can check offspring assignment against
  
  off.gens.hold = data.frame(matrix(vector(), 50, (length(ColumnData.Dup)+1))) ### have to create a data frame with the same number of columns as there are alleles, + 1 for the ID, and number of 
                                                            ### rows as there are simulated offpsring 
                                                              ### because it this data frame is specified here, it will thus blank it each time the loop goes through parent pairs
    for(j in 1:50){ ### for each parent pair, want to make a number of simulated offpsring. 
    
    segregate.mom = sample.int(n = 2, size = length(mom.who.bangs), replace = T) ### creates a vector of randomly generated 1 and 2, of length equal to the number of loci tested. 
                                                                ### This will allow to randomly choose which of the two potential alleles for each loci to include in offspring
    segregate.dad = sample.int(n = 2, size = length(ColumnData.Dup), replace = T)  ## repeat for father so his alleles are chosen independently as well
      position.alelle <- c((0:length(ColumnData)*2))
      position.alelle <- position.alelle[-length(position.alelle)]
   
      
                                                 #c(0,2,4,6,8,10,12,14) ## a vector of numbers, such that it, plus the number generated for the mother or father, will equal the row number in which the ## allele to be passed on to the simulated offpspring
      
      allele.hold = NULL ## creates an object to hold the offspring alleles - also nulls it each time it goes through the loop
         for(k in 1:length(position.alelle)){ ### this loop creates the simulated genotype for one offspring
              
              allele.mom = mom.who.bangs[(position.alelle[k] + segregate.mom[k])] ## holds the randomly chosen allele for the mother
              allele.dad = dad.who.bangs[(position.alelle[k] + segregate.dad[k])] ## holds the randomly chosen allele for the father
              colnames(allele.mom) = NULL ## remove the column names - this was causing an issue
              colnames(allele.dad) = NULL 
              seg.res = c(as.character(allele.mom[[1]]), as.character(allele.dad[[1]])) ### binds together the mothers and fathers allele that were passed on to the simulated offpsring for a given loci
              colnames(seg.res) = NULL ## again - this was causing problems
              allele.hold = c(allele.hold, seg.res) ## this collects the alleles for a simulated offspring for each loci, and holds them together
  #print(paste0("K progress = ", ((k/length(position.alelle))*100)))  
              
          } # K Closure
      
      off.name = paste(off.name.hold, j, sep = "_") ## generates an individual offspring name composed of its parents names (off.name.hold), and the number of the offspring
      allele.hold = data.frame(allele.hold) ### the allele.hold object needs to be coerced into a data frame format to work with it optimally
      colnames(allele.hold) = NULL ## remove the column names - makes it simpler to work with
    
    off.gens.hold[j,2:length(off.gens.hold)] = allele.hold ## there was an issue with trying to just add the offspring name to the front of the allele.hold object, then adding
                                                              ### this to the off.gens.hold data frame - for some reason, it would only add the level of the 'factor' R was
                                                                ## interpreting the name as - so this starts the genotypes at the second position where they should be, and keeps
                                                                  ## the first blank for now
      off.gens.hold[j,1] = off.name ## this now adds the offspring name to the front row
      colnames(off.gens.hold) = ColumnData.Dup
      allele.hold = NULL ## nulls out the allele.hold object so it can be recycled back through the loop
      #print(paste0("J progress = ", ((j/50)*100)))  
      } # J Closure
    gen.dat.file = paste("genotypes", off.name.hold,  sep = "_") ## makes an object name consisting of genotypes + the two parent names
    assign(x = gen.dat.file, value = off.gens.hold, envir = globalenv()) ## assigns the object name to the generated simulated offspring genotypes, and then assigns it to the global env.
  get.cross.names = rbind(get.cross.names, gen.dat.file)
} # I Closure

print(paste0("I progress = ", ((i/max(Nfamilies$Family))*100)))  
} # H Closure


for 
