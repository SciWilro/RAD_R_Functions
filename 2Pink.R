



### NOTE This only deals with inclusion thresholds -

## We can give them the ability to specify an output if we want - do we want to?

plink.r <- function(ped, map = NULL, hwe = NULL, maf = NULL, geno = NULL, 
  mind = NULL){
 
  require("stringr")
  
  ### Function to run Plink through R
  
  # For plink, at minimum should have ped file, map file, HWE value, maf value, geno 
  ## value, mind value and specify file output name
  ## Default values are from Plink
  # maf - minor allele frequency (0.01)
  # geno - maximum per-snp missing (0.1) - new version has no default - what should we do?
  # mind - maximum per-individual missing (0.1) - new versin has no default - what should we do?
  # hwe - Hardy-Weinberg disequilibrium p-value (exact) (0.001)
  
  
  ## Decided to remove the prompting for values because I realized that PLINK has defaults
  
  ped = as.character(ped) ## make the ped file name a character string if it is not
  
  ## if they have entered a map file name, make it a character string as well
  if(length(map) != 0){
    map = as.character(map)
  }
  
  
  ## check to see if a value has been entered for maf - if one has not, default to 0.01
  if(length(maf) == 0){
    maf = 0.01
    }
      
   if(is.numeric(maf) == FALSE){ # a non-numerical value has been entered for maf, stop the function and throw an error
     stop("Invalid maf value")
  }
  
    if(maf != 0 | maf != 1){ # make sure maf is between 0 and 1
      maf.checks.out = "Yes"
      }
        maf = as.numeric(maf)
  
   ## The code for maf is repeated for hwe, geno and mind
        
  if(length(hwe) == 0){
    hwe = 0.0001
    }
  
    if(is.numeric(hwe) == "FALSE"){
      stop("Invalid HWE value")
      }
  
      if(hwe >= 0 & hwe <= 1){
        hwe.checks.out = "Yes"
        }
      hwe = as.numeric(hwe)
  
  
   
    if(length(geno) == 0){
      geno = 0.01
    }
  
      if(is.numeric(geno) == "FALSE"){
        stop("Invalid geno value")
        }
      
      if(geno >= 0 & geno <= 1){
        geno.checks.out = "Yes"
        }
    geno = as.numeric(geno)
  
    
  if(length(mind) == 0){
    mind = 0.1
   
  }

    if(is.numeric(mind) == "FALSE"){
      stop("Invalid minor allele frequency")
        }
    mind = 0.01
    
      if(mind >= 0 & mind <= 1){
        mind.checks.out = "Yes"
    }

  mind = as.numeric(mind)
  
  
  
   ### check to see if they have entered a map file name
  if(length(map) == 0){
    map.yn = readline("Do the ped and map files have the same file name? y/n ") ## if they have not entered one, check that the name
      ## of the map file is the same as the ped file
   
    while(map.yn != "y" || map != "n"){ ## RYAN -> the while statement doesn't work yet, the ifs below do
       map.yn = readline("Do the ped and map files have the same file name? y/n ") ## queries if names are same, accepts input
        }
          if(map.yn == "n"){ ## if file names are not the same, queries user for filename and then accepts input
            map = readline("Please enter the map file name ")
            }
              if(map.yn == "y"){ ## if they hvae the same name, set the name of hte map file equal to the ped file
                map = ped
                }
  }  
  
  ## check to see how the unser has entered the ped and map file names - this seems perhaps redundant, but ...
  
  ## look for file names using sequences
  ped.seq = ".ped"
  map.seq = ".map"
  
  ## make NULL objects to hold the extracted data
  ped.str.extract = NULL
  map.str.extract = NULL
  
## string extract the .ped ending
  ped.str.extract = str_extract(string = ped, pattern = ped.seq)
  ## if there was no .ped, then there will be no ending, so add it
  if(length(ped.str.extract) != 3){
    ped = paste0(ped, ".ped")
  }
   ## repeat for map
  map.str.extract = str_extract(string = map, pattern = map.seq)
  
  if(length(map.str.extract != 3)){
    map = paste0(map, ".map")
  }

  ### since now know that the file names have been appended with .ped and .map, will remove them so it matches the PLINK command line
  ped.out = str_split(string = ped, pattern = ".ped")[[1]][1]
  map.out = str_split(string = map, pattern = ".map")[[1]][1]
  
  
  ### RYAN - here is where I need your help - I don't know how to tell R to look for different things in different places
 # move.file =
    
    
    
   # if(ped.checks.out == "Yes" & mind.checks.out == "Yes" & geno.checks.out == "Yes" & hwe.checks.out == "Yes" & maf.checks.out == "Yes"){
  #    should.I.run.it == "Yes"
   # }
  #  if(ped.checks.out != "Yes" | mind.checks.out != "Yes" | geno.checks.out != "Yes" & hwe.checks.out != "Yes" | maf.checks.out != "Yes"){
  #    should.I.run.it == "No"
   # }
  
    ## format the variables as they are to be entered in the command line
    mind.run <- paste0(" --mind ", mind)
    
    maf.run <- paste0(" --maf ", maf)
    
    geno.run <- paste0(" --geno ", geno)
  
    hwe.run = paste0(" --hwe ", hwe)
    
    ped.run = paste0(" --file ", ped.out)
    
    map.run = paste0(" --map ", map.out)
    
    
    ## create a file name for the exported data which indicates the levels of the variables used
    maf.out <- paste0("_maf_", maf)
    geno.out <- paste0("_geno_", geno)
    mind.out <- paste0("_mind_", mind)
    hwe.out = paste0("_hwe_", hwe)
    ## creates file name
    output = paste0(ped.out, maf.out, geno.out, mind.out, hwe.out)
      
    ## this is the code to be run on the command line
  to.run.code = paste0("./plink --noweb", ped.run, map.run, mind.run, geno.run, maf.run, hwe.run, 
    " --out ", output, "--recode")
  
  ## for now just prints until I can figure out hwo to get R to go to the command line.
  print(to.run.code)
 # copy.new.file = 
  
  #hold.working.dir = getwd()
  
  }

plink.r(ped = "test", map = "map.dummy")

## Works
