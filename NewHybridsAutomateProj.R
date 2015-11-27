### need to find a way for R to launch interative version of NewHybrids

## need to see how many cores the computer has

require(parallel)
require(stringr)
## create a directory to work in

dir.create

folder.data <- file.choose()

folder.data <- str_split(string = folder.data, pattern = "/")
folder.data.l <- length(folder.data[[1]])
folder.data.c <- c(folder.data[[1]])
folder.out <- paste(folder.data.c[2:(folder.data.l-1)],collapse = "/")



folder.data <- "/Users/brendanwringe/Desktop/DFO Aquaculture Interaction/NH Trial/"
where.NH <- "/Users/brendanwringe/newhybrids/"
files.anal <- list.files(path = folder.data)

dir.create(path = paste0(folder.data, "NH.Temp"))
where.temp <- paste0(folder.data, "NH.Temp", "/")


## 1 - need to know how many cores there are in the computer, so can figure out how many times to spawn NH

h.cores <- detectCores()

## 2 - need to know how many files are in the folder to be analyzed
#length
#list.files(path = )

## copy the 

#length(files.anal)

#if(files.anal > h.cores)(
  #have.repeat = "Yes"#
#)

#file.copy(from = where.NH, to = where.temp, recursive = TRUE)

#NH.rename <- paste0(NH.loc.vec, "_2")
#file.rename(from = NH.loc.vec, to = NH.rename)

## Copy the entire NH folder i times to the temp folder you made
  ## use a loop here because then can iteratively name the copies so that the file.copy doesn't
  ## lose it because you have things named the same - also - how else are you goign to ID things later?

 NH.loc.vec <- paste(where.temp, "newhybrids", sep = "/") ## 
for(i in 1:length(files.anal)){
  
  file.copy(from = where.NH, to = where.temp, recursive = TRUE) ## copies whole folder
  
  NH.rename <- paste(NH.loc.vec, i, sep="_") ## get vector of where the folder is, then add a number to it
  
  file.rename(from = NH.loc.vec, to = NH.rename) ## rename
  

}


## Copy the files to be analzyed into the NH folder copies - one file per copy
NH.copies <- list.files(where.temp) ## list of copies of NH in the temp folder

for(j in 1:length(files.anal)){ ## for each file to be analyzed

  
  to.file <- NH.copies[j] ## what copy of NH do you copy the file to be analyzed to?
  from.file <- files.anal[j] ## the file to be copied
  
  to.dir <- paste0(where.temp, to.file, "/") ## where is the NH copy - this is a directory
  from.dir <- paste0(folder.data, from.file) ### where is the data file to be copied
  
  file.copy(from = from.dir, to = to.dir) ## copy it on ovah
  
}

## once NH has run, need to extract (copy) the output

### get a list of all the copies of NH you have made etc.
NH.copy.list <- list.files(path = where.temp)
# NH.copy.list # for internal code checking

## create a new folder to put the results in - remember the inception thing - be smart
dir.create(path = paste0(folder.data, "NH.Results"))

res.path.make <- paste0(folder.data, "NH.Results") ## get the path to the new results folder

# k = 1 # for internal code checking

## get the results generated for each copy of NH/file to be analyzed
for(k in 1:length(NH.copy.list)){
  
  ## make a nice name for the folders of resutls we are goign to make - so so pretty
  res.name <- paste0(files.anal[k], "_Results") # file name + results - not as pretty as it could be 
  dir.create(path = paste0(res.path.make, "/", res.name)) ## make that folder to put the resutls in
  
  where.the.results.to <- paste0(res.path.make, "/", res.name, "/") ## where did you make that folder?
  
  NH.copy.to.get <- NH.copy.list[k] ## what copy of NH are we goign to look in?
  find.res.vec <- list.files(path = paste0(where.temp, NH.copy.to.get), pattern = "aa-") ## find all the files that
    ## begin with "aa-" inside that copy of NH - big shout out to Eric Anderson for using a regular string in the results name!
  
  find.res.vec.path <- paste0(paste0(where.temp, NH.copy.to.get, "/"), find.res.vec) ## find the path to all the files ID'd
  file.copy(from = find.res.vec.path, to = where.the.results.to) ## copy those files to the proper results folder
  
  ## keep making it pretty by renaming the results files so you know what data set they are associated with
  rename.vec <- gsub(x = find.res.vec, pattern = "aa-", replacement = paste0(files.anal[k], "_"))
  
  ## again - need a vectr that shows their computr postion
  rename.vec.path <- paste0(where.the.results.to, rename.vec) ## vector with new names
  old.name.vec.path <- paste0(where.the.results.to, find.res.vec) ## vector with old names
  
  file.rename(from = old.name.vec.path, to = rename.vec.path) ## check that vectorized functionality in file.rename - mad proops!
  
}

### so ... the whole moving

## this - this right here will launch NH if R is run from the command line/terminal - but the program wont' automaticall run
  ## because it demands user specified data - goes out of command line, into NH preventing further 'code' 
system("cd /Users/brendanwringe/Desktop/DFO\\ Aquaculture\\ Interaction/NH\\ Trial/NH.temp/newhybrids_1/; newhybs")


system(command =  "/Users/brendanwringe/Desktop/DFO\\ Aquaculture\\ Interaction/NH\\ Trial/NH.temp/newhybrids_1/newhybs")

system("ls")
system("cd /Users/brendanwringe/Desktop/DFO\\ Aquaculture\\ Interaction/NH\\ Trial/NH.temp/newhybrids_1/; newhybs; NewSalDat.top.48.NH.txt")

system("cd /Users/brendanwringe/Desktop/DFO\\ Aquaculture\\ Interaction/NH\\ Trial/NH.temp/newhybrids_1/newhybs; NewSalDat.top.48.NH.txt")
system("cd /Users/brendanwringe/Desktop/DFO\\ Aquaculture\\ Interaction/NH\\ Trial/NH.temp/newhybrids_1/; newhybs; NewSalDat.top.48.NH.txt; 0; 0; 1; 4; 0", ignore.stdout=T)

system("cd /Users/brendanwringe/Desktop/DFO\\ Aquaculture\\ Interaction/NH\\ Trial/NH.temp/newhybrids_1/; newhybs",input = "NewSalDat.top.48.NH.txt; 0; 0; 1; 4; 0", ignore.stdout=T)

system("cd /Users/brendanwringe/Desktop/DFO\\ Aquaculture\\ Interaction/NH\\ Trial/NH.temp/newhybrids_1/; newhybs", ignore.stdout = T)


system("cd /Users/brendanwringe/Desktop/DFO\\ Aquaculture\\ Interaction/NH\\ Trial/NH.temp/newhybrids_1/; newhybs", input = "NewSalDat.top.48.NH.txt")
system(print("NewSalDat.top.48.NH.txt"))

cd /Applications/R.app/Contents/MacOS

cd ./Desktop/DFO\ Aquaculture\ Interaction/NH\ Trial/NH.Temp/newhybrids_1


