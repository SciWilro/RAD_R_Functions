NH.auto.win -> function(folder.data, where.NH, burnin, sweeps){

options(scipen = 999) ## have to change the global R settings because NH doesn't understand scientific notation'

files.anal <- list.files(path = folder.data) ## get a list of the files to analyze

dir.create(path = paste0(folder.data, "NH.Temp"))

where.temp <- paste0(folder.data, "NH.Temp", "/") ## need to know where to put the files yo

NH.loc.vec <- paste(where.temp, "newhybrids", sep = "/") ## put some folders up in that temp dir

## get some copying going
for(i in 1:length(files.anal)){
file.copy(from = where.NH, to = where.temp, recursive = TRUE) ## copies whole folder -- not as necessery in Windows version
NH.rename <- paste(NH.loc.vec, i, sep="_") ## get vector of where the folder is, then add a number to it
file.rename(from = NH.loc.vec, to = NH.rename) ## rename
}

## Copy the files to be analzyed into the NH folder copies - one file per copy
NH.copies <- list.files(where.temp) ## list of copies of NH in the temp folder

## slide dem files.anal in
for(j in 1:length(files.anal)){ ## for each file to be analyzed
to.file <- NH.copies[j] ## what copy of NH do you copy the file to be analyzed to?
from.file <- files.anal[j] ## the file to be copied
to.dir <- paste0(where.temp, to.file, "/") ## where is the NH copy - this is a directory
from.dir <- paste0(folder.data, from.file) ### where is the data file to be copied
file.copy(from = from.dir, to = to.dir) ## copy it on ovah
}


NH.copy.list <- list.files(path = where.temp)

## I think you can see what these do
burnin.do <- paste("--burn-in", burnin, sep=" ")
sweeps.do <- paste("--num-sweeps", sweeps, sep=" ")

## you know I love them NULL vecssssss
jobs.vector <- NULL

### commmands: what files and how to anal them?

for(b in 1:length(NH.copy.list)){
b.copy <- NH.copy.list[b]
file.do <- paste("-d", files.anal[b])
what.temp <- paste0(where.temp, b.copy)
path.hold <- paste("cd", paste0(what.temp, " &"), "newhybrids.exe", file.do, burnin.do, sweeps.do, "--no-gui", sep = " ") ## must separate
## the get directory form teh command with & in Windows version
jobs.vector <- rbind(jobs.vector, path.hold) ### rbind - will use ROW apply later
}

## make slaves
mc.cores <- min(length(jobs.vector), detectCores())
makecl <- makeCluster( mc.cores, outfile="" )

### parallele row apply - apply each command in the vector to the shell =
parRapply(cl = makecl, x=jobs.vector, FUN=shell)

## gotta stop the slaves somehow
stopCluster(makecl) ### I think this works???

}