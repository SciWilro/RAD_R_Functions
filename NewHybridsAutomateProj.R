find.the.fucking.hybrids <- function(folder.data, where.NH, burnin, sweeps){
    
    files.anal <- list.files(path = folder.data)
    
    dir.create(path = paste0(folder.data, "NH.Temp"))
    where.temp <- paste0(folder.data, "NH.Temp", "/")
    
    h.cores <- detectCores()
    
    
    
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
    
    
    ## now to execute the command
    
    options(scipen = 999) ### if the number of digits is too big in either burnin or sweeps
    ## R will output in scientific notation, which is not interpreted correctly by NH
    
    burnin.do <- paste("--burn-in", burnin, sep=" ")
    sweeps.do <- paste("--num-sweeps", sweeps, sep=" ")
    ## replace with files.anal
    where.temp2 <- gsub(x = where.temp, pattern = " ", replacement = "\\ ", fixed = T)
    jobs.vector <- NULL
    
    for(b in 1:length(NH.copy.list)){
        
        b.copy <- NH.copy.list[b]
        file.do <- paste("-d", files.anal[b])
        what.temp <- paste0(where.temp2, b.copy)
        path.hold <- paste("cd", paste0(what.temp, ";"), "newhybsng", file.do, burnin.do, sweeps.do, "--no-gui", sep = " ")
        jobs.vector <- c(jobs.vector, path.hold)
        
    }
    
    
    if(length(jobs.vector)<h.cores){
        mclapply(X = jobs.vector, FUN = system, mc.cores = length(jobs.vector))
    }
    if(length(jobs.vector)>h.cores){
        mclapply(X = jobs.vector, FUN = system, mc.cores = h.cores)
    }
    
}