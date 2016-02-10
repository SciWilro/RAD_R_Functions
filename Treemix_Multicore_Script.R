where.treemix <- "/home/ian/treemix-1.12/src"
infile <- "NeutralTreemixInput.frq.gz"
maxtrees <- 10
treeoutput<-"NeutralOutput"

Treemix.Auto <- function(where.treemix, infile, maxtrees, treeoutput){
  require(parallel)
  h.cores <- detectCores()
  
  calltreemix<- paste0("cd ", where.treemix)
  
  input.text<- paste0("./treemix ", "-i ", infile)
  
  jobs.vector <- NULL
  for(b in 0:maxtrees){
    
    tree.do <- paste("-m", b)
    
    treeoutputhold <- paste0("-o ", treeoutput, "_", b)
    

    
      path.hold <- paste0(calltreemix, "; ", input.text, " ", tree.do, " ", treeoutputhold, " -se")
      jobs.vector <- c(jobs.vector, path.hold)
      
    }
    
  mclapply(X = jobs.vector, FUN = system, mc.cores = h.cores)
  