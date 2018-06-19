DecayFunc <- function(N0=10000,timeseries=NULL,years,nsim=1000,mean.decay=NULL,sd.decay=NULL,path=NULL,slopemodel=FALSE){

  require(ggplot2)
  
  #N0 - Initial population size
  #timeseries - data.frame with columns year and n corresponding to **sequential years** and population estimates.
  #years - this is the range of years in the model. If timeseries is NULL this must be specified.
  #nsim - number of times the simulation is to be repeated.
  #mean.decay - the mean decay rate (postive value) corresponding the the estimated % decline per year.
  #sd.decay - standard deviation estimated for the mean.decay rate.
  #path - if this value isn't NULL (default) the data used to generate the plot will be saved to this path. Note this must be a full file path ending in a .csv.
  #slopemodel - if timeseries information is provided the mean.decay and sd.decay will be calculated by function.
      # TRUE - mean.decay rate (and sd) for the time series will be estimated from an simple expodential linear model.
      # FALSE - mean.decay rate (and sd) will be estimated from the distribution of sequential differences (year to year changes in pop size).

  #estimate mean and sd of the decay rate using population data. 
  if(!is.null(timeseries)&!slopemodel){
    
    declines <- NULL
        for(i in 2:nrow(timeseries)){
          
          d <- 1-(timeseries[i,"n"]/timeseries[i-1,"n"])
          declines <- c(declines,d)
          
        }
    
    #set up model parameters
    mean.decay <- mean(declines,na.rm=T)
    sd.decay <- sd(declines,na.rm=T)
    years <- timeseries$year
   
  }
  
  if(!is.null(timeseries)&slopemodel){
    
    mod <- lm(log(n)~year,data=timeseries) # simple expodential model
    
    out <- summary(mod)
    
    mean.decay <- abs(coef(mod)[2])
    sd.decay <- coef(out)[2,2] # based on standard error of the slope
    years <- timeseries$year
    
    writeLines(paste("Estimated R2 of slope model =",round(out$adj.r.squared,2)))
  }

  #warming messages
  if(is.null(mean.decay)){stop("Must specify a mean.decay parameter if no timeseries information is provided")}
  if(is.null(sd.decay)){stop("Must specify a sd.decay parameter if no timeseries information is provided")}
  if(is.null(years)){stop("Must specify the range of years (e.g., years=2000:2018) if no timeseries information is provided")}
  
  #return information on decay rate
  writeLines(paste0("Estimated decay rate: ",abs(round(mean.decay,3))*100,"% per year"))
  writeLines(paste0("Estimated standard devation of decay rate: ",round(sd.decay,3)*100,"% per year"))

  #Recode the model parameters to fit existing code found online: https://sites.google.com/site/wild8390/spring-2017/lecture-lab-schedule/week-1-introduction/population-models-in-r
    r.mean<-mean.decay
    r.sd<-sd.decay
    T <- length(years)-1
  
    plotdata <- NULL # this will grow each loop -- this is a slow way of doing it but is probably fast enough and is at least easy to follow
    
    for (l in 1:nsim){
    
        r<-rnorm(T,r.mean,r.sd)
        
        #now use this in place of the constant r
        t<-N<-array(dim=T+1)
        
        #first element is initial value
        N[1]<-N0
        t[1]<-0
        for (i in 1:(T))
        {
          N[i+1]<-N[i]*exp(-r[i])
          t[i+1]=t[i]+1
        }
        #make it pretty
        temp <- data.frame(year=unique(years),N=N,nsim=l)
        
        plotdata <- rbind(plotdata,temp)
        
    }

  mediandat <- plotdata%>%
               group_by(year)%>%
              summarise(mean.n = mean(N, na.rm = TRUE),
                        sd.n = sd(N, na.rm = TRUE),
                        n.temp = n(),
                        median.n = median(N,na.rm=T),
                        max.n = max(N,na.rm=T),
                        min.n = min(N,na.rm=T)) %>%
              mutate(se.n = sd.n / sqrt(n.temp),
                     lower.ci.n = mean.n - qt(1 - (0.05 / 2), n.temp - 1) * se.n,
                     upper.ci.n = mean.n + qt(1 - (0.05 / 2), n.temp - 1) * se.n)%>%
                ungroup()%>%data.frame()
  
  p1 <- ggplot()+
  geom_line(data=plotdata,aes(x=year,y=N,group=nsim),lty=2,size=0.5,col="grey50")+
  geom_line(data=mediandat,aes(x=year,y=mean.n),size=1.2)+
  geom_line(data=mediandat,aes(x=year,y=lower.ci.n),lty=2,size=1.2)+
  geom_line(data=mediandat,aes(x=year,y=upper.ci.n),lty=2,size=1.2)+
  theme_bw()+
  labs(x="Year",y=expression(paste("Population size" %+-% "95% CI",sep="")))
  
  #Return data if path provided
  if(!is.null(path)){write.csv(x = plotdata,file = path,row.names=F)}

  return(p1)

}
