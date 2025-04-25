f <- function(par){
    "[<-" <- ADoverload("[<-")
    getAll(par, data)
    sdObs<-exp(logsdObs)
    sigmaRW<-exp(logSigmaRW)
    tau <- exp(logTau)
    kappa <- exp(logKappa)
    Q <- tau*(kappa^4*spde$c0+2*kappa^2*spde$g1+spde$g2)

    ntime<-length(mymu)
    nll <- -sum(dnorm(mymu[2:ntime],mymu[1:(ntime-1)],sigmaRW*sqrt(utime[2:ntime]-utime[1:(ntime-1)]),TRUE))
    PredObs <- projObs%*%omega 
    nll<- nll-dgmrf(omega, Q=Q, log = TRUE) #hack
    
    nll <-nll -sum(dnorm(obs,PredObs+mymu[timeid],sdObs,TRUE))
  
    sdField <- sqrt(1/(4*pi*tau^2*kappa^2))
    range <- sqrt(8)/kappa
    ADREPORT(sdField)
    ADREPORT(range)
    nll
}



##' Reconstruct missing geoid signal and a water level time series
##' @param dat Input data set; must at least contain columns names: time (decimal years), lon (longitude in decimal degrees ), lat (latitude in decimal degrees), height (meters)  
##' @param maxEdge Maximum side length of the triangles in the mesh (created by the function "fm_mesh_2d" from the package "fmesher").
##' @param myshape a shapefile/polygon of the class "sf" defining the boundary of the lake where the model will be reconstructed. 
##' @param UTM logic variable to specify if the coordinates should be projected to UTM coordinates (the default is TRUE)  
##' @description The function getLakeGeoid reconstruct a static spatial water level signal and a water level time series. The static spatial signal is modeled as a Gaussian Markov Random Field on a triangular mesh where the nodes specifies the neighbor structure. The water level time series is modeled with an Random Walk as the underlying process and the observation error is here assumed to be Gaussian. The Static spatial field is here intended to models potential missing geoid model signals.       
##' @import RTMB
##' @import fmesher
##' @import terra
##' @importFrom  sf st_transform
##' @export

getLakeGeoid<-function(dat,maxEdge,myshape,UTM=TRUE){
    xy<-cbind(dat$lon,dat$lat)
    if(UTM){
        xy<-ll2utm(dat)
        myproj<-getProjUTM(dat$lon,dat$lat)
        # project shapefile
        myshape <- st_transform (myshape, crs = myproj)
    }
    
    mymesh<-fm_mesh_2d(boundary =myshape$geometry,max.edge=maxEdge)
    spde <- fm_fem(mymesh)
    
    
    projObs <- fm_basis(mymesh, loc = xy)
    #trim observations on the border of the polygon 
    idx<-which(rowSums(projObs)>0.5)
    dat<-dat[idx,]
    xy<-xy[idx,]
    projObs <- fm_basis(mymesh , loc = xy)

    timeid<-as.integer(as.factor(dat$time))
    utime<-unique(dat$time)
# data list
    data<-list()
    data$mesh <- mymesh
    data$spde <- spde
    data$projObs <- projObs
    data$obs <- dat$wse
    data$timeid<-timeid
    data$utime<-utime

# define parameters
    par <- list()
    par$mymu <- rep(mean(data$obs),length(utime))
    par$logsdObs<-log(0.2)
    par$logTau <- 0  
    par$logKappa <- 0  
    par$omega <- numeric(data$mesh$n)
    par$logSigmaRW<-0

    data <- local({data<-data; environment()})
    environment(f) <- data

    obj <- MakeADFun(f, par, random=c("omega","mymu"), silent=FALSE)
    fit <- nlminb(obj$par, obj$fn, obj$gr,control=list(eval.max=2000, iter.max=2000))

    sdr <- sdreport(obj)
    pl <- as.list(sdr,"Est")
    plsd <- as.list(sdr,"Std")
    plr <- as.list(sdr,"Est", report=TRUE)  
    plrsd <- as.list(sdr,"Std", report=TRUE)
    mymu<-pl$mymu
    omega<-pl$omega
    out<-list(omega=omega,wl=mymu, time=utime,mymesh=mymesh,fit=fit)
    class(out)<-"LakeGeoid"
    out
}



#' Plot an object returned by the function getLakeGeoid()
#' @param fit Object returned by getLakeGeoid()
#' @param zlim vector with the zlim interval c(low,high), defined by the field "omega" if not specified 
#' @param dat The raw water level data can be added to the plot
#' @param doSave plot is saved to a pdf file
#' @importFrom fields image.plot
#' @import viridis
#' @keywords plot
#' @method plot LakeGeoid
#' @examples plot(fit)
#' @export
#' 

plot.LakeGeoid <-function(fit,zlim=NULL,dat=NULL,doSave=FALSE){
    colS<-'blue3'
    mymesh<-fit$mymesh
    omega<-fit$omega
    wl<-fit$wl
    time<-fit$time
    if(is.null(zlim)){
        eps <- 1e-9
        zlim <-  c(min(omega,na.rm=TRUE)-eps, max(omega,na.rm=TRUE)+eps)

    }
    if(doSave)pdf('out.pdf',12,6)
    par(mfrow=c(1,2),mar=c(4,4,1,6))
    cc<-viridis::mako(500)
    tc<-apply(col2rgb(cc)/255, 2, function(x)rgb(x[1],x[2],x[3],1))
    plot(mymesh, col="gray", lwd=.01,xlab='Easting',ylab='Northing')
    nt<-nrow(mymesh$graph$tv)
    lamv<-apply(mymesh$graph$tv,1,function(idx)mean(omega[idx]))
    mybreaks<-seq(zlim[1],zlim[2],length=501)
    lamc<-tc[as.numeric(cut(lamv,breaks = mybreaks))]
    dummy<-sapply(1:nrow(mymesh$graph$tv),function(i){idx<-mymesh$graph$tv[i,];polygon(mymesh$loc[idx,1], mymesh$loc[idx,2], col=lamc[i], border=NA)})
    fields::image.plot(as.matrix(mybreaks), col=tc, type="n", legend.only=TRUE)
    box()
    axis(1)
    axis(2)
    par(mar=c(4,4,1,1))
    plot(time,wl,t='b',xlab='Time in decimal years',ylab='Elevation w.r.t. EGM2008',pch=3,col=colS)
    if(!is.null(dat)){points(dat$time,dat$wse,col='gray',pch='.')
        lines(time,wl,t='b',pch=3,col=colS,lwd=2)
    }

    if(doSave)dev.off()
}



##' Predict the geoid correction error   
##' @param fit Object returned by getLakeGeoid()
##' @param newdat a data set at least including the columns named "lon" and "lat" with the coordinates where the corretion should be predicted. 
##' @param UTM logic variable to indicate if the estimated field is in a UTM projection (FALSE means the field is in degrees) 
##' @import terra
##' @import sf
##' @import fmesher
##' @keywords predict
##' @export

getGeoidCorr<-function(fit,newdat,UTM=TRUE){

    omega<-fit$omega
    mymesh<-fit$mymesh
    xy<-cbind(newdat$lon,newdat$lat)
    if(UTM) xy<-ll2utm(newdat)
    projObs <- fm_basis(mymesh, loc = xy)
    PredObs <- as.vector(projObs%*%omega)
    gCorr<-rep(NA,length(PredObs))
    # identify observations outside model area
    idx<-which(rowSums(as.array(projObs))>0.5)
    gCorr[idx]<-PredObs[idx]
    newdat<-cbind(newdat,gCorr)
    
    #val <-PredObs[idx]
    #nCol <- 40
#tc <- viridis::mako(nCol)
#eps <- 1e-9
#breaks <-  seq(min(val,na.rm=TRUE)-eps, max(val,na.rm=TRUE)+eps, length=nCol+1)
#valc<-tc[as.integer(cut(val,breaks = breaks))]
#plot(C2$lon,C2$lat, col=valc, pch=16,xlab='Longitude',ylab='Latitude')
#ticks <- pretty(val)
#fields::image.plot(as.matrix(val), nlevel=nCol, col=tc, type="n", legend.only=TRUE, horizontal=FALSE, breaks=breaks)
    newdat

}




#' Simple print  of convergence statutus
#'
#' This function presents a summary of the output
#' @param x An object of class "LakeGeoid"
#' @return Print the objective function and state convergence
#' @export
#' 

print.LakeGeoid<-function(x){
    cat("tsHydro fit: negative log likelihood is",x$fit$objective, "Convergence", ifelse(x$fit$convergence==0, "OK", "failed"),"\n ")
}



#' Summary of output
#'
#' This function presents a summary of the output
#' @param x An object of class "LakeGeoid"
#' @return Summary of output
#' @method summary LakeGeoid
#' @export

summary.LakeGeoid <- function(x){
	npar <- length(x$fit$par)
	logLik <- x$fit$objective
	conv <- x$fit$convergence == 0
        mypar<-as.numeric(exp(x$fit$par))
        parnames<-sub("log","",names(x$fit$par))

        #SigmaRW<-as.numeric(exp(fit$opt$par[2]))
	#res<-list(numpar = npar,
	#	nlogLik = logLik,
	#	converged = conv,
        cat("\n-------------------------\n")
	cat("Summary for getLakeGeoid\n")
        cat("-------------------------\n")
	cat(paste(ifelse(conv,"Converged","Not converged"),"with a negative log likelihood of",round(logLik,3),"\n\n"))
	cat(paste("Number of parameters:",npar,"\n\n"))
        for(i in 1:length(mypar)){
            cat(paste("Par", i,":", parnames[i]," = ",round(mypar[i],3),"\n\n"))
        }
}



    
