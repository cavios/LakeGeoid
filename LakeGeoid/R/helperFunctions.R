# convert to utm coordinates





getProjUTM<-function(lon,lat){
    hem<-ifelse(mean(lat)< 0,'south','north')
    myutm<-floor((mean(lon) + 180) / 6) + 1
    myProj<-paste0("+proj=utm +zone=",myutm," +",hem," +ellps=WGS84")
    myProj
}



ll2utm<-function(dat){
    myProj<-getProjUTM(dat$lon,dat$lat)
    datLL<-cbind(dat$lon,dat$lat)
    v <- vect(datLL, crs="+proj=longlat")
    u <- terra::project(v, myProj)
    xy <- crds(u)
}

