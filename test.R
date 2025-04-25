library(LakeGeoid)
library(sf)
data(kivu) #SWOT Raster data (250m) for lake Kivu 
data(myshape) # shapefile for lake kivu
data(C2kivu) # CryoSat-2 data for lake Kivu

#myshape <- st_read("myshape.shp")
#simplify shapefile to reduce the amount of triangles in the mesh 
shpsimple<-st_simplify(myshape,dTolerance=300)
fit<-getLakeGeoid(kivu,maxEdge=3000,myshape=shpsimple,UTM=TRUE)
plot(fit)
summary(fit)
print(fit)
#get geoid correction at the coordinates of CryoSat-2 
C2new<-getGeoidCorr(fit,C2kivu)
