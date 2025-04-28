# LakeGeoid


The package can be installed from r-universe by:
```R
install.packages("LakeGeoid", repos = c("https://cavios.r-universe.dev", "https://cloud.r-project.org"))
```

Please find a description of all functions [here](https://cavios.r-universe.dev/LakeGeoid/doc/manual.html)


## A quick example for the African lake Kivu

The package contain 3 test data sets



```R
library(LakeGeoid)
library(sf)
data(kivu) #SWOT Raster data (250m) for lake Kivu 
data(myshape) # shapefile for lake kivu
data(C2kivu) # CryoSat-2 data for lake Kivu

#myshape <- st_read("myshape.shp")
#simplify shapefile to reduce the amount of triangles in the mesh 
shpsimple<-st_simplify(myshape,dTolerance=300)
fit<-getLakeGeoid(kivu,maxEdge=3000,myshape=shpsimple,UTM=TRUE)

summary(fit)

-------------------------
Summary for getLakeGeoid
-------------------------
Converged with a negative log likelihood of -177695.689 

Number of parameters: 4 

Par 1 : sdObs  =  0.072 

Par 2 : Tau  =  1100142845.643 

Par 3 : Kappa  =  0 

Par 4 : SigmaRW  =  2.639 

```

Plotting the model output

```R
plot(fit)
```

<p align="center">
  <img src="figs/out-1.png?raw=true">
</p>


We can get geoid correction at the coordinates of CryoSat-2 

```R
C2new<-getGeoidCorr(fit,C2kivu)
head(C2new)
      time       lat      lon   height satid ocval sat     gCorr
1 2010.559 -2.332386 29.12126 1462.901     9    97 C2E 0.1911031
2 2010.559 -2.329602 29.12361 1462.984     9    97 C2E 0.1859167
3 2010.559 -2.326954 29.12208 1462.767     9    98 C2E 0.1833408
4 2010.559 -2.318522 29.13137 1462.518     9    98 C2E 0.1705378
5 2010.559 -2.315892 29.12931 1461.234     9    97 C2E 0.1637882
6 2010.559 -2.313156 29.13028 1463.349     9    97 C2E 0.1558789
```

The column ```gCorr``` contains the geoid correction which sould be added to the heights





