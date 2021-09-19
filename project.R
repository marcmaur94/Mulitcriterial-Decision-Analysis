setwd("C:/MCDA/Zurich_GIS")
install.packages("rgdal","raster","tmap","reshape2","ggplot2","geoR","gstat","sf","geosphere","tiff")
install.packages('lpSolve',repos='http://cran.us.r-project.org')
library(sf)
library(gstat)
library(geoR)
library(rgdal) 
library(raster) 
library(tmap)
library(reshape2)
library(ggplot2)
library(geosphere)
library(lpSolve)
library(tiff)

#loading all necessary datasets
hospitals = read.csv("hospitals_zh_df.csv",sep=";")
pop_dens= readOGR("other_data_zh\\Population_density\\gd-b-00.03-13-vz2013statpopg\\STATPOP2013G.shp")
pop_dens_df=as.data.frame(pop_dens)
zurich_hp_shape = readOGR("zurich_data.shp")
zurich_hp_df=data.frame(zurich_hp_shape)
border_zh_shape = readOGR("vector25_zh\\Kanton_zh.shp") 
roads_zh = readOGR("vector25_zh\\hauptstrassen_zh.shp")
roads_zh_df=as.data.frame(roads_zh)
boden_shape = readOGR("vector25_zh\\boden_zh.shp")
zentren = readOGR("other_data_zh\\Kant_Richtplan\\SLA\\Zentrum_polygon.shp") #polygons of cantonal/regional centres
no2 = raster("other_data_zh\\NO2_Modellierung\\grids\\no2_immi_2010.tif")
oev= readOGR("ch.are.gueteklassen_oev\\LV03\\Oev_Gueteklassen_ARE.shp")
as09_rast = raster("arealstatistik_zh\\as09_27.tif") 
zh_hs_shape = readOGR("vector25_zh\\hauptstrassen_zh.shp")




###criterias for the MADM

##distance to highway access
highway_access=roads_zh[which(roads_zh$OBJECTVAL=="Ein_Ausf"|roads_zh$OBJECTVAL=="A_Zufahrt"),] #extracting all the highway accesses out of the roads_zh dataset
#map with the highway access polygons
highway_access_raster=rasterize(highway_access, as09_rast,"LENGTH") #The argument LENGTH doesn't really matter since we run the distance function on the layer anyway
dist_highway_access=distance(highway_access_raster)
#now we are setting the utility function such that the scores is 1 if it is below 300 meters, 0 if it is above 5000 meters and it should decrease linearly in between that.
uf_autobahn= function(x) {  
  y = vector() 
  y[x<=300] = 1
  y[x>300&x<5000] = -1/4700*x[x>300&x<5000]+(1+300/4700)
  y[x>=5000] = 0  
  return(y)  
} 
plot(uf_autobahn, xlim = c(0,5500),main= 'Utility function for highway access',ylab='Score',xlab='Distance to highway access in meters')
highway_access_criterion = calc(dist_highway_access, uf_autobahn) #calculating the scores of our layer with the given utility function
highway_access_criterion = mask(highway_access_criterion,border_zh_shape)
tm_shape(highway_access_criterion)+
  tm_raster(title="Score")+
  tm_shape(roads_zh[roads_zh$OBJECTVAL=="Ein_Ausf"|roads_zh$OBJECTVAL=="Autobahn"|roads_zh$OBJECTVAL=="Autostr"|roads_zh$OBJECTVAL=="Autob_RI"|roads_zh$OBJECTVAL=="A_Zufahrt",])+
  tm_lines()+
  tm_shape(border_zh_shape, is.master = TRUE)+
  tm_polygons(lwd = 1,alpha=0)+
  tm_layout(legend.position = c(0.75,0.71),inner.margins = c(0.03,0.03,0.03,0.08))+
  tm_scale_bar(position = c(0.525,0.02))+
  tm_compass(position = c(0.85,0.03))


##n02 emmissions
#here we are setting the utility function such that the score is 1 if the no2 measurement is below 0 (non-existent), 0 if it has the maximum score of all data and it should decrease linearly in between that.
uf_no2= function(x) {  
  y = vector() 
  y[x<0] = 1
  y[x>=0] = -2/5/21741*x[x>=0]+1
  return(y)  
} 
plot(uf_no2,xlim=c(0,21741),main= 'Utility function for NO2 emissions',ylab='Score',xlab='NO2 emissions in [ng / m3]')
crs(no2)="+init=epsg:21781" #setting the right coordinate reference system
no2_criterion=mask(no2,border_zh_shape) #clipping the layer to the cantonal borders
no2_criterion = resample(no2_criterion, as09_rast, method="bilinear") #This is the only criterion which does not have a raster size of 100x100m2. Therefore it needs to be adjusted.
no2_criterion$score = calc(no2_criterion, uf_no2) #calculating the scores of our layer with the given utility function
tm_shape(no2_criterion$score)+
  tm_raster(title = "Score")+
  tm_shape(border_zh_shape, is.master = TRUE)+
  tm_polygons(lwd = 1,alpha=0)+
  tm_layout(legend.position = c(0.75,0.71),inner.margins = c(0.03,0.03,0.03,0.08))+
  tm_scale_bar(position = c(0.525,0.02))+
  tm_compass(position = c(0.85,0.03))


##accessibility of public transport (ÖV-Güteklasse)
#for each class ("ÖV-Güteklasse") from A (best) to E (worst) a score is assigned: A is 1, B is 0.8, C is 0.6, D is 0.4 and E would be 0.2 (There are no fields with score E).
for (i in 1:nrow(oev)) { 
  if (oev$KLASSE[i]=="A"){
    oev$score[i]=1
  }
    else 
      if (oev$KLASSE[i]=="B"){
        oev$score[i]=0.75
      }
        else
          if (oev$KLASSE[i]=="C"){
        oev$score[i]=0.5
          }
        else
          if (oev$KLASSE[i]=="D"){
            oev$score[i]=0.25
          }
        else {
            oev$score[i]=0
          }
}

oev_criterion=rasterize(oev, as09_rast, "score") #making a raster out of the dataset
oev_criterion[is.na(oev_criterion)]=0 #setting empty fields to 0 
oev_criterion = mask(oev_criterion, border_zh_shape) #clipping the map to the area of the canton of Zurich
m=as.data.frame(matrix(c("A","B","C","D","E",1,0.75,0.5,0.25,0),5,2))
ggplot(m,aes(x=m$V1,y=m$V2))+
  geom_col()+
  xlab('Category of "ÖV-Güteklasse"')+
  ggtitle("Utility function for accessibility of public transport")+
  ylab("Score")+
  theme_classic() 

tm_shape(oev_criterion)+
  tm_raster(title = "Score")+
  tm_shape(border_zh_shape, is.master = TRUE)+
  tm_polygons(lwd = 1,alpha=0)+
  tm_layout(legend.position = c(0.75,0.71),inner.margins = c(0.03,0.03,0.03,0.08))+
  tm_scale_bar(position = c(0.525,0.02))+
  tm_compass(position = c(0.85,0.03))


##distance to main road
zh_hs_rast = rasterize(zh_hs_shape, as09_rast, "LENGTH") #making a raster out of the main road dataset
zh_hs_dist_rast = distance(zh_hs_rast) #calculating the distances to the main roads
zh_hs_dist_rast = mask(zh_hs_dist_rast, border_zh_shape) #clipping the map to the area of the canton of Zurich
#now we are creating a utility function which assigns the score of 1 to places closer than 500 meters to a main road and 0 to places further away. This will be used as a constraint in the MODM.
uf_roads= function(x) {  
  y = vector() 
  y[x<=1000] = 1-1/1000*x[x<=1000]
  y[x>1000] = 0  
  return(y) 
}
plot(uf_roads,xlim=c(0,1000),main= 'Utility function for the distance to main roads',ylab='Score',xlab='Distance to main street  in meters')
criterion_hs=calc(zh_hs_dist_rast,uf_roads) #calculation of the criteria dataset
crs(criterion_hs)="+init=epsg:21781" #setting the right coordinate reference system

tm_shape(criterion_hs)+
  tm_raster(title = "Score")+
  tm_shape(border_zh_shape, is.master = TRUE)+
  tm_polygons(lwd = 1,alpha=0)+
  tm_layout(legend.position = c(0.75,0.71),inner.margins = c(0.03,0.03,0.03,0.08))+
  tm_scale_bar(position = c(0.525,0.02))+
  tm_compass(position = c(0.85,0.03))


##central areas
center_rast=rasterize(zentren,as09_rast,"ZENTRUM_ID")
center_rast[center_rast>0]=1 #here, we are setting the score to 1 if the pixel is in the central area...
center_rast[is.na(center_rast)]=0 #...and to 0 if it is outside.
center_rast= mask(center_rast, border_zh_shape) #clipping the map to the area of the canton of Zurich
crs(center_rast)="+init=epsg:21781" #setting the right coordinate reference system
m=as.data.frame(matrix(c("inside","outside",1,0),2,2))
ggplot(m,aes(x=m$V1,y=m$V2))+
    geom_col()+
    xlab("")+
    ggtitle("Utility function for regional centers")+
    ylab("Score")+ 
    theme_classic() 

tm_shape(center_rast)+
  tm_raster(col="layer",title = "Score",n=2,labels = c("0","1"))+
  tm_shape(border_zh_shape, is.master = TRUE)+
  tm_polygons(lwd = 1,alpha=0)+
  tm_layout(legend.position = c(0.75,0.75),inner.margins = c(0.03,0.03,0.03,0.08))+
  tm_scale_bar(position = c(0.525,0.02))+
  tm_compass(position = c(0.85,0.03))


##housing prices
#doing kriging on existent housing price data with the price per sqm
zurich_hp_shape_rd=remove.duplicates(zurich_hp_shape)
zurich_hp_df_rd = zurich_hp_shape_rd@data
zurich_hp_df = zurich_hp_shape_rd@data 
rs = sample(c(1:nrow(zurich_hp_df_rd)), 5000, replace = FALSE) #Randomly sample 5000 rows from the rows in zurich_hp_df_rd 
rs_ts = rs[1:4000] #4000 rows are in the training set 
rs_vs = rs[4001:5000] #1000 rows are in the validation set 
zurich_hp_df_rd$PPSQM=zurich_hp_df$P/zurich_hp_df$SQM  #calculation of the price per square meter 
zurich_hp_df_ts= zurich_hp_df_rd[rs_ts,] # training set data frame 
zurich_hp_df_vs= zurich_hp_df_rd[rs_vs,] # validation 
new_rast_extent = extent(669000, 717000, 223000, 284000) #These coordinates represent the area that the raster will cover and enclose the canton of Zurich. 
new_rast = raster(ext = new_rast_extent,resolution = 500) #We set the spatial resolution of the map to 500 m raster cells. You may change this, but kriging may give an error if you have too small cells. 
crs(new_rast)="+init=epsg:21781" #Here we set the coordinate system of the new raster. 
new_rast_coords = coordinates(new_rast) #Make a dataframe with the coordinates of all the raster cells for later use
zurich_data_ts_geo=as.geodata(zurich_hp_df_ts,coords.col=c("X","Y"),data.col="PPSQM") # property rents denoted as "P" in the data set 
lags = seq(length=1000,from=100,by=100) #number and size of lags 
zurich_ts_sv=variog(coords=zurich_data_ts_geo$coords,data=zurich_data_ts_geo$data,max.dist=10000,uvec=lags) #Calculates the semivariance for each distance lag 
plot(zurich_ts_sv,main= "Semivariogram of housing prices data")+
lines.variomodel(cov.model = "exponential", cov.pars = c(750000,850), nug =400000, max.dist=10000)
kc = krige.control(type.krige="ok", cov.model="exponential",nugget=800,cov.pars=c(850000,800)) #Sets the kriging parameters, which you can copy from the function above.
zh_hp_OK=krige.conv(zurich_data_ts_geo,locations=new_rast_coords,krige=kc) #Perform the kriging on all the coordinates in the new raster.
zh_hp_OK_rast = new_rast
hist(zh_hp_OK$predict,xlim=c(10,35),main="Histogram of housing price data",xlab='Price for housing in [CHF/m2]')
values(zh_hp_OK_rast) = zh_hp_OK$predict
crs(zh_hp_OK_rast)="+init=epsg:21781" #setting the right coordinate reference system

#Following, we create a utility function for the housing prices, which assigns scores between 0 and 1 to the raster data. 
#The chosen function especially pays attention that fields with a low price have a high score and such with a high price have a low score. Therefore, the scores in the middle have a big variety.
uf_hp = function(x) {  
  y = -1/pi*atan(0.5*(x-median(zh_hp_OK$predict)))+0.5
  return(y) 
} 
plot(uf_hp,xlim=c(0,40),ylim=c(0,1),main= 'Utility function for the housing price',ylab='Score',xlab='Price for housing in [CHF/m2]')
housing_price_criterion = calc(zh_hp_OK_rast, uf_hp) #calculation of the criteria dataset
crs(zh_hp_OK_rast)="+init=epsg:21781" #setting the right coordinate reference system
housing_price_criterion = mask(housing_price_criterion,border_zh_shape)
  
tm_shape(housing_price_criterion)+
  tm_raster(title = "Score")+
  tm_shape(border_zh_shape, is.master = TRUE)+
  tm_polygons(lwd = 1,alpha=0)+
  tm_layout(legend.position = c(0.75,0.71),inner.margins = c(0.03,0.03,0.03,0.08))+
  tm_scale_bar(position = c(0.525,0.02))+
  tm_compass(position = c(0.85,0.03))


##population catchment area
pop_dens_rast=rasterize(pop_dens,as09_rast,"B13BTOT")#rasterizing the population density dataset using the Arealstatistik
size=121 #population catchment area with 12km radius => 120 cells (+1 to make it an odd number)
m=matrix(0,nrow=size,ncol=size) #creating a matrix with zeros and with the size of the catchment area
middle= ceiling(size/2) #calculating the radius of the catchment area
#Here we want to make a focal matrix where the score is 1 if the field is in a circle of 12km radius (=middle) around the centre of the matrix (focal point)
#To check if the field is in that radius or not, we use Pythagoras.
for (i in 1:nrow(m)){
  for (j in 1:ncol(m)){
    if(sqrt((middle-i)^2+(middle-j)^2)<middle){
      m[i,j]=1
      }
  }
}
catchment_area = focal(pop_dens_rast, w=m,fun = sum, na.rm=TRUE,pad=TRUE)#creating the focal statistics using our matrix created above (m) as a moving window.
hist(catchment_area,xlab="Population",main="Histogram of population catchment after focal statistics",xlim=c(0,500000))

#utility function for population density
uf_pop_dens = function(x) {  
  y = vector() 
  y[x<150000] = 1/150000*x[x<150000]
  y[x>=150000] = 1  
  return(y) 
} 
plot(uf_pop_dens,xlim=c(0,200000),main= 'Utility function for the population density',ylab='Score',xlab='population density')
pop_dens_criterion = calc(catchment_area, uf_pop_dens)
crs(zh_hp_OK_rast)="+init=epsg:21781" #setting the right coordinate reference system
pop_dens_criterion = mask(pop_dens_criterion,border_zh_shape)
tm_shape(pop_dens_criterion)+
  tm_raster(title = "Score")+
  tm_shape(border_zh_shape, is.master = TRUE)+
  tm_polygons(lwd = 1,alpha=0)+
  tm_layout(legend.position = c(0.75,0.71),inner.margins = c(0.03,0.03,0.03,0.08))+
  tm_scale_bar(position = c(0.525,0.02))+
  tm_compass(position = c(0.85,0.03))


##matching origin and resolution
housing_price_criterion = resample(housing_price_criterion, as09_rast, method="bilinear") #This is the only criterion which does not have a raster size of 100x100m2. Therefore it needs to be adjusted.
housing_price_criterion= mask(housing_price_criterion, border_zh_shape) #clipping the map to the area of the canton of Zurich



###criterias weighting
#the following matrix (criterias2) is the result of an analytical hierarchical process.
criterias2=read.csv("criterias.csv",sep=";")
sum_col=colSums(criterias2)
size=7
N=matrix(0,size,size)
#here we calculate the matrix N according to the slides.
for (i in 1:nrow(criterias2)) {
  for (j in 1:ncol(criterias2)) {
    N[i,j]=criterias2[i,j]/sum_col[j]
  }
}
x=rowSums(N)/size
criterias_new=matrix(0,size,size)
for (i in 1:nrow(criterias2)) {
  for (j in 1:ncol(criterias2)) {
    criterias_new[i,j]=criterias2[i,j]*x[j]
  }
}
means=rowSums(criterias_new)/x
lamda=mean(means)
ci=(lamda-size)/(size-1)
ri=1.32
cr=ci/ri
print(cr<0.1) #prints TRUE if the set weights for our criterias are consistent and FALSE if they are not.



###final equation to sum up all layers

total=x[1]*highway_access_criterion+x[2]*no2_criterion+x[3]*oev_criterion+x[4]*criterion_hs+x[5]*center_rast+x[6]*housing_price_criterion+x[7]*pop_dens_criterion
total=total$layer.2 #From the NO2 file, somehow we get the two layers, but we use only the second one

tm_shape(total)+
  tm_raster(title = "Score")+
  tm_layout(inner.margins = c(0.03,0.25,0.15,0.03))+
  tm_layout(legend.position = c(0.75,0.71),inner.margins = c(0.03,0.03,0.03,0.08))+
  tm_scale_bar(position = c(0.525,0.02))+
  tm_compass(position = c(0.85,0.03))


#For sensitivity analysis, we increased each criterion score by 0.1 and decreased all the other scores by 0.1/6 (to make it 1 in the sum)

total_sens_anal=(x[1]-0.1/6)*highway_access_criterion+(x[2]+0.1)*no2_criterion+(x[3]-0.1/6)*oev_criterion+(x[4]-0.1/6)*criterion_hs+(x[5]-0.1/6)*center_rast+(x[6]-0.1/6)*housing_price_criterion+(x[7]-0.1/6)*pop_dens_criterion
crs(total)="+init=epsg:21781" #setting the right coordinate reference system


#prepare hospitals dataset for mapping
x_coords=hospitals$coords.x1
y_coords=hospitals$coords.x2
hospitals_coords=cbind(x_coords,y_coords)
hospital_points=SpatialPoints(hospitals_coords[,1:2])
crs(hospital_points)="+init=epsg:21781"
scores_hospitals=extract(total,hospital_points) #extracting the scores from the raster file in our hospital locations 
hospital_points$scores=scores_hospitals #adding scores to hospital points file
hospitals$score=scores_hospitals #adding scores to hospitals data frame
hospital_points$Hospital="Hospital"

tm_shape(total)+
  tm_raster(title = "Score")+
  tm_layout(legend.position = c(0.75,0.57),inner.margins = c(0.03,0.03,0.03,0.15))+
  tm_shape(hospital_points)+
  tm_dots("score",size=0.1,col="black")+
  tm_scale_bar(position = c(0.525,0.02))+
  tm_compass(position = c(0.85,0.03))

best_locations=(total>max(hospitals$score))
tm_shape(border_zh_shape)+
  tm_polygons()+
  tm_shape(best_locations)+
  tm_raster(palette = c("FALSE"="gainsboro","TRUE"="red"))+
  tm_layout(legend.show = FALSE,inner.margins = c(0.03,0.03,0.03,0.1),title.position = c(0.03,0.96))+
  tm_scale_bar(position = c(0.525,0.02))+
  tm_compass(position = c(0.85,0.03))

tm_shape(border_zh_shape)+
  tm_polygons()+
  tm_shape(hospital_points)+
  tm_dots("scores",size=0.3,title="Hospital scores")+
  tm_layout(inner.margins = c(0.03,0.03,0.03,0.25),legend.bg.color = "lightgrey",legend.frame = "black",legend.frame.lwd = 0.3)+
  tm_scale_bar(position = c(0.525,0.02))+
  tm_compass(position = c(0.85,0.03))

hospitals$score_round=round(hospitals$score,2)
ggplot(hospitals,aes(x=name,y=score_round))+
  geom_text(label=hospitals$score_round,nudge_y = 0.1)+
  geom_bar(stat="identity",width = 0.8)+
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab(" ")+
  ylab("Score")
hist(total,xlab="Score",main="Histogram of scores",xlim=c(0,1),ylim=c(0,25000))



###Linear Programming with the hospital scores.
#With linear programming, we try to figure out which hospitals should be shut down to reduce the capacity 
#in the whole canton by 25%. Therefore, we tried 3 different approaches.

#setting variables needed for the objective function and constraints
max_beds=sum(hospitals$Betten)*0.75 #maximum number of beds in the canton
beds=hospitals$Betten #capacity of the hospitals, measured by the number of beds
beds_mat = matrix(beds, nrow = 1) #bringing the beds dataset into the right dimensions
scores=hospitals$score #scores extracted at the hospital locations in the score map


#Approach 1: the objective function maximizes the scores of the hospitals.
#Variables can only have the score 0 (being deconstructed) or 1 (staying)
#As a result, the Universitätsspital Zürich with the highest of all capacities is shut down.

output1=lp("max",scores,beds_mat,"<",max_beds,all.bin = TRUE)
hospital_points$neuverteilung=output1$solution
hospitals_choice=hospital_points[hospital_points$neuverteilung==1,]
hospitals$lp_only_scores=output1$solution

tm_shape(total)+
  tm_raster(title = "Score")+
  tm_shape(hospitals_choice,projection=crs("+init=epsg:21781"))+
  tm_dots(size=0.2)+
  tm_layout(inner.margins = c(0.02,0.02,0.12,0.12),legend.position = c(0.71,0.64),title = "Chosen hospitals considering only scores",title.position = c(0.02,0.96))

hospitals[,c("name","Betten","lp_only_scores")]


#Approach 2: the objective function maximizes an abstract score consisting of the scores of the hospitals times the number of beds.
#Variables can still only have the score 0 (being deconstructed) or 1 (staying)
#Compared to Approach 1, hospitals with little capacity but a good score gain a higher weight.
#Also, hospitals with high capacity and a bad score loose weight.

output2=lp("max",scores*beds,beds_mat,c("<"),max_beds,all.bin = TRUE)
hospital_points$neuverteilung=output2$solution
hospitals_choice=hospital_points[hospital_points$neuverteilung==1,]
hospitals$lp_capacity_and_score=output2$solution

tm_shape(total)+
  tm_raster(title = "score")+
  tm_shape(hospitals_choice,projection=crs("+init=epsg:21781"))+
  tm_dots(size=0.2)+
  tm_layout(inner.margins = c(0.02,0.02,0.12,0.12),legend.position = c(0.71,0.64),title = "Chosen hospitals considering capacity and score",title.position = c(0.02,0.96))

hospitals[,c("name","Betten","lp_capacity_and_score")]


#Approach 3: As an alternative, we try here to set a new capacity (number of beds) for each hospital.
#The objective function only maximizes the scores again, but this time the variables can be any natural number (representing the number of beds).
#As constraints, the scores can not exceed the present capacity, meaning that the current hospitals should not be extended. 
#Further, the number of beds overall has to be 75% of the existing number of beds again.

output3=lp(direction="max",scores,rbind(rep(1,26) ,diag(1,nrow(hospitals))),rep("<=",nrow(hospitals)+1),c(max_beds,beds))
hospital_points$neuverteilung=output3$solution
hospitals_choice=hospital_points[hospital_points$neuverteilung!=0,]
hospitals$lp_only_capacity=output3$solution

tm_shape(total$layer.2)+
  tm_raster(title = "Score")+
  tm_shape(hospitals_choice,projection=crs("+init=epsg:21781"))+
  tm_dots(size=0.2)+
  tm_layout(inner.margins = c(0.02,0.02,0.12,0.12),legend.position = c(0.71,0.64),title = "Chosen hospitals considering only capacity",title.position = c(0.02,0.96))

hospitals[,c("name","Betten","lp_only_capacity")]


results_overview=hospitals[,c("name","Betten","score_round","lp_only_scores","lp_capacity_and_score","lp_only_capacity")]

 