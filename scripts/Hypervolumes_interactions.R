

# Name : Hypervolume analysis
# Authors : Pauline + Diane
# Objectives ## creating hypervolumes and computing distances between them
#----------------------------------------------------------------------------


# Hypervolume analyses can take a long time to run (hours - days). 

# This script is set up for several hypervolume analyses, depending on the way data
# is pooled : interactions of altitude and morphotype levels, or altitude x sex levels.
# You can find the beginning of each one by searching for "wantAnalysisByAltitude"
# (or -ByMorphotypeL or -ByMorphotypeD, etc).
# The first part of the script prepares the data and should be run in any case.





# Required packages -------------------------------------------------------

library(ggplot2)
library(GGally)
library(viridis)
library(hypervolume)
library(readxl)
library(xlsx)
library(cati)
library(mice)
library(e1071)
library(tidyverse)
library(mgcv)
library(fmsb)
library(cowplot)
library(reshape2)



# Z_test function, p_value adjusted ---------------------------------------------
ztest_function <- function(x,y){
  z_stat <- (mean(x) - mean(y)) / sqrt(sd(x)^2 / ind + sd(y)^2 / ind)
  p_value <- 2*pnorm(-abs(z_stat)) 
  p_adjust <- p.adjust(p_value, method = "bonferroni", n=num_comp)  
  return(p_adjust)
}

# defining colors manually
viridis_col=c("#39568CFF","#1F968BFF","#29AF7FFF","#73D055FF","#B8DE29FF")
viridis_col2=c("#1F968BFF","#73D055FF")
viridis_col3=c("#B8DE29FF","#39568CFF") #light in green, dark in blue

# Open dataL  -------------------------------------------------------------

DATA <- read.csv2("data/Altitudinal transects Amblystogenium pacificum.csv")

data<- na.omit(DATA)
data=data[,c(3,5:7,15,16,18,20)] #remove useless columns


# dataL preparation --------------------

# Select traits
trait_axes <- c("Mass","Body_size_index_corr","Proteines/Mass","Lipids/Mass","Sugars/Mass")
num_trait = length(trait_axes)
print(paste0("The number of considered traits is : ", num_trait))

# Remove missing values
data <- na.omit(data)

dataL=subset(data,data$Morphotype=="L")
dataD=subset(data,data$Morphotype=="D")
dataF=subset(data,data$Sex=="F")
dataM=subset(data,data$Sex=="M")

# Center and scale numerical dataL
for (i in 4:dim(dataL)[2]){ #scaling by columns [2]
  dataL[,i]<-scale(dataL[,i], center=TRUE, scale=TRUE)
}

for (i in 4:dim(dataD)[2]){ #scaling by columns [2]
  dataD[,i]<-scale(dataD[,i], center=TRUE, scale=TRUE)
}

for (i in 4:dim(dataF)[2]){ #scaling by columns [2]
  dataF[,i]<-scale(dataF[,i], center=TRUE, scale=TRUE)
}

for (i in 4:dim(dataM)[2]){ #scaling by columns [2]
  dataM[,i]<-scale(dataM[,i], center=TRUE, scale=TRUE)
}

# Hypervolumes analysis 

wantAnalysisByAltitudeL=F
if(wantAnalysisByAltitudeL==TRUE){
  
  # Check the number of individual per hypervolume (i.e. per altitude)
  table(dataL$Altitude)
  
  
  # Building hypervolumes ---------------------------------------------------
  
  # list the hypervolumes to be calculated
  cat_list_altitudeL = sort(as.character(unique(dataL$Altitude)))
  num_cat_altitudeL = length(cat_list_altitudeL)
  print(paste0("The number of altitude categories is : ", num_cat_altitudeL))
  
  # Calculate hypervolumes (1 per altitude) and join all of them
  hv_list_altitudeL = new("HypervolumeList") # list that will contain all the hypervolumes
  hv_list_altitudeL@HVList = vector(mode="list",length=num_cat_altitudeL)
  for (i in 1:num_cat_altitudeL){ 
    # select dataL
    dataL_ce_cat_altitudeL = dataL[dataL$Altitude==cat_list_altitudeL[i],trait_axes]
    # Calculate hypervolume for each altitude
    hv_list_altitudeL@HVList[[i]] <- hypervolume_gaussian(dataL_ce_cat_altitudeL,name=as.character(cat_list_altitudeL[i]),
                                                          verbose=FALSE,kde.bandwidth = estimate_bandwidth(dataL_ce_cat_altitudeL))
  }
  
  NOMpng="results/Hypervolume_Traits_by_altitude_Lmorpho.png"
  png(file = NOMpng,width=400,height = 400)
  plot(hypervolume_join(hv_list_altitudeL@HVList[[1]], hv_list_altitudeL@HVList[[2]], hv_list_altitudeL[[3]], 
                        hv_list_altitudeL@HVList[[4]], hv_list_altitudeL@HVList[[5]]),
       show.density = T, show.dataL = T, show.contour = F, cex.legend = 1.6, cex.names = 1.6, col= viridis_col, 
       names = trait_axes)
  dev.off()
  
   # Assessing similatity between hypervolumes ---------------------------------------------------
  
  # Sorensen similarity index
  overlap_altitudeL = matrix(NA, nrow=num_cat_altitudeL, ncol=num_cat_altitudeL, dimnames = list(cat_list_altitudeL,cat_list_altitudeL)) #empty matrix
  #overlap_altitudeL_CI = matrix(NA, nrow=num_cat_altitudeL, ncol=num_cat_altitudeL, dimnames = list(cat_list_altitudeL,cat_list_altitudeL)) #empty matrix
  for (i in 1:num_cat_altitudeL){
    for (j in 1:i){
      if (i!=j){
        # Assemble the hypervolumes two by two
        this_set = hypervolume_set(hv_list_altitudeL@HVList[[i]], hv_list_altitudeL@HVList[[j]], check.memory=FALSE)
        # Calculate Sorensen
        overlap_altitudeL[i,j] = hypervolume_overlap_statistics(this_set)["sorensen"]
        #overlap_altitudeL_CI[i,j] = hypervolume_overlap_confidence(this_set)["sorensen"]
      }  }  } #computing
  
  NOMpng="results/Hypervolume_Traits_by_altitude_Lmorpho_Sorensen_index.png"
  png(file = NOMpng)
  plot(hclust(as.dist(1-overlap_altitudeL)), main = "(a) Sorensen's similarity index", ylab = "1-Similarity", hang = -1)
  dev.off()
  
  
  # Centroid distances
  distance_altitudeL = matrix(NA, nrow = num_cat_altitudeL,ncol = num_cat_altitudeL, dimnames = list(cat_list_altitudeL,cat_list_altitudeL))
  for(i in 1:num_cat_altitudeL){
    for(j in 1:i){
      if (i!=j){
        distance_altitudeL[i,j]<-hypervolume_distance(hv_list_altitudeL@HVList[[i]], hv_list_altitudeL@HVList[[j]],check.memory = FALSE)
      } } }
  
  NOMpng="results/Hypervolume_Traits_by_altitude_Lmorpho_Centroid_distance.png"
  png(file = NOMpng)
  plot(hclust(as.dist(distance_altitudeL)), main = "(b) Centroid distance",ylab = "Distance", hang = -1)
  dev.off()
  
  
  # Volumes, contributions, centro?ds observed #### ---------------------------------------------------
  
  # Get volume 
  # "Volume is a measure of hypervolume size and represents the width of the multidimensional,
  # trait space of the altitudes (i.e. the variability of all traits shaping the hypervolume simultaneously)."
  volume_altitudeL<-as.data.frame(get_volume(hv_list_altitudeL))
  
  # Get profile trait contribution and plot them 
  forme_altitudeL=matrix(0,num_cat_altitudeL,num_trait, dimnames = list(cat_list_altitudeL,trait_axes))
  for (i in 1:num_cat_altitudeL){
    forme_altitudeL[i,]<-hypervolume_variable_importance(hv_list_altitudeL@HVList[[i]])
  }
  
  plotforme_altitudeL<-as.data.frame(forme_altitudeL)
  plotforme_altitudeL<- rbind(rep(round(max(plotforme_altitudeL)),num_trait) , rep(0,num_trait) ,plotforme_altitudeL) # Set plot limits
  
  NOMpng="results/Hypervolume_Traits_by_altitude_Lmorpho_profile.png"
  png(file = NOMpng)
  par(mar=rep(0.8,4))
  par(mfrow=c(2,3)) # To adapt 
  for(i in 1:num_cat_altitudeL){ 
    radarchart(plotforme_altitudeL[c(1,2,i+2),], axistype=2, title=cat_list_altitudeL[i], maxmin = T,
               plwd=3, plty=1, cglcol="grey", cglty=1, cglwd=1, axislabcol="black", vlcex=0.9)
  }
  dev.off()
  
  
  # Get centroids 
  centroides_altitudeL<-as.data.frame(get_centroid(hv_list_altitudeL))
  centroides_altitudeL$Altitude <- rownames(centroides_altitudeL)
  
  
  #  Comparing volumes (i.e. volume of each cloud) between altitudes ---------------------------------------------------
  
  # Bootstraping volumes within each altitude to compare them
  
  # "Because we could only calculate 1 hypervolume/altitude,
  # (hypervolume calculation requires at least one individual per each trait involved in hypervolume shaping), 
  # we lacked variability for comparing volumes between levels.
  # To solve this problem, we simulated 100 hypervolumes per altitude using bootstrapping"
  
  table(dataL$Altitude) # see number of individuals by altitude
  sim=100 # Number of simulation
  ind=9 # Number of individual = Number of individuals for the case where there are the least
  # "Because the volume depends on the number of values used to build the hypervolume (Lamanna et al.,2014), 
  # we set the number of values used to build the hypervolume (i.e. number of individuals) at 21 for each simulation, 
  # which was the number of values obtained in the considered factor with the least sampled level
  # (i.e., no-competition)."
  
  volumeBS_altitudeL=matrix(0,num_cat_altitudeL,sim) # Row = Altitude Col = resampling
  rownames(volumeBS_altitudeL)=cat_list_altitudeL
  for (i in 1:num_cat_altitudeL){  
    dataL_ce_cat_altitudeL = dataL[dataL$Altitude==cat_list_altitudeL[i],trait_axes] # dataL for this altitude
    for(j in 1:sim){
      
      dataLBS_altitudeL = dataL_ce_cat_altitudeL[sample(nrow(dataL_ce_cat_altitudeL),ind,replace=T),] # sampling within the altitude
      hyperBS_altitudeL = hypervolume_gaussian(dataLBS_altitudeL,verbose=FALSE, kde.bandwidth = estimate_bandwidth(dataLBS_altitudeL))
      volumeBS_altitudeL[i,j]<-get_volume(hyperBS_altitudeL)
    }}
  
  ## Boxplot of resampled (i.e. boostraped) volumes 
#   "we compared the volume of different sites, with all the hypervolumes built from the same number of values. 
#    we compared the variability of the hypervolumes built for each altitude against the variability of a hypervolume built using individuals from
# all sites at the same time"
  
  plotvolumeBS_altitudeL<-as.data.frame(volumeBS_altitudeL)
  plotvolumeBS_altitudeL$Altitude <- rownames(plotvolumeBS_altitudeL)
  plotvolumeBS_altitudeL<-melt(plotvolumeBS_altitudeL,  id="Altitude")
  
  NOMpng="results/Hypervolume_Traits_by_altitude_Lmorpho_volumeBS.png"
  png(file = NOMpng,width = 400, height =300)
  ggplot(plotvolumeBS_altitudeL, aes(Altitude, value, fill=Altitude))+
    geom_boxplot() + scale_fill_viridis_d(begin=0.3,end=0.9)+ theme_bw()
  dev.off()
  
   
  ## Ztest to test significant difference between volumes
  ztest_volume_altitudeL <- matrix(NA,num_cat_altitudeL,num_cat_altitudeL, dimnames = list(cat_list_altitudeL,cat_list_altitudeL))
  num_comp = num_cat_altitudeL*(num_cat_altitudeL-1)/2 # Number of comparisons 
  for(i in 1:num_cat_altitudeL){ 
    for(j in 1:i){
      if (i!=j){
        ztest_volume_altitudeL[i,j]<- ztest_function(plotvolumeBS_altitudeL[plotvolumeBS_altitudeL$Altitude==cat_list_altitudeL[i],"value"],
                                                     plotvolumeBS_altitudeL[plotvolumeBS_altitudeL$Altitude==cat_list_altitudeL[j],"value"])
      } } }
  
  # Plot resampled volumes (i.e. boostraped volumes) : mean, sd, letters for differences
  volumeztest_altitudeL<-as.data.frame(tapply(plotvolumeBS_altitudeL$value, plotvolumeBS_altitudeL$Altitude, mean))
  colnames(volumeztest_altitudeL)="mean"
  volumeztest_altitudeL$sd<-tapply(plotvolumeBS_altitudeL$value, plotvolumeBS_altitudeL$Altitude, sd)
  volumeztest_altitudeL$Altitude<-rownames(volumeztest_altitudeL)
  volumeztest_altitudeL$lowersd <- volumeztest_altitudeL$mean-volumeztest_altitudeL$sd
  volumeztest_altitudeL$uppersd <- volumeztest_altitudeL$mean+volumeztest_altitudeL$sd
  volumeztest_altitudeL$letter <- c("a","a","b","a","a")## Put the letters manually according to the significant differences in ztest_volume_altitudeL
  
  NOMpng="results/Hypervolume_Traits_by_altitude_Lmorpho_volumeZtest.png"
  png(file = NOMpng,width = 400, height =300)
  ggplot(volumeztest_altitudeL,aes(x= Altitude, y=mean, color=Altitude,label = letter)) +
    geom_point(shape  = 16, size   = 6) + 
    geom_errorbar(aes(ymin  = lowersd, ymax  =  uppersd), width =  0.1,size  =  1.2) +
    theme_bw() + ylab(expression(Volume~~(SD^5))) + ggtitle ("") +
    scale_color_viridis_d(begin=0.3,end=0.9) + 
    geom_text(aes(y = 150), color  = "black") +
    ylim(0,500)+
    theme(text= element_text(size=18), axis.title.x = element_blank(), 
          plot.title = element_text(size=18, face="italic"),legend.position = "none")
  
  dev.off()
  
  
  
  # Comparing centroids between altitude ####---------------------------------------------------
  
  # Boostraping centroids
  table(dataL$Altitude) 
  sim=100
  indi=9
  
  centroBS_altitudeL=matrix(0,sim,num_trait) # Resampled hypervolume centro?ds Row = simulation, Col = traits
  formeBS_altitudeL=matrix(0,sim,num_trait)  # Resampled hypervolume trait contributions Row = simulation, Col = traits
  for(j in 1:sim){
    
    dataLBS_altitudeL = dataL[sample(nrow(dataL),indi,replace=F),trait_axes] 
    hyperBS_altitudeL = hypervolume_gaussian(dataLBS_altitudeL,verbose=FALSE, kde.bandwidth = estimate_bandwidth(dataLBS_altitudeL))
    centroBS_altitudeL[j,(1:num_trait)]<-get_centroid(hyperBS_altitudeL)
    formeBS_altitudeL[j,(1:num_trait)] <- hypervolume_variable_importance(hyperBS_altitudeL) 
  }
  
  
  #### Comparison of centroids between altitudes and the null model over all altitudes ####
  
  #Create a matrix containing the lower bounds, a matrix containing the upper bounds of the CI and 
  # then compare with the observed values
  Centroides_bornes_inf_altitudeL=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Centroides_bornes_inf_altitudeL[1,i]<-quantile(centroBS_altitudeL[,i], 0.025)
  }
  
  Centroides_bornes_sup_altitudeL=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Centroides_bornes_sup_altitudeL[1,i]<-quantile(centroBS_altitudeL[,i], 0.975)
  }
  
  # Plot the results
  NOMpng="results/Hypervolume_Traits_by_altitude_Lmorpho_plotgrid_centroides.png"
  png(file = NOMpng)
  plot_grid(
    ggplot(centroides_altitudeL, aes(x = Altitude, y = Mass, color=Altitude))+
      ggtitle("") + ylab("Mass")  + xlab("Altitude") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_altitudeL[1,1], 
                    ymax = Centroides_bornes_sup_altitudeL[1,1]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.3,end=0.9) + theme_bw() +
      theme(axis.title.x = element_blank(), legend.position = "none", plot.title = element_text(size=12, face="bold.italic")),
    
    ggplot(centroides_altitudeL, aes(x = Altitude, y = Body_size_index_corr, color=Altitude))+
      ggtitle("") + ylab("Body_size_index")  + xlab("Altitude") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_altitudeL[1,2], 
                    ymax = Centroides_bornes_sup_altitudeL[1,2]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.3,end=0.9) + theme_bw() +
      theme(axis.title.x = element_blank(), legend.position = "none", plot.title = element_text(size=12, face="bold.italic")),
    
    ggplot(centroides_altitudeL, aes(x = Altitude, y = Proteines.Mass, color=Altitude))+
      ggtitle("")+ ylab("Proteins:Mass") + xlab("Altitude") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_altitudeL[1,3], 
                    ymax = Centroides_bornes_sup_altitudeL[1,3]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.3,end=0.9) + theme_bw() +
      theme(axis.title.x = element_blank(),legend.position = "none"),
    
    ggplot(centroides_altitudeL, aes(x = Altitude, y = Lipids.Mass, color=Altitude))+
      ggtitle("")+ ylab("Lipids:Mass") + xlab("Altitude") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_altitudeL[1,4], 
                    ymax = Centroides_bornes_sup_altitudeL[1,4]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.3,end=0.9)+ theme_bw() +
      theme(axis.title.x = element_blank(),legend.position = "none"),
    
    ggplot(centroides_altitudeL, aes(x = Altitude, y = Sugars.Mass, color=Altitude))+
      ggtitle("")+ ylab("Sugars:Mass") + xlab("Altitude") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_altitudeL[1,5], 
                    ymax = Centroides_bornes_sup_altitudeL[1,5]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.3,end=0.9)+ theme_bw() +
      theme(axis.title.x = element_blank(),legend.position = "none"),
    
    ncol = 3, nrow = 2, byrow = F)
  dev.off()
  

  
    ### Comparison of trait contributions between altitudes and the null model over all altitudes #### ---------------------------------------------------------------------
  
  ##Create a matrix containing the lower bounds, a matrix containing the upper bounds of the CI and then compare with the observed values
  
  Forme_bornes_inf_altitudeL=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Forme_bornes_inf_altitudeL[1,i]<-quantile(formeBS_altitudeL[,i], 0.025)
  }
  
  Forme_bornes_sup_altitudeL=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Forme_bornes_sup_altitudeL[1,i]<-quantile(formeBS_altitudeL[,i], 0.975)
  }
  
  Forme_IC_altitudeL=matrix(NA,num_cat_altitudeL,num_trait, dimnames = list(cat_list_altitudeL,trait_axes))
  for (i in 1:num_cat_altitudeL){
    for (j in 1:num_trait){
      Forme_IC_altitudeL[i,j]<- if (forme_altitudeL[i,j]<Forme_bornes_inf_altitudeL[1,j]){"-"} 
      else if (forme_altitudeL[i,j]>Forme_bornes_sup_altitudeL[1,j]){"+"} else {"ns"}
    }}
  
  # There are certainly simpler than all these lines but it works...
  plotformeIC_altitudeL<-as.data.frame(forme_altitudeL)
  plotformeIC_altitudeL$Altitude<-rownames(plotformeIC_altitudeL)
  plotformeIC_altitudeL$type<-rep("obs",num_cat_altitudeL)
  Forme_IC_altitudeL=as.data.frame(Forme_IC_altitudeL)
  Forme_IC_altitudeL$Altitude<-rownames(Forme_IC_altitudeL)
  Forme_bornes_inf_altitudeL=as.data.frame(Forme_bornes_inf_altitudeL)
  Forme_bornes_inf_altitudeL$Altitude<-rownames(Forme_bornes_inf_altitudeL)
  Forme_bornes_inf_altitudeL$type<-"inf"
  Forme_bornes_sup_altitudeL=as.data.frame(Forme_bornes_sup_altitudeL)
  Forme_bornes_sup_altitudeL$Altitude<-rownames(Forme_bornes_sup_altitudeL)
  Forme_bornes_sup_altitudeL$type<-"sup"
  plotformeIC_altitudeL<-bind_rows(plotformeIC_altitudeL,Forme_bornes_inf_altitudeL,Forme_bornes_sup_altitudeL)
  plotformeIC_altitudeL <-rbind(plotformeIC_altitudeL, rep(0,num_trait)) # Limit inf
  plotformeIC_altitudeL <-rbind(plotformeIC_altitudeL, rep(round(max(plotformeIC_altitudeL[,1:num_trait])),num_trait)) 
  
  colnames(plotformeIC_altitudeL)<- c(trait_axes,"Altitude","type")
  
  NOMpng="results/Hypervolume_Traits_by_altitude_Lmorpho_radarchart.png"
  png(file = NOMpng)
  par(mar=rep(0.8,4))
  par(mfrow=c(1,1)) #Adapt
  coul <- c("grey", "grey", viridis_col) 
  coul_leg <- c("white", "white", viridis_col)
  coul_in <- c("grey","white",rep(alpha("white",0.1),num_cat_altitudeL))
  radarchart(plotformeIC_altitudeL[c(9,8,7,6,1:5),1:num_trait], axistype=6, title="", maxmin = T,
             pcol=coul, plwd=2, plty=1, pfcol = coul_in, cglcol="grey", cglty=2, cglwd=1, axislabcol="black", vlcex=1)
  legend(x="bottomleft", legend = c("","",cat_list_altitudeL), 
         bty = "n", pch=20 , col= coul_leg , text.col = "black", cex=1, pt.cex=3, horiz = F)
  
  dev.off()
  
 
wantAnalysisByAltitudeD=F
if(wantAnalysisByAltitudeD==TRUE){
  
  # Check the number of individual per hypervolume (i.e. per altitude)
  table(dataD$Altitude)
  
  
  # Building hypervolumes ---------------------------------------------------
  
  # list the hypervolumes to be calculated
  cat_list_altitudeD = sort(as.character(unique(dataD$Altitude)))
  num_cat_altitudeD = length(cat_list_altitudeD)
  print(paste0("The number of altitude categories is : ", num_cat_altitudeD))
  
  # Calculate hypervolumes (1 per altitude) and join all of them
  hv_list_altitudeD = new("HypervolumeList") # list that will contain all the hypervolumes
  hv_list_altitudeD@HVList = vector(mode="list",length=num_cat_altitudeD)
  for (i in 1:num_cat_altitudeD){ 
    # select dataD
    dataD_ce_cat_altitudeD = dataD[dataD$Altitude==cat_list_altitudeD[i],trait_axes]
    # Calculate hypervolume for each altitude
    hv_list_altitudeD@HVList[[i]] <- hypervolume_gaussian(dataD_ce_cat_altitudeD,name=as.character(cat_list_altitudeD[i]),
                                                          verbose=FALSE,kde.bandwidth = estimate_bandwidth(dataD_ce_cat_altitudeD))
  }
  
  NOMpng="results/Hypervolume_Traits_by_altitude_Dmorpho.png"
  png(file = NOMpng,width=400,height = 400)
  plot(hypervolume_join(hv_list_altitudeD@HVList[[1]], hv_list_altitudeD@HVList[[2]], hv_list_altitudeD[[3]], 
                        hv_list_altitudeD@HVList[[4]], hv_list_altitudeD@HVList[[5]]),
       show.density = T, show.dataD = T, show.contour = F, cex.legend = 1.6, cex.names = 1.6, col= viridis_col, 
       names = trait_axes)
  dev.off()
  
 
  # Assessing similatity between hypervolumes ---------------------------------------------------
  
  # Sorensen similarity index
  overlap_altitudeD = matrix(NA, nrow=num_cat_altitudeD, ncol=num_cat_altitudeD, dimnames = list(cat_list_altitudeD,cat_list_altitudeD)) #empty matrix
  #overlap_altitudeD_CI = matrix(NA, nrow=num_cat_altitudeD, ncol=num_cat_altitudeD, dimnames = list(cat_list_altitudeD,cat_list_altitudeD)) #empty matrix
  for (i in 1:num_cat_altitudeD){
    for (j in 1:i){
      if (i!=j){
        # Assemble the hypervolumes two by two
        this_set = hypervolume_set(hv_list_altitudeD@HVList[[i]], hv_list_altitudeD@HVList[[j]], check.memory=FALSE)
        # Calculate Sorensen
        overlap_altitudeD[i,j] = hypervolume_overlap_statistics(this_set)["sorensen"]
        #overlap_altitudeD_CI[i,j] = hypervolume_overlap_confidence(this_set)["sorensen"]
      }  }  } #computing
  
  NOMpng="results/Hypervolume_Traits_by_altitude_Dmorpho_Sorensen_index.png"
  png(file = NOMpng)
  plot(hclust(as.dist(1-overlap_altitudeD)), main = "(a) Sorensen's similarity index", ylab = "1-Similarity", hang = -1)
  dev.off()
  

  
  # Centroid distances
  distance_altitudeD = matrix(NA, nrow = num_cat_altitudeD,ncol = num_cat_altitudeD, dimnames = list(cat_list_altitudeD,cat_list_altitudeD))
  for(i in 1:num_cat_altitudeD){
    for(j in 1:i){
      if (i!=j){
        distance_altitudeD[i,j]<-hypervolume_distance(hv_list_altitudeD@HVList[[i]], hv_list_altitudeD@HVList[[j]],check.memory = FALSE)
      } } }
  
  NOMpng="results/Hypervolume_Traits_by_altitude_Dmorpho_Centroid_distance.png"
  png(file = NOMpng)
  plot(hclust(as.dist(distance_altitudeD)), main = "(b) Centroid distance",ylab = "Distance", hang = -1)
  dev.off()
  
  # Volumes, contributions, centro?ds observed #### ---------------------------------------------------
  
  # Get volume 
  # "Volume is a measure of hypervolume size and represents the width of the multidimensional,
  # trait space of the altitudes (i.e. the variability of all traits shaping the hypervolume simultaneously)."
  volume_altitudeD<-as.data.frame(get_volume(hv_list_altitudeD))
  
  # Get profile trait contribution and plot them 
  forme_altitudeD=matrix(0,num_cat_altitudeD,num_trait, dimnames = list(cat_list_altitudeD,trait_axes))
  for (i in 1:num_cat_altitudeD){
    forme_altitudeD[i,]<-hypervolume_variable_importance(hv_list_altitudeD@HVList[[i]])
  }
  
  plotforme_altitudeD<-as.data.frame(forme_altitudeD)
  plotforme_altitudeD<- rbind(rep(round(max(plotforme_altitudeD)),num_trait) , rep(0,num_trait) ,plotforme_altitudeD) # Set plot limits
  
  NOMpng="results/Hypervolume_Traits_by_altitude_Dmorpho_profile.png"
  png(file = NOMpng)
  par(mar=rep(0.8,4))
  par(mfrow=c(2,3)) # To adapt 
  for(i in 1:num_cat_altitudeD){ 
    radarchart(plotforme_altitudeD[c(1,2,i+2),], axistype=2, title=cat_list_altitudeD[i], maxmin = T,
               plwd=3, plty=1, cglcol="grey", cglty=1, cglwd=1, axislabcol="black", vlcex=0.9)
  }
  dev.off()
  

  
  # Get centroids 
  centroides_altitudeD<-as.data.frame(get_centroid(hv_list_altitudeD))
  centroides_altitudeD$Altitude <- rownames(centroides_altitudeD)
  
  
  #  Comparing volumes (i.e. volume of each cloud) between altitudes ---------------------------------------------------
  
  # Bootstraping volumes within each altitude to compare them
  
  # "Because we could only calculate 1 hypervolume/altitude,
  # (hypervolume calculation requires at least one individual per each trait involved in hypervolume shaping), 
  # we lacked variability for comparing volumes between levels.
  # To solve this problem, we simulated 100 hypervolumes per altitude using bootstrapping"
  
  table(dataD$Altitude) # see number of individuals by altitude
  sim=100 # Number of simulation
  ind=11 # Number of individual = Number of individuals for the case where there are the least
  # "Because the volume depends on the number of values used to build the hypervolume (Lamanna et al.,2014), 
  # we set the number of values used to build the hypervolume (i.e. number of individuals) at 21 for each simulation, 
  # which was the number of values obtained in the considered factor with the least sampled level
  # (i.e., no-competition)."
  
  volumeBS_altitudeD=matrix(0,num_cat_altitudeD,sim) # Row = Altitude Col = resampling
  rownames(volumeBS_altitudeD)=cat_list_altitudeD
  for (i in 1:num_cat_altitudeD){  
    dataD_ce_cat_altitudeD = dataD[dataD$Altitude==cat_list_altitudeD[i],trait_axes] # dataD for this altitude
    for(j in 1:sim){
      
      dataDBS_altitudeD = dataD_ce_cat_altitudeD[sample(nrow(dataD_ce_cat_altitudeD),ind,replace=T),] # sampling within the altitude
      hyperBS_altitudeD = hypervolume_gaussian(dataDBS_altitudeD,verbose=FALSE, kde.bandwidth = estimate_bandwidth(dataDBS_altitudeD))
      volumeBS_altitudeD[i,j]<-get_volume(hyperBS_altitudeD)
    }}
  
  ## Boxplot of resampled (i.e. boostraped) volumes 
#   "we compared the volume of different sites, with all the hypervolumes built from the same number of values. 
#    we compared the variability of the hypervolumes built for each altitude against the variability of a hypervolume built using individuals from
# all sites at the same time"
  
  plotvolumeBS_altitudeD<-as.data.frame(volumeBS_altitudeD)
  plotvolumeBS_altitudeD$Altitude <- rownames(plotvolumeBS_altitudeD)
  plotvolumeBS_altitudeD<-melt(plotvolumeBS_altitudeD,  id="Altitude")
  
  NOMpng="results/Hypervolume_Traits_by_altitude_Dmorpho_volumeBS.png"
  png(file = NOMpng,width = 400, height =300)
  ggplot(plotvolumeBS_altitudeD, aes(Altitude, value, fill=Altitude))+
    geom_boxplot() + scale_fill_viridis_d(begin=0.3,end=0.9)+ theme_bw()
  dev.off()
  
 
  ## Ztest to test significant difference between volumes
  ztest_volume_altitudeD <- matrix(NA,num_cat_altitudeD,num_cat_altitudeD, dimnames = list(cat_list_altitudeD,cat_list_altitudeD))
  num_comp = num_cat_altitudeD*(num_cat_altitudeD-1)/2 # Number of comparisons 
  for(i in 1:num_cat_altitudeD){ 
    for(j in 1:i){
      if (i!=j){
        ztest_volume_altitudeD[i,j]<- ztest_function(plotvolumeBS_altitudeD[plotvolumeBS_altitudeD$Altitude==cat_list_altitudeD[i],"value"],
                                                     plotvolumeBS_altitudeD[plotvolumeBS_altitudeD$Altitude==cat_list_altitudeD[j],"value"])
      } } }
  
  # Plot resampled volumes (i.e. boostraped volumes) : mean, sd, letters for differences
  volumeztest_altitudeD<-as.data.frame(tapply(plotvolumeBS_altitudeD$value, plotvolumeBS_altitudeD$Altitude, mean))
  colnames(volumeztest_altitudeD)="mean"
  volumeztest_altitudeD$sd<-tapply(plotvolumeBS_altitudeD$value, plotvolumeBS_altitudeD$Altitude, sd)
  volumeztest_altitudeD$Altitude<-rownames(volumeztest_altitudeD)
  volumeztest_altitudeD$lowersd <- volumeztest_altitudeD$mean-volumeztest_altitudeD$sd
  volumeztest_altitudeD$uppersd <- volumeztest_altitudeD$mean+volumeztest_altitudeD$sd
  volumeztest_altitudeD$letter <- c("ab","ab","a","b","ab")## Put the letters manually according to the significant differences in ztest_volume_altitudeD
  
  NOMpng="results/Hypervolume_Traits_by_altitude_Dmorpho_volumeZtest.png"
  png(file = NOMpng,width = 400, height =300)
  ggplot(volumeztest_altitudeD,aes(x= Altitude, y=mean, color=Altitude,label = letter)) +
    geom_point(shape  = 16, size   = 6) + 
    geom_errorbar(aes(ymin  = lowersd, ymax  =  uppersd), width =  0.1,size  =  1.2) +
    theme_bw() + ylab(expression(Volume~~(SD^5))) + ggtitle ("") +
    scale_color_viridis_d(begin=0.3,end=0.9) + 
    geom_text(aes(y = 150), color  = "black") +
    ylim(0,500)+
    theme(text= element_text(size=18), axis.title.x = element_blank(), 
          plot.title = element_text(size=18, face="italic"),legend.position = "none")
  
  dev.off()
  
  
  
  
  # Comparing centroids between altitude ####---------------------------------------------------
  
  # Boostraping centroids
  table(dataD$Altitude) 
  sim=100
  indi=11
  
  centroBS_altitudeD=matrix(0,sim,num_trait) # Resampled hypervolume centro?ds Row = simulation, Col = traits
  formeBS_altitudeD=matrix(0,sim,num_trait)  # Resampled hypervolume trait contributions Row = simulation, Col = traits
  for(j in 1:sim){
    
    dataDBS_altitudeD = dataD[sample(nrow(dataD),indi,replace=F),trait_axes] 
    hyperBS_altitudeD = hypervolume_gaussian(dataDBS_altitudeD,verbose=FALSE, kde.bandwidth = estimate_bandwidth(dataDBS_altitudeD))
    centroBS_altitudeD[j,(1:num_trait)]<-get_centroid(hyperBS_altitudeD)
    formeBS_altitudeD[j,(1:num_trait)] <- hypervolume_variable_importance(hyperBS_altitudeD) 
  }
  
  
  #### Comparison of centroids between altitudes and the null model over all altitudes ####
  
  #Create a matrix containing the lower bounds, a matrix containing the upper bounds of the CI and 
  # then compare with the observed values
  Centroides_bornes_inf_altitudeD=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Centroides_bornes_inf_altitudeD[1,i]<-quantile(centroBS_altitudeD[,i], 0.025)
  }
  
  Centroides_bornes_sup_altitudeD=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Centroides_bornes_sup_altitudeD[1,i]<-quantile(centroBS_altitudeD[,i], 0.975)
  }
  
  # Plot the results
  NOMpng="results/Hypervolume_Traits_by_altitude_Dmorpho_plotgrid_centroides.png"
  png(file = NOMpng)
  plot_grid(
    ggplot(centroides_altitudeD, aes(x = Altitude, y = Mass, color=Altitude))+
      ggtitle("") + ylab("Mass")  + xlab("Altitude") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_altitudeD[1,1], 
                    ymax = Centroides_bornes_sup_altitudeD[1,1]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.3,end=0.9) + theme_bw() +
      theme(axis.title.x = element_blank(), legend.position = "none", plot.title = element_text(size=12, face="bold.italic")),
    
    ggplot(centroides_altitudeD, aes(x = Altitude, y = Body_size_index_corr, color=Altitude))+
      ggtitle("") + ylab("Body_size_index")  + xlab("Altitude") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_altitudeD[1,2], 
                    ymax = Centroides_bornes_sup_altitudeD[1,2]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.3,end=0.9) + theme_bw() +
      theme(axis.title.x = element_blank(), legend.position = "none", plot.title = element_text(size=12, face="bold.italic")),
    
    ggplot(centroides_altitudeD, aes(x = Altitude, y = Proteines.Mass, color=Altitude))+
      ggtitle("")+ ylab("Proteins:Mass") + xlab("Altitude") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_altitudeD[1,3], 
                    ymax = Centroides_bornes_sup_altitudeD[1,3]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.3,end=0.9) + theme_bw() +
      theme(axis.title.x = element_blank(),legend.position = "none"),
    
    ggplot(centroides_altitudeD, aes(x = Altitude, y = Lipids.Mass, color=Altitude))+
      ggtitle("")+ ylab("Lipids:Mass") + xlab("Altitude") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_altitudeD[1,4], 
                    ymax = Centroides_bornes_sup_altitudeD[1,4]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.3,end=0.9)+ theme_bw() +
      theme(axis.title.x = element_blank(),legend.position = "none"),
    
    ggplot(centroides_altitudeD, aes(x = Altitude, y = Sugars.Mass, color=Altitude))+
      ggtitle("")+ ylab("Sugars:Mass") + xlab("Altitude") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_altitudeD[1,5], 
                    ymax = Centroides_bornes_sup_altitudeD[1,5]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.3,end=0.9)+ theme_bw() +
      theme(axis.title.x = element_blank(),legend.position = "none"),
    
    ncol = 3, nrow = 2, byrow = F)
  dev.off()
  

  
  ### Comparison of trait contributions between altitudes and the null model over all altitudes #### ---------------------------------------------------------------------
  
  ##Create a matrix containing the lower bounds, a matrix containing the upper bounds of the CI and then compare with the observed values
  
  Forme_bornes_inf_altitudeD=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Forme_bornes_inf_altitudeD[1,i]<-quantile(formeBS_altitudeD[,i], 0.025)
  }
  
  Forme_bornes_sup_altitudeD=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Forme_bornes_sup_altitudeD[1,i]<-quantile(formeBS_altitudeD[,i], 0.975)
  }
  
  Forme_IC_altitudeD=matrix(NA,num_cat_altitudeD,num_trait, dimnames = list(cat_list_altitudeD,trait_axes))
  for (i in 1:num_cat_altitudeD){
    for (j in 1:num_trait){
      Forme_IC_altitudeD[i,j]<- if (forme_altitudeD[i,j]<Forme_bornes_inf_altitudeD[1,j]){"-"} 
      else if (forme_altitudeD[i,j]>Forme_bornes_sup_altitudeD[1,j]){"+"} else {"ns"}
    }}
  
  # There are certainly simpler than all these lines but it works...
  plotformeIC_altitudeD<-as.data.frame(forme_altitudeD)
  plotformeIC_altitudeD$Altitude<-rownames(plotformeIC_altitudeD)
  plotformeIC_altitudeD$type<-rep("obs",num_cat_altitudeD)
  Forme_IC_altitudeD=as.data.frame(Forme_IC_altitudeD)
  Forme_IC_altitudeD$Altitude<-rownames(Forme_IC_altitudeD)
  Forme_bornes_inf_altitudeD=as.data.frame(Forme_bornes_inf_altitudeD)
  Forme_bornes_inf_altitudeD$Altitude<-rownames(Forme_bornes_inf_altitudeD)
  Forme_bornes_inf_altitudeD$type<-"inf"
  Forme_bornes_sup_altitudeD=as.data.frame(Forme_bornes_sup_altitudeD)
  Forme_bornes_sup_altitudeD$Altitude<-rownames(Forme_bornes_sup_altitudeD)
  Forme_bornes_sup_altitudeD$type<-"sup"
  plotformeIC_altitudeD<-bind_rows(plotformeIC_altitudeD,Forme_bornes_inf_altitudeD,Forme_bornes_sup_altitudeD)
  plotformeIC_altitudeD <-rbind(plotformeIC_altitudeD, rep(0,num_trait)) # Limit inf
  plotformeIC_altitudeD <-rbind(plotformeIC_altitudeD, rep(round(max(plotformeIC_altitudeD[,1:num_trait])),num_trait)) 
  
  colnames(plotformeIC_altitudeD)<- c(trait_axes,"Altitude","type")
  
  NOMpng="results/Hypervolume_Traits_by_altitude_Dmorpho_radarchart.png"
  png(file = NOMpng)
  par(mar=rep(0.8,4))
  par(mfrow=c(1,1)) #Adapt
  coul <- c("grey", "grey", viridis_col) 
  coul_leg <- c("white", "white", viridis_col)
  coul_in <- c("grey","white",rep(alpha("white",0.1),num_cat_altitudeD))
  radarchart(plotformeIC_altitudeD[c(9,8,7,6,1:5),1:num_trait], axistype=6, title="", maxmin = T,
             pcol=coul, plwd=2, plty=1, pfcol = coul_in, cglcol="grey", cglty=2, cglwd=1, axislabcol="black", vlcex=1)
  legend(x="bottomleft", legend = c("","",cat_list_altitudeD), 
         bty = "n", pch=20 , col= coul_leg , text.col = "black", cex=1, pt.cex=3, horiz = F)
  
  dev.off()
  
  save.image(paste0(save_path,"insect_hypervolumes_altitude_Dmorpho.Rdata"))
}


wantAnalysisByAltitudeF=F
if(wantAnalysisByAltitudeF==TRUE){
  
  # Check the number of individual per hypervolume (i.e. per altitude)
  table(dataF$Altitude)
  
  
  # Building hypervolumes ---------------------------------------------------
  
  # list the hypervolumes to be calculated
  cat_list_altitudeF = sort(as.character(unique(dataF$Altitude)))
  num_cat_altitudeF = length(cat_list_altitudeF)
  print(paste0("The number of altitude categories is : ", num_cat_altitudeF))
  
  # Calculate hypervolumes (1 per altitude) and join all of them
  hv_list_altitudeF = new("HypervolumeList") # list that will contain all the hypervolumes
  hv_list_altitudeF@HVList = vector(mode="list",length=num_cat_altitudeF)
  for (i in 1:num_cat_altitudeF){ 
    # select dataF
    dataF_ce_cat_altitudeF = dataF[dataF$Altitude==cat_list_altitudeF[i],trait_axes]
    # Calculate hypervolume for each altitude
    hv_list_altitudeF@HVList[[i]] <- hypervolume_gaussian(dataF_ce_cat_altitudeF,name=as.character(cat_list_altitudeF[i]),
                                                          verbose=FALSE,kde.bandwidth = estimate_bandwidth(dataF_ce_cat_altitudeF))
  }
  
  NOMpng="results/Hypervolume_Traits_by_altitude_Fsex.png"
  png(file = NOMpng,width=400,height = 400)
  plot(hypervolume_join(hv_list_altitudeF@HVList[[1]], hv_list_altitudeF@HVList[[2]], hv_list_altitudeF[[3]], 
                        hv_list_altitudeF@HVList[[4]], hv_list_altitudeF@HVList[[5]]),
       show.density = T, show.dataF = T, show.contour = F, cex.legend = 1.6, cex.names = 1.6, col= viridis_col, 
       names = trait_axes)
  dev.off()
  
  # Assessing similatity between hypervolumes ---------------------------------------------------
  
  # Sorensen similarity index
  overlap_altitudeF = matrix(NA, nrow=num_cat_altitudeF, ncol=num_cat_altitudeF, dimnames = list(cat_list_altitudeF,cat_list_altitudeF)) #empty matrix
  #overlap_altitudeF_CI = matrix(NA, nrow=num_cat_altitudeF, ncol=num_cat_altitudeF, dimnames = list(cat_list_altitudeF,cat_list_altitudeF)) #empty matrix
  for (i in 1:num_cat_altitudeF){
    for (j in 1:i){
      if (i!=j){
        # Assemble the hypervolumes two by two
        this_set = hypervolume_set(hv_list_altitudeF@HVList[[i]], hv_list_altitudeF@HVList[[j]], check.memory=FALSE)
        # Calculate Sorensen
        overlap_altitudeF[i,j] = hypervolume_overlap_statistics(this_set)["sorensen"]
        #overlap_altitudeF_CI[i,j] = hypervolume_overlap_confidence(this_set)["sorensen"]
      }  }  } #computing
  
  NOMpng="results/Hypervolume_Traits_by_altitude_Fsex_Sorensen_index.png"
  png(file = NOMpng)
  plot(hclust(as.dist(1-overlap_altitudeF)), main = "(a) Sorensen's similarity index", ylab = "1-Similarity", hang = -1)
  dev.off()
  
  
  
  # Centro?d distances
  distance_altitudeF = matrix(NA, nrow = num_cat_altitudeF,ncol = num_cat_altitudeF, dimnames = list(cat_list_altitudeF,cat_list_altitudeF))
  for(i in 1:num_cat_altitudeF){
    for(j in 1:i){
      if (i!=j){
        distance_altitudeF[i,j]<-hypervolume_distance(hv_list_altitudeF@HVList[[i]], hv_list_altitudeF@HVList[[j]],check.memory = FALSE)
      } } }
  
  NOMpng="results/Hypervolume_Traits_by_altitude_Fsex_Centroid_distance.png"
  png(file = NOMpng)
  plot(hclust(as.dist(distance_altitudeF)), main = "(b) Centroid distance",ylab = "Distance", hang = -1)
  dev.off()
  
  
  # Volumes, contributions, centro?ds observed #### ---------------------------------------------------
  
  # Get volume 
  # "Volume is a measure of hypervolume size and represents the width of the multidimensional,
  # trait space of the altitudes (i.e. the variability of all traits shaping the hypervolume simultaneously)."
  volume_altitudeF<-as.data.frame(get_volume(hv_list_altitudeF))
  
  # Get profile trait contribution and plot them 
  forme_altitudeF=matrix(0,num_cat_altitudeF,num_trait, dimnames = list(cat_list_altitudeF,trait_axes))
  for (i in 1:num_cat_altitudeF){
    forme_altitudeF[i,]<-hypervolume_variable_importance(hv_list_altitudeF@HVList[[i]])
  }
  
  plotforme_altitudeF<-as.data.frame(forme_altitudeF)
  plotforme_altitudeF<- rbind(rep(round(max(plotforme_altitudeF)),num_trait) , rep(0,num_trait) ,plotforme_altitudeF) # Set plot limits
  
  NOMpng="results/Hypervolume_Traits_by_altitude_Fsex_profile.png"
  png(file = NOMpng)
  par(mar=rep(0.8,4))
  par(mfrow=c(2,3)) # To adapt 
  for(i in 1:num_cat_altitudeF){ 
    radarchart(plotforme_altitudeF[c(1,2,i+2),], axistype=2, title=cat_list_altitudeF[i], maxmin = T,
               plwd=3, plty=1, cglcol="grey", cglty=1, cglwd=1, axislabcol="black", vlcex=0.9)
  }
  dev.off()
  
  
  
  # Get centroids 
  centroides_altitudeF<-as.data.frame(get_centroid(hv_list_altitudeF))
  centroides_altitudeF$Altitude <- rownames(centroides_altitudeF)
  
  
  #  Comparing volumes (i.e. volume of each cloud) between altitudes ---------------------------------------------------
  
  # Bootstraping volumes within each altitude to compare them
  
  table(dataF$Altitude) # see number of individuals by altitude
  sim=100 # Number of simulation
  ind=7 # Number of individual = Number of individuals for the case where there are the least

  volumeBS_altitudeF=matrix(0,num_cat_altitudeF,sim) # Row = Altitude Col = resampling
  rownames(volumeBS_altitudeF)=cat_list_altitudeF
  for (i in 1:num_cat_altitudeF){  
    dataF_ce_cat_altitudeF = dataF[dataF$Altitude==cat_list_altitudeF[i],trait_axes] # dataF for this altitude
    for(j in 1:sim){
      
      dataFBS_altitudeF = dataF_ce_cat_altitudeF[sample(nrow(dataF_ce_cat_altitudeF),ind,replace=T),] # sampling within the altitude
      hyperBS_altitudeF = hypervolume_gaussian(dataFBS_altitudeF,verbose=FALSE, kde.bandwidth = estimate_bandwidth(dataFBS_altitudeF))
      volumeBS_altitudeF[i,j]<-get_volume(hyperBS_altitudeF)
    }}
  
  ## Boxplot of resampled (i.e. boostraped) volumes 
  plotvolumeBS_altitudeF<-as.data.frame(volumeBS_altitudeF)
  plotvolumeBS_altitudeF$Altitude <- rownames(plotvolumeBS_altitudeF)
  plotvolumeBS_altitudeF<-melt(plotvolumeBS_altitudeF,  id="Altitude")
  
  NOMpng="results/Hypervolume_Traits_by_altitude_Fsex_volumeBS.png"
  png(file = NOMpng,width = 400, height =300)
  ggplot(plotvolumeBS_altitudeF, aes(Altitude, value, fill=Altitude))+
    geom_boxplot() + scale_fill_viridis_d(begin=0.3,end=0.9)+ theme_bw()
  dev.off()
  

  
  ## Ztest to test significant difference between volumes
  ztest_volume_altitudeF <- matrix(NA,num_cat_altitudeF,num_cat_altitudeF, dimnames = list(cat_list_altitudeF,cat_list_altitudeF))
  num_comp = num_cat_altitudeF*(num_cat_altitudeF-1)/2 # Number of comparisons 
  for(i in 1:num_cat_altitudeF){ 
    for(j in 1:i){
      if (i!=j){
        ztest_volume_altitudeF[i,j]<- ztest_function(plotvolumeBS_altitudeF[plotvolumeBS_altitudeF$Altitude==cat_list_altitudeF[i],"value"],
                                                     plotvolumeBS_altitudeF[plotvolumeBS_altitudeF$Altitude==cat_list_altitudeF[j],"value"])
      } } }
  
  # Plot resampled volumes (i.e. boostraped volumes) : mean, sd, letters for differences
  volumeztest_altitudeF<-as.data.frame(tapply(plotvolumeBS_altitudeF$value, plotvolumeBS_altitudeF$Altitude, mean))
  colnames(volumeztest_altitudeF)="mean"
  volumeztest_altitudeF$sd<-tapply(plotvolumeBS_altitudeF$value, plotvolumeBS_altitudeF$Altitude, sd)
  volumeztest_altitudeF$Altitude<-rownames(volumeztest_altitudeF)
  volumeztest_altitudeF$lowersd <- volumeztest_altitudeF$mean-volumeztest_altitudeF$sd
  volumeztest_altitudeF$uppersd <- volumeztest_altitudeF$mean+volumeztest_altitudeF$sd
  volumeztest_altitudeF$letter <- c("ab","a","ab","ab","b")## Put the letters manually according to the significant differences in ztest_volume_altitudeF
  
  NOMpng="results/Hypervolume_Traits_by_altitude_Fsex_volumeZtest.png"
  png(file = NOMpng,width = 400, height =300)
  ggplot(volumeztest_altitudeF,aes(x= Altitude, y=mean, color=Altitude,label = letter)) +
    geom_point(shape  = 16, size   = 6) + 
    geom_errorbar(aes(ymin  = lowersd, ymax  =  uppersd), width =  0.1,size  =  1.2) +
    theme_bw() + ylab(expression(Volume~~(SD^5))) + ggtitle ("") +
    scale_color_viridis_d(begin=0.3,end=0.9) + 
    geom_text(aes(y = 150), color  = "black") +
    ylim(0,500)+
    theme(text= element_text(size=18), axis.title.x = element_blank(), 
          plot.title = element_text(size=18, face="italic"),legend.position = "none")
  
  dev.off()
  
  
  
  
  
  # Comparing centroids between altitude ####---------------------------------------------------
  
  # Boostraping centroids
  table(dataF$Altitude) 
  sim=100
  indi=7
  
  centroBS_altitudeF=matrix(0,sim,num_trait) # Resampled hypervolume centro?ds Row = simulation, Col = traits
  formeBS_altitudeF=matrix(0,sim,num_trait)  # Resampled hypervolume trait contributions Row = simulation, Col = traits
  for(j in 1:sim){
    
    dataFBS_altitudeF = dataF[sample(nrow(dataF),indi,replace=F),trait_axes] 
    hyperBS_altitudeF = hypervolume_gaussian(dataFBS_altitudeF,verbose=FALSE, kde.bandwidth = estimate_bandwidth(dataFBS_altitudeF))
    centroBS_altitudeF[j,(1:num_trait)]<-get_centroid(hyperBS_altitudeF)
    formeBS_altitudeF[j,(1:num_trait)] <- hypervolume_variable_importance(hyperBS_altitudeF) 
  }
  
  
  #### Comparison of centroids between altitudes and the null model over all altitudes ####
  
  #Create a matrix containing the lower bounds, a matrix containing the upper bounds of the CI and 
  # then compare with the observed values
  Centroides_bornes_inf_altitudeF=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Centroides_bornes_inf_altitudeF[1,i]<-quantile(centroBS_altitudeF[,i], 0.025)
  }
  
  Centroides_bornes_sup_altitudeF=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Centroides_bornes_sup_altitudeF[1,i]<-quantile(centroBS_altitudeF[,i], 0.975)
  }
  
  # Plot the results
  NOMpng="results/Hypervolume_Traits_by_altitude_Fsex_plotgrid_centroides.png"
  png(file = NOMpng)
  plot_grid(
    ggplot(centroides_altitudeF, aes(x = Altitude, y = Mass, color=Altitude))+
      ggtitle("") + ylab("Mass")  + xlab("Altitude") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_altitudeF[1,1], 
                    ymax = Centroides_bornes_sup_altitudeF[1,1]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.3,end=0.9) + theme_bw() +
      theme(axis.title.x = element_blank(), legend.position = "none", plot.title = element_text(size=12, face="bold.italic")),
    
    ggplot(centroides_altitudeF, aes(x = Altitude, y = Body_size_index_corr, color=Altitude))+
      ggtitle("") + ylab("Body_size_index")  + xlab("Altitude") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_altitudeF[1,2], 
                    ymax = Centroides_bornes_sup_altitudeF[1,2]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.3,end=0.9) + theme_bw() +
      theme(axis.title.x = element_blank(), legend.position = "none", plot.title = element_text(size=12, face="bold.italic")),
    
    ggplot(centroides_altitudeF, aes(x = Altitude, y = Proteines.Mass, color=Altitude))+
      ggtitle("")+ ylab("Proteins:Mass") + xlab("Altitude") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_altitudeF[1,3], 
                    ymax = Centroides_bornes_sup_altitudeF[1,3]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.3,end=0.9) + theme_bw() +
      theme(axis.title.x = element_blank(),legend.position = "none"),
    
    ggplot(centroides_altitudeF, aes(x = Altitude, y = Lipids.Mass, color=Altitude))+
      ggtitle("")+ ylab("Lipids:Mass") + xlab("Altitude") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_altitudeF[1,4], 
                    ymax = Centroides_bornes_sup_altitudeF[1,4]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.3,end=0.9)+ theme_bw() +
      theme(axis.title.x = element_blank(),legend.position = "none"),
    
    ggplot(centroides_altitudeF, aes(x = Altitude, y = Sugars.Mass, color=Altitude))+
      ggtitle("")+ ylab("Sugars:Mass") + xlab("Altitude") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_altitudeF[1,5], 
                    ymax = Centroides_bornes_sup_altitudeF[1,5]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.3,end=0.9)+ theme_bw() +
      theme(axis.title.x = element_blank(),legend.position = "none"),
    
    ncol = 3, nrow = 2, byrow = F)
  dev.off()
  
  
  
  ### Comparison of trait contributions between altitudes and the null model over all altitudes #### ---------------------------------------------------------------------
  
  ##Create a matrix containing the lower bounds, a matrix containing the upper bounds of the CI and then compare with the observed values
  
  Forme_bornes_inf_altitudeF=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Forme_bornes_inf_altitudeF[1,i]<-quantile(formeBS_altitudeF[,i], 0.025)
  }
  
  Forme_bornes_sup_altitudeF=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Forme_bornes_sup_altitudeF[1,i]<-quantile(formeBS_altitudeF[,i], 0.975)
  }
  
  Forme_IC_altitudeF=matrix(NA,num_cat_altitudeF,num_trait, dimnames = list(cat_list_altitudeF,trait_axes))
  for (i in 1:num_cat_altitudeF){
    for (j in 1:num_trait){
      Forme_IC_altitudeF[i,j]<- if (forme_altitudeF[i,j]<Forme_bornes_inf_altitudeF[1,j]){"-"} 
      else if (forme_altitudeF[i,j]>Forme_bornes_sup_altitudeF[1,j]){"+"} else {"ns"}
    }}
  
  # There are certainly simpler than all these lines but it works...
  plotformeIC_altitudeF<-as.data.frame(forme_altitudeF)
  plotformeIC_altitudeF$Altitude<-rownames(plotformeIC_altitudeF)
  plotformeIC_altitudeF$type<-rep("obs",num_cat_altitudeF)
  Forme_IC_altitudeF=as.data.frame(Forme_IC_altitudeF)
  Forme_IC_altitudeF$Altitude<-rownames(Forme_IC_altitudeF)
  Forme_bornes_inf_altitudeF=as.data.frame(Forme_bornes_inf_altitudeF)
  Forme_bornes_inf_altitudeF$Altitude<-rownames(Forme_bornes_inf_altitudeF)
  Forme_bornes_inf_altitudeF$type<-"inf"
  Forme_bornes_sup_altitudeF=as.data.frame(Forme_bornes_sup_altitudeF)
  Forme_bornes_sup_altitudeF$Altitude<-rownames(Forme_bornes_sup_altitudeF)
  Forme_bornes_sup_altitudeF$type<-"sup"
  plotformeIC_altitudeF<-bind_rows(plotformeIC_altitudeF,Forme_bornes_inf_altitudeF,Forme_bornes_sup_altitudeF)
  plotformeIC_altitudeF <-rbind(plotformeIC_altitudeF, rep(0,num_trait)) # Limit inf
  plotformeIC_altitudeF <-rbind(plotformeIC_altitudeF, rep(round(max(plotformeIC_altitudeF[,1:num_trait])),num_trait)) 
  
  colnames(plotformeIC_altitudeF)<- c(trait_axes,"Altitude","type")
  
  NOMpng="results/Hypervolume_Traits_by_altitude_Fsex_radarchart.png"
  png(file = NOMpng)
  par(mar=rep(0.8,4))
  par(mfrow=c(1,1)) #Adapt
  coul <- c("grey", "grey", viridis_col) 
  coul_leg <- c("white", "white", viridis_col)
  coul_in <- c("grey","white",rep(alpha("white",0.1),num_cat_altitudeF))
  radarchart(plotformeIC_altitudeF[c(9,8,7,6,1:5),1:num_trait], axistype=6, title="", maxmin = T,
             pcol=coul, plwd=2, plty=1, pfcol = coul_in, cglcol="grey", cglty=2, cglwd=1, axislabcol="black", vlcex=1)
  legend(x="bottomleft", legend = c("","",cat_list_altitudeF), 
         bty = "n", pch=20 , col= coul_leg , text.col = "black", cex=1, pt.cex=3, horiz = F)
  
  dev.off()
  
  save.image("results/insect_hypervolumes_altitude_Fsex.Rdata")
}

wantAnalysisByAltitudeM=F
if(wantAnalysisByAltitudeM==TRUE){
  
  # Check the number of individual per hypervolume (i.e. per altitude)
  table(dataM$Altitude)
  
  
  # Building hypervolumes ---------------------------------------------------
  
  # list the hypervolumes to be calculated
  cat_list_altitudeM = sort(as.character(unique(dataM$Altitude)))
  num_cat_altitudeM = length(cat_list_altitudeM)
  print(paste0("The number of altitude categories is : ", num_cat_altitudeM))
  
  # Calculate hypervolumes (1 per altitude) and join all of them
  hv_list_altitudeM = new("HypervolumeList") # list that will contain all the hypervolumes
  hv_list_altitudeM@HVList = vector(mode="list",length=num_cat_altitudeM)
  for (i in 1:num_cat_altitudeM){ 
    # select dataM
    dataM_ce_cat_altitudeM = dataM[dataM$Altitude==cat_list_altitudeM[i],trait_axes]
    # Calculate hypervolume for each altitude
    hv_list_altitudeM@HVList[[i]] <- hypervolume_gaussian(dataM_ce_cat_altitudeM,name=as.character(cat_list_altitudeM[i]),
                                                          verbose=FALSE,kde.bandwidth = estimate_bandwidth(dataM_ce_cat_altitudeM))
  }
  
  NOMpng="results/Hypervolume_Traits_by_altitude_Msex.png"
  png(file = NOMpng,width=400,height = 400)
  plot(hypervolume_join(hv_list_altitudeM@HVList[[1]], hv_list_altitudeM@HVList[[2]], hv_list_altitudeM[[3]], 
                        hv_list_altitudeM@HVList[[4]], hv_list_altitudeM@HVList[[5]]),
       show.density = T, show.dataM = T, show.contour = F, cex.legend = 1.6, cex.names = 1.6, col= viridis_col, 
       names = trait_axes)
  dev.off()
  
  
  # Assessing similatity between hypervolumes ---------------------------------------------------
  
  # Sorensen similarity index
  overlap_altitudeM = matrix(NA, nrow=num_cat_altitudeM, ncol=num_cat_altitudeM, dimnames = list(cat_list_altitudeM,cat_list_altitudeM)) #empty matrix
  #overlap_altitudeM_CI = matrix(NA, nrow=num_cat_altitudeM, ncol=num_cat_altitudeM, dimnames = list(cat_list_altitudeM,cat_list_altitudeM)) #empty matrix
  for (i in 1:num_cat_altitudeM){
    for (j in 1:i){
      if (i!=j){
        # Assemble the hypervolumes two by two
        this_set = hypervolume_set(hv_list_altitudeM@HVList[[i]], hv_list_altitudeM@HVList[[j]], check.memory=FALSE)
        # Calculate Sorensen
        overlap_altitudeM[i,j] = hypervolume_overlap_statistics(this_set)["sorensen"]
        #overlap_altitudeM_CI[i,j] = hypervolume_overlap_confidence(this_set)["sorensen"]
      }  }  } #computing
  
  NOMpng="results/Hypervolume_Traits_by_altitude_Msex_Sorensen_index.png"
  png(file = NOMpng)
  plot(hclust(as.dist(1-overlap_altitudeM)), main = "(a) Sorensen's similarity index", ylab = "1-Similarity", hang = -1)
  dev.off()
  
  
  
  # Centroid distances
  distance_altitudeM = matrix(NA, nrow = num_cat_altitudeM,ncol = num_cat_altitudeM, dimnames = list(cat_list_altitudeM,cat_list_altitudeM))
  for(i in 1:num_cat_altitudeM){
    for(j in 1:i){
      if (i!=j){
        distance_altitudeM[i,j]<-hypervolume_distance(hv_list_altitudeM@HVList[[i]], hv_list_altitudeM@HVList[[j]],check.memory = FALSE)
      } } }
  
  NOMpng="results/Hypervolume_Traits_by_altitude_Msex_Centroid_distance.png"
  png(file = NOMpng)
  plot(hclust(as.dist(distance_altitudeM)), main = "(b) Centroid distance",ylab = "Distance", hang = -1)
  dev.off()
  
  
  # Volumes, contributions, centro?ds observed #### ---------------------------------------------------
  
  # Get volume 
  volume_altitudeM<-as.data.frame(get_volume(hv_list_altitudeM))
  
  # Get profile trait contribution and plot them 
  forme_altitudeM=matrix(0,num_cat_altitudeM,num_trait, dimnames = list(cat_list_altitudeM,trait_axes))
  for (i in 1:num_cat_altitudeM){
    forme_altitudeM[i,]<-hypervolume_variable_importance(hv_list_altitudeM@HVList[[i]])
  }
  
  plotforme_altitudeM<-as.data.frame(forme_altitudeM)
  plotforme_altitudeM<- rbind(rep(round(max(plotforme_altitudeM)),num_trait) , rep(0,num_trait) ,plotforme_altitudeM) # Set plot limits
  
  NOMpng="results/Hypervolume_Traits_by_altitude_Msex_profile.png"
  png(file = NOMpng)
  par(mar=rep(0.8,4))
  par(mfrow=c(2,3)) # To adapt 
  for(i in 1:num_cat_altitudeM){ 
    radarchart(plotforme_altitudeM[c(1,2,i+2),], axistype=2, title=cat_list_altitudeM[i], maxmin = T,
               plwd=3, plty=1, cglcol="grey", cglty=1, cglwd=1, axislabcol="black", vlcex=0.9)
  }
  dev.off()
  
  
  
  
  # Get centro?ds 
  centroides_altitudeM<-as.data.frame(get_centroid(hv_list_altitudeM))
  centroides_altitudeM$Altitude <- rownames(centroides_altitudeM)
  
  
  #  Comparing volumes (i.e. volume of each cloud) between altitudes ---------------------------------------------------
  
  # Bootstraping volumes within each altitude to compare them
  
  table(dataM$Altitude) # see number of individuals by altitude
  sim=100 # Number of simulation
  ind=10 # Number of individual = Number of individuals for the case where there are the least

  volumeBS_altitudeM=matrix(0,num_cat_altitudeM,sim) # Row = Altitude Col = resampling
  rownames(volumeBS_altitudeM)=cat_list_altitudeM
  for (i in 1:num_cat_altitudeM){  
    dataM_ce_cat_altitudeM = dataM[dataM$Altitude==cat_list_altitudeM[i],trait_axes] # dataM for this altitude
    for(j in 1:sim){
      
      dataMBS_altitudeM = dataM_ce_cat_altitudeM[sample(nrow(dataM_ce_cat_altitudeM),ind,replace=T),] # sampling within the altitude
      hyperBS_altitudeM = hypervolume_gaussian(dataMBS_altitudeM,verbose=FALSE, kde.bandwidth = estimate_bandwidth(dataMBS_altitudeM))
      volumeBS_altitudeM[i,j]<-get_volume(hyperBS_altitudeM)
    }}
  
  ## Boxplot of resampled (i.e. boostraped) volumes 

  plotvolumeBS_altitudeM<-as.data.frame(volumeBS_altitudeM)
  plotvolumeBS_altitudeM$Altitude <- rownames(plotvolumeBS_altitudeM)
  plotvolumeBS_altitudeM<-melt(plotvolumeBS_altitudeM,  id="Altitude")
  
  NOMpng="results/Hypervolume_Traits_by_altitude_Msex_volumeBS.png"
  png(file = NOMpng,width = 400, height =300)
  ggplot(plotvolumeBS_altitudeM, aes(Altitude, value, fill=Altitude))+
    geom_boxplot() + scale_fill_viridis_d(begin=0.3,end=0.9)+ theme_bw()
  dev.off()
  
  
  ## Ztest to test significant difference between volumes
  ztest_volume_altitudeM <- matrix(NA,num_cat_altitudeM,num_cat_altitudeM, dimnames = list(cat_list_altitudeM,cat_list_altitudeM))
  num_comp = num_cat_altitudeM*(num_cat_altitudeM-1)/2 # Number of comparisons 
  for(i in 1:num_cat_altitudeM){ 
    for(j in 1:i){
      if (i!=j){
        ztest_volume_altitudeM[i,j]<- ztest_function(plotvolumeBS_altitudeM[plotvolumeBS_altitudeM$Altitude==cat_list_altitudeM[i],"value"],
                                                     plotvolumeBS_altitudeM[plotvolumeBS_altitudeM$Altitude==cat_list_altitudeM[j],"value"])
      } } }
  
  # Plot resampled volumes (i.e. boostraped volumes) : mean, sd, letters for differences
  volumeztest_altitudeM<-as.data.frame(tapply(plotvolumeBS_altitudeM$value, plotvolumeBS_altitudeM$Altitude, mean))
  colnames(volumeztest_altitudeM)="mean"
  volumeztest_altitudeM$sd<-tapply(plotvolumeBS_altitudeM$value, plotvolumeBS_altitudeM$Altitude, sd)
  volumeztest_altitudeM$Altitude<-rownames(volumeztest_altitudeM)
  volumeztest_altitudeM$lowersd <- volumeztest_altitudeM$mean-volumeztest_altitudeM$sd
  volumeztest_altitudeM$uppersd <- volumeztest_altitudeM$mean+volumeztest_altitudeM$sd
  volumeztest_altitudeM$letter <- c("ab","ab","a","ab","b")## Put the letters manually according to the significant differences in ztest_volume_altitudeM
  
  NOMpng="results/Hypervolume_Traits_by_altitude_Msex_volumeZtest.png"
  png(file = NOMpng,width = 400, height =300)
  ggplot(volumeztest_altitudeM,aes(x= Altitude, y=mean, color=Altitude,label = letter)) +
    geom_point(shape  = 16, size   = 6) + 
    geom_errorbar(aes(ymin  = lowersd, ymax  =  uppersd), width =  0.1,size  =  1.2) +
    theme_bw() + ylab(expression(Volume~~(SD^5))) + ggtitle ("") +
    scale_color_viridis_d(begin=0.3,end=0.9) + 
    geom_text(aes(y = 150), color  = "black") +
    ylim(0,500)+
    theme(text= element_text(size=18), axis.title.x = element_blank(), 
          plot.title = element_text(size=18, face="italic"),legend.position = "none")
  
  dev.off()
  
  
  
  # Comparing centroids between altitude ####---------------------------------------------------
  
  # Boostraping centroids
  table(dataM$Altitude) 
  sim=100
  indi=10
  
  centroBS_altitudeM=matrix(0,sim,num_trait) # Resampled hypervolume centro?ds Row = simulation, Col = traits
  formeBS_altitudeM=matrix(0,sim,num_trait)  # Resampled hypervolume trait contributions Row = simulation, Col = traits
  for(j in 1:sim){
    
    dataMBS_altitudeM = dataM[sample(nrow(dataM),indi,replace=F),trait_axes] 
    hyperBS_altitudeM = hypervolume_gaussian(dataMBS_altitudeM,verbose=FALSE, kde.bandwidth = estimate_bandwidth(dataMBS_altitudeM))
    centroBS_altitudeM[j,(1:num_trait)]<-get_centroid(hyperBS_altitudeM)
    formeBS_altitudeM[j,(1:num_trait)] <- hypervolume_variable_importance(hyperBS_altitudeM) 
  }
  
  
  #### Comparison of centroids between altitudes and the null model over all altitudes ####
  
  #Create a matrix containing the lower bounds, a matrix containing the upper bounds of the CI and 
  # then compare with the observed values
  Centroides_bornes_inf_altitudeM=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Centroides_bornes_inf_altitudeM[1,i]<-quantile(centroBS_altitudeM[,i], 0.025)
  }
  
  Centroides_bornes_sup_altitudeM=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Centroides_bornes_sup_altitudeM[1,i]<-quantile(centroBS_altitudeM[,i], 0.975)
  }
  
  # Plot the results
  NOMpng="results/Hypervolume_Traits_by_altitude_Msex_plotgrid_centroides.png"
  png(file = NOMpng)
  plot_grid(
    ggplot(centroides_altitudeM, aes(x = Altitude, y = Mass, color=Altitude))+
      ggtitle("") + ylab("Mass")  + xlab("Altitude") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_altitudeM[1,1], 
                    ymax = Centroides_bornes_sup_altitudeM[1,1]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.3,end=0.9) + theme_bw() +
      theme(axis.title.x = element_blank(), legend.position = "none", plot.title = element_text(size=12, face="bold.italic")),
    
    ggplot(centroides_altitudeM, aes(x = Altitude, y = Body_size_index_corr, color=Altitude))+
      ggtitle("") + ylab("Body_size_index")  + xlab("Altitude") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_altitudeM[1,2], 
                    ymax = Centroides_bornes_sup_altitudeM[1,2]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.3,end=0.9) + theme_bw() +
      theme(axis.title.x = element_blank(), legend.position = "none", plot.title = element_text(size=12, face="bold.italic")),
    
    ggplot(centroides_altitudeM, aes(x = Altitude, y = Proteines.Mass, color=Altitude))+
      ggtitle("")+ ylab("Proteins:Mass") + xlab("Altitude") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_altitudeM[1,3], 
                    ymax = Centroides_bornes_sup_altitudeM[1,3]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.3,end=0.9) + theme_bw() +
      theme(axis.title.x = element_blank(),legend.position = "none"),
    
    ggplot(centroides_altitudeM, aes(x = Altitude, y = Lipids.Mass, color=Altitude))+
      ggtitle("")+ ylab("Lipids:Mass") + xlab("Altitude") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_altitudeM[1,4], 
                    ymax = Centroides_bornes_sup_altitudeM[1,4]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.3,end=0.9)+ theme_bw() +
      theme(axis.title.x = element_blank(),legend.position = "none"),
    
    ggplot(centroides_altitudeM, aes(x = Altitude, y = Sugars.Mass, color=Altitude))+
      ggtitle("")+ ylab("Sugars:Mass") + xlab("Altitude") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_altitudeM[1,5], 
                    ymax = Centroides_bornes_sup_altitudeM[1,5]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.3,end=0.9)+ theme_bw() +
      theme(axis.title.x = element_blank(),legend.position = "none"),
    
    ncol = 3, nrow = 2, byrow = F)
  dev.off()
  

  
  ### Comparison of trait contributions between altitudes and the null model over all altitudes #### ---------------------------------------------------------------------
  
  ##Create a matrix containing the lower bounds, a matrix containing the upper bounds of the CI and then compare with the observed values
  
  Forme_bornes_inf_altitudeM=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Forme_bornes_inf_altitudeM[1,i]<-quantile(formeBS_altitudeM[,i], 0.025)
  }
  
  Forme_bornes_sup_altitudeM=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Forme_bornes_sup_altitudeM[1,i]<-quantile(formeBS_altitudeM[,i], 0.975)
  }
  
  Forme_IC_altitudeM=matrix(NA,num_cat_altitudeM,num_trait, dimnames = list(cat_list_altitudeM,trait_axes))
  for (i in 1:num_cat_altitudeM){
    for (j in 1:num_trait){
      Forme_IC_altitudeM[i,j]<- if (forme_altitudeM[i,j]<Forme_bornes_inf_altitudeM[1,j]){"-"} 
      else if (forme_altitudeM[i,j]>Forme_bornes_sup_altitudeM[1,j]){"+"} else {"ns"}
    }}
  
  # There are certainly simpler than all these lines but it works...
  plotformeIC_altitudeM<-as.data.frame(forme_altitudeM)
  plotformeIC_altitudeM$Altitude<-rownames(plotformeIC_altitudeM)
  plotformeIC_altitudeM$type<-rep("obs",num_cat_altitudeM)
  Forme_IC_altitudeM=as.data.frame(Forme_IC_altitudeM)
  Forme_IC_altitudeM$Altitude<-rownames(Forme_IC_altitudeM)
  Forme_bornes_inf_altitudeM=as.data.frame(Forme_bornes_inf_altitudeM)
  Forme_bornes_inf_altitudeM$Altitude<-rownames(Forme_bornes_inf_altitudeM)
  Forme_bornes_inf_altitudeM$type<-"inf"
  Forme_bornes_sup_altitudeM=as.data.frame(Forme_bornes_sup_altitudeM)
  Forme_bornes_sup_altitudeM$Altitude<-rownames(Forme_bornes_sup_altitudeM)
  Forme_bornes_sup_altitudeM$type<-"sup"
  plotformeIC_altitudeM<-bind_rows(plotformeIC_altitudeM,Forme_bornes_inf_altitudeM,Forme_bornes_sup_altitudeM)
  plotformeIC_altitudeM <-rbind(plotformeIC_altitudeM, rep(0,num_trait)) # Limit inf
  plotformeIC_altitudeM <-rbind(plotformeIC_altitudeM, rep(round(max(plotformeIC_altitudeM[,1:num_trait])),num_trait)) 
  
  colnames(plotformeIC_altitudeM)<- c(trait_axes,"Altitude","type")
  
  NOMpng="results/Hypervolume_Traits_by_altitude_Msex_radarchart.png"
  png(file = NOMpng)
  par(mar=rep(0.8,4))
  par(mfrow=c(1,1)) #Adapt
  coul <- c("grey", "grey", viridis_col) 
  coul_leg <- c("white", "white", viridis_col)
  coul_in <- c("grey","white",rep(alpha("white",0.1),num_cat_altitudeM))
  radarchart(plotformeIC_altitudeM[c(9,8,7,6,1:5),1:num_trait], axistype=6, title="", maxmin = T,
             pcol=coul, plwd=2, plty=1, pfcol = coul_in, cglcol="grey", cglty=2, cglwd=1, axislabcol="black", vlcex=1)
  legend(x="bottomleft", legend = c("","",cat_list_altitudeM), 
         bty = "n", pch=20 , col= coul_leg , text.col = "black", cex=1, pt.cex=3, horiz = F)
  
  dev.off()
  
   save.image("results/insect_hypervolumes_altitude_Msex.Rdata"))
}


#
OverlapFile="results/overlap_altitudeL.csv"
write.table(overlap_altitudeL,OverlapFile,sep=";",dec = ".")

OverlapFile="results/overlap_altitudeD.csv"
write.table(overlap_altitudeD,OverlapFile,sep=";",dec = ".")

OverlapFile="results/overlap_altitudeF.csv"
write.table(overlap_altitudeF,OverlapFile,sep=";",dec = ".")

OverlapFile="results/overlap_altitudeM.csv"
write.table(overlap_altitudeM,OverlapFile,sep=";",dec = ".")

distanceFile="results/distance_altitudeL.csv"
write.table(distance_altitudeL,distanceFile,sep=";",dec = ".")

distanceFile="results/distance_altitudeD.csv"
write.table(distance_altitudeD,distanceFile,sep=";",dec = ".")

distanceFile="results/distance_altitudeF.csv"
write.table(distance_altitudeF,distanceFile,sep=";",dec = ".")

distanceFile="results/distance_altitudeM.csv"
write.table(distance_altitudeM,distanceFile,sep=";",dec = ".")




wantAnalysisBySexL=T
if(wantAnalysisBySexL==TRUE){
  
  # Check the number of individual per hypervolume (i.e. per morphotype)
  table(dataL$Sex)
  
  
  # Building hypervolumes ---------------------------------------------------
  
  # list the hypervolumes to be calculated
  cat_list_sexL = sort(as.character(unique(dataL$Sex)))
  num_cat_sexL = length(cat_list_sexL)
  print(paste0("The number of altitude categories is : ", num_cat_sexL))
  
  # Calculate hypervolumes (1 per morphotype) and join all of them
  hv_list_sexL = new("HypervolumeList") # list that will contain all the hypervolumes
  hv_list_sexL@HVList = vector(mode="list",length=num_cat_sexL)
  for (i in 1:num_cat_sexL){ 
    # select data
    data_ce_cat_sexL = dataL[dataL$Sex==cat_list_sexL[i],trait_axes]
    # Calculate hypervolume for each altitude
    hv_list_sexL@HVList[[i]] <- hypervolume_gaussian(data_ce_cat_sexL,name=as.character(cat_list_sexL[i]),
                                                     verbose=FALSE,kde.bandwidth = estimate_bandwidth(data_ce_cat_sexL))
  }
  
  NOMpng="results/Hypervolume_Traits_by_sex_Lmorpho.png"
  png(file = NOMpng,width=400,height=400)
  plot(hypervolume_join(hv_list_sexL@HVList[[1]], hv_list_sexL@HVList[[2]]),
       show.density = T, show.data = T, show.contour = F, cex.legend = 1.6, cex.names = 1.6, col= viridis_col2, 
       names = trait_axes)
  dev.off()
  

  
  # Assessing similatity between hypervolumes ---------------------------------------------------
  
  # Sorensen similarity index
  overlap_sexL = matrix(NA, nrow=num_cat_sexL, ncol=num_cat_sexL, dimnames = list(cat_list_sexL,cat_list_sexL)) #empty matrix
  #overlap_sexL_CI = matrix(NA, nrow=num_cat_sexL, ncol=num_cat_sexL, dimnames = list(cat_list_sexL,cat_list_sexL)) #empty matrix
  for (i in 1:num_cat_sexL){
    for (j in 1:i){
      if (i!=j){
        # Assemble the hypervolumes two by two
        this_set = hypervolume_set(hv_list_sexL@HVList[[i]], hv_list_sexL@HVList[[j]], check.memory=FALSE)
        # Calculate Sorensen
        overlap_sexL[i,j] = hypervolume_overlap_statistics(this_set)["sorensen"]
        #overlap_sexL_CI[i,j] = hypervolume_overlap_confidence(this_set)["sorensen"]
      }  }  } #computing
  
  
  # Centro?d distances
  distance_sexL = matrix(NA, nrow = num_cat_sexL,ncol = num_cat_sexL, dimnames = list(cat_list_sexL,cat_list_sexL))
  for(i in 1:num_cat_sexL){
    for(j in 1:i){
      if (i!=j){
        distance_sexL[i,j]<-hypervolume_distance(hv_list_sexL@HVList[[i]], hv_list_sexL@HVList[[j]],check.memory = FALSE)
      } } }
  
  
  
  # Volumes, contributions, centro?ds observed #### ---------------------------------------------------
  
  # Get volume 
  volume_sexL<-as.data.frame(get_volume(hv_list_sexL))
  
  # Get profile trait contribution and plot them 
  forme_sexL=matrix(0,num_cat_sexL,num_trait, dimnames = list(cat_list_sexL,trait_axes))
  for (i in 1:num_cat_sexL){
    forme_sexL[i,]<-hypervolume_variable_importance(hv_list_sexL@HVList[[i]])
  }
  
  plotforme_sexL<-as.data.frame(forme_sexL)
  plotforme_sexL<- rbind(rep(round(max(plotforme_sexL)),num_trait) , rep(0,num_trait) ,plotforme_sexL) # Set plot limits
  
  NOMpng="results/Hypervolume_Traits_by_sex Lmorpho_profile.png"
  png(file = NOMpng)
  par(mar=rep(0.8,4))
  par(mfrow=c(1,2)) # To adapt 
  for(i in 1:num_cat_sexL){ 
    radarchart(plotforme_sexL[c(1,2,i+2),], axistype=2, title=cat_list_sexL[i], maxmin = T,
               plwd=3, plty=1, cglcol="grey", cglty=1, cglwd=1, axislabcol="black", vlcex=0.9)
  }
  dev.off()
  

  
  # Get centro?ds 
  centroides_sexL<-as.data.frame(get_centroid(hv_list_sexL))
  centroides_sexL$Sex <- rownames(centroides_sexL)
  
  
  #  Comparing volumes (i.e. volume of each cloud) between morphotypes ---------------------------------------------------
  
  # Bootstraping volumes within each altitude to compare them

  table(dataL$Sex) # see number of individuals by altitude
  sim=100 # Number of simulation
  ind=29 # Number of individual = Number of individuals for the case where there are the least

  volumeBS_sexL=matrix(0,num_cat_sexL,sim) # Row = Altitude Col = resampling
  rownames(volumeBS_sexL)=cat_list_sexL
  for (i in 1:num_cat_sexL){  
    data_ce_cat_sexL = dataL[dataL$Sex==cat_list_sexL[i],trait_axes] # Data for this altitude
    for(j in 1:sim){
      
      dataBS_sexL = data_ce_cat_sexL[sample(nrow(data_ce_cat_sexL),ind,replace=T),] # sampling within the altitude
      hyperBS_sexL = hypervolume_gaussian(dataBS_sexL,verbose=FALSE, kde.bandwidth = estimate_bandwidth(dataBS_sexL))
      volumeBS_sexL[i,j]<-get_volume(hyperBS_sexL)
    }}
  
  ## Boxplot of resampled (i.e. boostraped) volumes 

  plotvolumeBS_sexL<-as.data.frame(volumeBS_sexL)
  plotvolumeBS_sexL$Sex <- rownames(plotvolumeBS_sexL)
  plotvolumeBS_sexL<-melt(plotvolumeBS_sexL,  id="Sex")
  
  NOMpng="results/Hypervolume_Traits_by_sex_Lmorpho_volumeBS.png"
  png(file = NOMpng,width = 400, height =300)
  ggplot(plotvolumeBS_sexL, aes(Sex, value, fill=Sex))+
    geom_boxplot() + scale_fill_viridis_d(begin=0.3,end=0.9)+ theme_bw()
  dev.off()
  
  
  ## Ztest to test significant difference between volumes
  ztest_volume_sexL <- matrix(NA,num_cat_sexL,num_cat_sexL, dimnames = list(cat_list_sexL,cat_list_sexL))
  num_comp = num_cat_sexL*(num_cat_sexL-1)/2 # Number of comparisons 
  for(i in 1:num_cat_sexL){ 
    for(j in 1:i){
      if (i!=j){
        ztest_volume_sexL[i,j]<- ztest_function(plotvolumeBS_sexL[plotvolumeBS_sexL$Sex==cat_list_sexL[i],"value"],
                                                plotvolumeBS_sexL[plotvolumeBS_sexL$Sex==cat_list_sexL[j],"value"])
      } } }
  
  # Plot resampled volumes (i.e. boostraped volumes) : mean, sd, letters for differences
  volumeztest_sexL<-as.data.frame(tapply(plotvolumeBS_sexL$value, plotvolumeBS_sexL$Sex, mean))
  colnames(volumeztest_sexL)="mean"
  volumeztest_sexL$sd<-tapply(plotvolumeBS_sexL$value, plotvolumeBS_sexL$Sex, sd)
  volumeztest_sexL$Sex<-rownames(volumeztest_sexL)
  volumeztest_sexL$lowersd <- volumeztest_sexL$mean-volumeztest_sexL$sd
  volumeztest_sexL$uppersd <- volumeztest_sexL$mean+volumeztest_sexL$sd
  volumeztest_sexL$letter <- c("a","b")## Put the letters manually according to the significant differences in ztest_volume_sexL
  
  NOMpng="results/Hypervolume_Traits_by_sex_Lmorpho_volumeZtest.png"
  png(file = NOMpng,width = 400, height =300)
  ggplot(volumeztest_sexL,aes(x=Sex, y=mean, color=Sex,label = letter)) +
    geom_point(shape  = 16, size   = 6) + 
    geom_errorbar(aes(ymin  = lowersd, ymax  =  uppersd), width =  0.1,size  =  1.2) +
    theme_bw() + ylab(expression(Volume~~(SD^5))) + ggtitle ("") +
    scale_color_viridis_d(begin=0.3,end=0.9) + 
    geom_text(aes(y = 150), color  = "black") +
    theme(text= element_text(size=18), axis.title.x = element_blank(), 
          plot.title = element_text(size=18, face="italic"),legend.position = "none")
  
  dev.off()
  
  
 
  
  
  # Comparing centroids between altitude ####---------------------------------------------------
  
  # Boostraping centroids
  table(dataL$Sex) 
  sim=100
  indi=29
  
  centroBS_sexL=matrix(0,sim,num_trait) # Resampled hypervolume centro?ds Row = simulation, Col = traits
  formeBS_sexL=matrix(0,sim,num_trait)  # Resampled hypervolume trait contributions Row = simulation, Col = traits
  for(j in 1:sim){
    
    dataBS_sexL = dataL[sample(nrow(dataL),indi,replace=F),trait_axes] 
    hyperBS_sexL = hypervolume_gaussian(dataBS_sexL,verbose=FALSE, kde.bandwidth = estimate_bandwidth(dataBS_sexL))
    centroBS_sexL[j,(1:num_trait)]<-get_centroid(hyperBS_sexL)
    formeBS_sexL[j,(1:num_trait)] <- hypervolume_variable_importance(hyperBS_sexL) 
  }
  
  #### Comparison of centroids between morphotypes and the null model over all morphotypes ####
  
  #Create a matrix containing the lower bounds, a matrix containing the upper bounds of the CI and 
  # then compare with the observed values
  Centroides_bornes_inf_sexL=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Centroides_bornes_inf_sexL[1,i]<-quantile(centroBS_sexL[,i], 0.025)
  }
  
  Centroides_bornes_sup_sexL=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Centroides_bornes_sup_sexL[1,i]<-quantile(centroBS_sexL[,i], 0.975)
  }
  
  
  # Plot the results
  library(cowplot)
  NOMpng="results/Hypervolume_Traits_by_sex_Lmorpho_plotgrid_centroides.png"
  png(file = NOMpng)
  plot_grid(
    ggplot(centroides_sexL, aes(x = Sex, y = Mass, color=Sex))+
      ggtitle("") + ylab("Mass")  + xlab("Sex") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_sexL[1,1], 
                    ymax = Centroides_bornes_sup_sexL[1,1]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.3,end=0.9) + theme_bw() +
      theme(axis.title.x = element_blank(), legend.position = "none", plot.title = element_text(size=12, face="bold.italic")),
    
    ggplot(centroides_sexL, aes(x = Sex, y = Body_size_index_corr, color=Sex))+
      ggtitle("") + ylab("Body_size_index")  + xlab("Sex") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_sexL[1,2], 
                    ymax = Centroides_bornes_sup_sexL[1,2]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.3,end=0.9) + theme_bw() +
      theme(axis.title.x = element_blank(), legend.position = "none", plot.title = element_text(size=12, face="bold.italic")),
    
    ggplot(centroides_sexL, aes(x = Sex, y = Proteines.Mass, color=Sex))+
      ggtitle("")+ ylab("Proteins:Mass") + xlab("Sex") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_sexL[1,3], 
                    ymax = Centroides_bornes_sup_sexL[1,3]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.3,end=0.9) + theme_bw() +
      theme(axis.title.x = element_blank(),legend.position = "none"),
    
    ggplot(centroides_sexL, aes(x = Sex, y = Lipids.Mass, color=Sex))+
      ggtitle("")+ ylab("Lipids:Mass") + xlab("Sex") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_sexL[1,4], 
                    ymax = Centroides_bornes_sup_sexL[1,4]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.3,end=0.9)+ theme_bw() +
      theme(axis.title.x = element_blank(),legend.position = "none"),
    
    ggplot(centroides_sexL, aes(x = Sex, y = Sugars.Mass, color=Sex))+
      ggtitle("")+ ylab("Sugars:Mass") + xlab("Sex") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_sexL[1,5], 
                    ymax = Centroides_bornes_sup_sexL[1,5]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.3,end=0.9)+ theme_bw() +
      theme(axis.title.x = element_blank(),legend.position = "none"),
    
    ncol = 3, nrow = 2, byrow = F)
  dev.off()
  
 
  
  ### Comparison of trait contributions between morphotypes and the null model over all morphotypes #### ---------------------------------------------------------------------
  
  ##Create a matrix containing the lower bounds, a matrix containing the upper bounds of the CI and then compare with the observed values
  
  Forme_bornes_inf_sexL=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Forme_bornes_inf_sexL[1,i]<-quantile(formeBS_sexL[,i], 0.025)
  }
  
  Forme_bornes_sup_sexL=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Forme_bornes_sup_sexL[1,i]<-quantile(formeBS_sexL[,i], 0.975)
  }
  
  Forme_IC_sexL=matrix(NA,num_cat_sexL,num_trait, dimnames = list(cat_list_sexL,trait_axes))
  for (i in 1:num_cat_sexL){
    for (j in 1:num_trait){
      Forme_IC_sexL[i,j]<- if (forme_sexL[i,j]<Forme_bornes_inf_sexL[1,j]){"-"} 
      else if (forme_sexL[i,j]>Forme_bornes_sup_sexL[1,j]){"+"} else {"ns"}
    }}
  
  # There are certainly simpler than all these lines but it works...
  plotformeIC_sexL<-as.data.frame(forme_sexL)
  plotformeIC_sexL$Sex<-rownames(plotformeIC_sexL)
  plotformeIC_sexL$type<-rep("obs",num_cat_sexL)
  Forme_IC_sexL=as.data.frame(Forme_IC_sexL)
  Forme_IC_sexL$Sex<-rownames(Forme_IC_sexL)
  Forme_bornes_inf_sexL=as.data.frame(Forme_bornes_inf_sexL)
  Forme_bornes_inf_sexL$Sex<-rownames(Forme_bornes_inf_sexL)
  Forme_bornes_inf_sexL$type<-"inf"
  Forme_bornes_sup_sexL=as.data.frame(Forme_bornes_sup_sexL)
  Forme_bornes_sup_sexL$Sex<-rownames(Forme_bornes_sup_sexL)
  Forme_bornes_sup_sexL$type<-"sup"
  plotformeIC_sexL<-bind_rows(plotformeIC_sexL,Forme_bornes_inf_sexL,Forme_bornes_sup_sexL)
  plotformeIC_sexL <-rbind(plotformeIC_sexL, rep(0,num_trait)) # Limit inf
  plotformeIC_sexL <-rbind(plotformeIC_sexL, rep(round(max(plotformeIC_sexL[,1:num_trait])),num_trait)) 
  
  colnames(plotformeIC_sexL)<- c(trait_axes,"Sex","type")
  
  NOMpng="results/Hypervolume_Traits_by_sex_Lmorpho_radarchart.png"
  png(file = NOMpng)
  par(mar=rep(0.8,4))
  par(mfrow=c(1,1)) #Adapt
  coul <- c("grey", "grey", viridis_col2) 
  coul_leg <- c("white", "white", viridis_col2)
  coul_in <- c("grey","white",rep(alpha("white",0.1),num_cat_sexL))
  radarchart(plotformeIC_sexL[c(6,5,4,3, 1:2),1:num_trait], axistype=6, title="", maxmin = T,
             pcol=coul, plwd=2, plty=1, pfcol = coul_in, cglcol="grey", cglty=2, cglwd=1, axislabcol="black", vlcex=1)
  legend(x="bottomleft", legend = c("","",cat_list_sexL), 
         bty = "n", pch=20 , col= coul_leg , text.col = "black", cex=1, pt.cex=3, horiz = F)
  
  dev.off()
  
 
  save.image("results/insect_hypervolumes_sex_Lmorpho.RData")
  
}

wantAnalysisBysexD=T
if(wantAnalysisBysexD==TRUE){
  
  # Check the number of individual per hypervolume (i.e. per morphotype)
  table(dataD$Sex)
  
  
  # Building hypervolumes ---------------------------------------------------
  
  # list the hypervolumes to be calculated
  cat_list_sexD = sort(as.character(unique(dataD$Sex)))
  num_cat_sexD = length(cat_list_sexD)
  print(paste0("The number of altitude categories is : ", num_cat_sexD))
  
  # Calculate hypervolumes (1 per morphotype) and join all of them
  hv_list_sexD = new("HypervolumeList") # list that will contain all the hypervolumes
  hv_list_sexD@HVList = vector(mode="list",length=num_cat_sexD)
  for (i in 1:num_cat_sexD){ 
    # select data
    data_ce_cat_sexD = dataD[dataD$Sex==cat_list_sexD[i],trait_axes]
    # Calculate hypervolume for each altitude
    hv_list_sexD@HVList[[i]] <- hypervolume_gaussian(data_ce_cat_sexD,name=as.character(cat_list_sexD[i]),
                                                     verbose=FALSE,kde.bandwidth = estimate_bandwidth(data_ce_cat_sexD))
  }
  
  NOMpng="results/Hypervolume_Traits_by_sex_Dmorpho.png"
  png(file = NOMpng,width=400,height=400)
  plot(hypervolume_join(hv_list_sexD@HVList[[1]], hv_list_sexD@HVList[[2]]),
       show.density = T, show.data = T, show.contour = F, cex.legend = 1.6, cex.names = 1.6, col= viridis_col2, 
       names = trait_axes)
  dev.off()
  

  # Assessing similatity between hypervolumes ---------------------------------------------------
  
  # Sorensen similarity index
  overlap_sexD = matrix(NA, nrow=num_cat_sexD, ncol=num_cat_sexD, dimnames = list(cat_list_sexD,cat_list_sexD)) #empty matrix
  #overlap_sexD_CI = matrix(NA, nrow=num_cat_sexD, ncol=num_cat_sexD, dimnames = list(cat_list_sexD,cat_list_sexD)) #empty matrix
  for (i in 1:num_cat_sexD){
    for (j in 1:i){
      if (i!=j){
        # Assemble the hypervolumes two by two
        this_set = hypervolume_set(hv_list_sexD@HVList[[i]], hv_list_sexD@HVList[[j]], check.memory=FALSE)
        # Calculate Sorensen
        overlap_sexD[i,j] = hypervolume_overlap_statistics(this_set)["sorensen"]
        #overlap_sexD_CI[i,j] = hypervolume_overlap_confidence(this_set)["sorensen"]
      }  }  } #computing
  
  
  # Centro?d distances
  distance_sexD = matrix(NA, nrow = num_cat_sexD,ncol = num_cat_sexD, dimnames = list(cat_list_sexD,cat_list_sexD))
  for(i in 1:num_cat_sexD){
    for(j in 1:i){
      if (i!=j){
        distance_sexD[i,j]<-hypervolume_distance(hv_list_sexD@HVList[[i]], hv_list_sexD@HVList[[j]],check.memory = FALSE)
      } } }
  
  
  
  # Volumes, contributions, centro?ds observed #### ---------------------------------------------------
  
  # Get volume 
  volume_sexD<-as.data.frame(get_volume(hv_list_sexD))
  
  # Get profile trait contribution and plot them 
  forme_sexD=matrix(0,num_cat_sexD,num_trait, dimnames = list(cat_list_sexD,trait_axes))
  for (i in 1:num_cat_sexD){
    forme_sexD[i,]<-hypervolume_variable_importance(hv_list_sexD@HVList[[i]])
  }
  
  plotforme_sexD<-as.data.frame(forme_sexD)
  plotforme_sexD<- rbind(rep(round(max(plotforme_sexD)),num_trait) , rep(0,num_trait) ,plotforme_sexD) # Set plot limits
  
  NOMpng="results/Hypervolume_Traits_by_sex_Dmorpho_profile.png"
  png(file = NOMpng)
  par(mar=rep(0.8,4))
  par(mfrow=c(1,2)) # To adapt 
  for(i in 1:num_cat_sexD){ 
    radarchart(plotforme_sexD[c(1,2,i+2),], axistype=2, title=cat_list_sexD[i], maxmin = T,
               plwd=3, plty=1, cglcol="grey", cglty=1, cglwd=1, axislabcol="black", vlcex=0.9)
  }
  dev.off()
  

  
  # Get centro?ds 
  centroides_sexD<-as.data.frame(get_centroid(hv_list_sexD))
  centroides_sexD$Sex <- rownames(centroides_sexD)
  
  
  #  Comparing volumes (i.e. volume of each cloud) between morphotypes ---------------------------------------------------
  
  # Bootstraping volumes within each altitude to compare them

  table(dataD$Sex) # see number of individuals by altitude
  sim=100 # Number of simulation
  ind=26 # Number of individual = Number of individuals for the case where there are the least

  volumeBS_sexD=matrix(0,num_cat_sexD,sim) # Row = Altitude Col = resampling
  rownames(volumeBS_sexD)=cat_list_sexD
  for (i in 1:num_cat_sexD){  
    data_ce_cat_sexD = dataD[dataD$Sex==cat_list_sexD[i],trait_axes] # Data for this altitude
    for(j in 1:sim){
      
      dataBS_sexD = data_ce_cat_sexD[sample(nrow(data_ce_cat_sexD),ind,replace=T),] # sampling within the altitude
      hyperBS_sexD = hypervolume_gaussian(dataBS_sexD,verbose=FALSE, kde.bandwidth = estimate_bandwidth(dataBS_sexD))
      volumeBS_sexD[i,j]<-get_volume(hyperBS_sexD)
    }}
  
  ## Boxplot of resampled (i.e. boostraped) volumes 

  plotvolumeBS_sexD<-as.data.frame(volumeBS_sexD)
  plotvolumeBS_sexD$Sex <- rownames(plotvolumeBS_sexD)
  plotvolumeBS_sexD<-melt(plotvolumeBS_sexD,  id="Sex")
  
  NOMpng="results/Hypervolume_Traits_by_sex_Dmorpho_volumeBS.png"
  png(file = NOMpng,width = 400, height =300)
  ggplot(plotvolumeBS_sexD, aes(Sex, value, fill=Sex))+
    geom_boxplot() + scale_fill_viridis_d(begin=0.3,end=0.9)+ theme_bw()
  dev.off()
  

  ## Ztest to test significant difference between volumes
  ztest_volume_sexD <- matrix(NA,num_cat_sexD,num_cat_sexD, dimnames = list(cat_list_sexD,cat_list_sexD))
  num_comp = num_cat_sexD*(num_cat_sexD-1)/2 # Number of comparisons 
  for(i in 1:num_cat_sexD){ 
    for(j in 1:i){
      if (i!=j){
        ztest_volume_sexD[i,j]<- ztest_function(plotvolumeBS_sexD[plotvolumeBS_sexD$Sex==cat_list_sexD[i],"value"],
                                                plotvolumeBS_sexD[plotvolumeBS_sexD$Sex==cat_list_sexD[j],"value"])
      } } }
  
  # Plot resampled volumes (i.e. boostraped volumes) : mean, sd, letters for differences
  volumeztest_sexD<-as.data.frame(tapply(plotvolumeBS_sexD$value, plotvolumeBS_sexD$Sex, mean))
  colnames(volumeztest_sexD)="mean"
  volumeztest_sexD$sd<-tapply(plotvolumeBS_sexD$value, plotvolumeBS_sexD$Sex, sd)
  volumeztest_sexD$Sex<-rownames(volumeztest_sexD)
  volumeztest_sexD$lowersd <- volumeztest_sexD$mean-volumeztest_sexD$sd
  volumeztest_sexD$uppersd <- volumeztest_sexD$mean+volumeztest_sexD$sd
  volumeztest_sexD$letter <- c("a","b")## Put the letters manually according to the significant differences in ztest_volume_sexD
  
  NOMpng="results/Hypervolume_Traits_by_sex_Dmorpho_volumeZtest.png"
  png(file = NOMpng,width = 400, height =300)
  ggplot(volumeztest_sexD,aes(x=Sex, y=mean, color=Sex,label = letter)) +
    geom_point(shape  = 16, size   = 6) + 
    geom_errorbar(aes(ymin  = lowersd, ymax  =  uppersd), width =  0.1,size  =  1.2) +
    theme_bw() + ylab(expression(Volume~~(SD^5))) + ggtitle ("") +
    scale_color_viridis_d(begin=0.3,end=0.9) + 
    geom_text(aes(y = 150), color  = "black") +
    theme(text= element_text(size=18), axis.title.x = element_blank(), 
          plot.title = element_text(size=18, face="italic"),legend.position = "none")
  
  dev.off()
  
  
 
  
  # Comparing centroids between altitude ####---------------------------------------------------
  
  # Boostraping centroids
  table(dataD$Sex) 
  sim=100
  indi=26
  
  centroBS_sexD=matrix(0,sim,num_trait) # Resampled hypervolume centro?ds Row = simulation, Col = traits
  formeBS_sexD=matrix(0,sim,num_trait)  # Resampled hypervolume trait contributions Row = simulation, Col = traits
  for(j in 1:sim){
    
    dataBS_sexD = dataD[sample(nrow(dataD),indi,replace=F),trait_axes] 
    hyperBS_sexD = hypervolume_gaussian(dataBS_sexD,verbose=FALSE, kde.bandwidth = estimate_bandwidth(dataBS_sexD))
    centroBS_sexD[j,(1:num_trait)]<-get_centroid(hyperBS_sexD)
    formeBS_sexD[j,(1:num_trait)] <- hypervolume_variable_importance(hyperBS_sexD) 
  }
  
  #### Comparison of centroids between morphotypes and the null model over all morphotypes ####
  
  #Create a matrix containing the lower bounds, a matrix containing the upper bounds of the CI and 
  # then compare with the observed values
  Centroides_bornes_inf_sexD=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Centroides_bornes_inf_sexD[1,i]<-quantile(centroBS_sexD[,i], 0.025)
  }
  
  Centroides_bornes_sup_sexD=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Centroides_bornes_sup_sexD[1,i]<-quantile(centroBS_sexD[,i], 0.975)
  }
  
  
  # Plot the results
  library(cowplot)
  NOMpng="results/Hypervolume_Traits_by_sex_Dmorpho_plotgrid_centroides.png"
  png(file = NOMpng)
  plot_grid(
    ggplot(centroides_sexD, aes(x = Sex, y = Mass, color=Sex))+
      ggtitle("") + ylab("Mass")  + xlab("Sex") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_sexD[1,1], 
                    ymax = Centroides_bornes_sup_sexD[1,1]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.3,end=0.9) + theme_bw() +
      theme(axis.title.x = element_blank(), legend.position = "none", plot.title = element_text(size=12, face="bold.italic")),
    
    ggplot(centroides_sexD, aes(x = Sex, y = Body_size_index_corr, color=Sex))+
      ggtitle("") + ylab("Body_size_index")  + xlab("Sex") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_sexD[1,2], 
                    ymax = Centroides_bornes_sup_sexD[1,2]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.3,end=0.9) + theme_bw() +
      theme(axis.title.x = element_blank(), legend.position = "none", plot.title = element_text(size=12, face="bold.italic")),
    
    ggplot(centroides_sexD, aes(x = Sex, y = Proteines.Mass, color=Sex))+
      ggtitle("")+ ylab("Proteins:Mass") + xlab("Sex") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_sexD[1,3], 
                    ymax = Centroides_bornes_sup_sexD[1,3]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.3,end=0.9) + theme_bw() +
      theme(axis.title.x = element_blank(),legend.position = "none"),
    
    ggplot(centroides_sexD, aes(x = Sex, y = Lipids.Mass, color=Sex))+
      ggtitle("")+ ylab("Lipids:Mass") + xlab("Sex") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_sexD[1,4], 
                    ymax = Centroides_bornes_sup_sexD[1,4]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.3,end=0.9)+ theme_bw() +
      theme(axis.title.x = element_blank(),legend.position = "none"),
    
    ggplot(centroides_sexD, aes(x = Sex, y = Sugars.Mass, color=Sex))+
      ggtitle("")+ ylab("Sugars:Mass") + xlab("Sex") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_sexD[1,5], 
                    ymax = Centroides_bornes_sup_sexD[1,5]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.3,end=0.9)+ theme_bw() +
      theme(axis.title.x = element_blank(),legend.position = "none"),
    
    ncol = 3, nrow = 2, byrow = F)
  dev.off()
  
  
  ### Comparison of trait contributions between morphotypes and the null model over all morphotypes #### ---------------------------------------------------------------------
  
  ##Create a matrix containing the lower bounds, a matrix containing the upper bounds of the CI and then compare with the observed values
  
  Forme_bornes_inf_sexD=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Forme_bornes_inf_sexD[1,i]<-quantile(formeBS_sexD[,i], 0.025)
  }
  
  Forme_bornes_sup_sexD=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Forme_bornes_sup_sexD[1,i]<-quantile(formeBS_sexD[,i], 0.975)
  }
  
  Forme_IC_sexD=matrix(NA,num_cat_sexD,num_trait, dimnames = list(cat_list_sexD,trait_axes))
  for (i in 1:num_cat_sexD){
    for (j in 1:num_trait){
      Forme_IC_sexD[i,j]<- if (forme_sexD[i,j]<Forme_bornes_inf_sexD[1,j]){"-"} 
      else if (forme_sexD[i,j]>Forme_bornes_sup_sexD[1,j]){"+"} else {"ns"}
    }}
  
  # There are certainly simpler than all these lines but it works...
  plotformeIC_sexD<-as.data.frame(forme_sexD)
  plotformeIC_sexD$Sex<-rownames(plotformeIC_sexD)
  plotformeIC_sexD$type<-rep("obs",num_cat_sexD)
  Forme_IC_sexD=as.data.frame(Forme_IC_sexD)
  Forme_IC_sexD$Sex<-rownames(Forme_IC_sexD)
  Forme_bornes_inf_sexD=as.data.frame(Forme_bornes_inf_sexD)
  Forme_bornes_inf_sexD$Sex<-rownames(Forme_bornes_inf_sexD)
  Forme_bornes_inf_sexD$type<-"inf"
  Forme_bornes_sup_sexD=as.data.frame(Forme_bornes_sup_sexD)
  Forme_bornes_sup_sexD$Sex<-rownames(Forme_bornes_sup_sexD)
  Forme_bornes_sup_sexD$type<-"sup"
  plotformeIC_sexD<-bind_rows(plotformeIC_sexD,Forme_bornes_inf_sexD,Forme_bornes_sup_sexD)
  plotformeIC_sexD <-rbind(plotformeIC_sexD, rep(0,num_trait)) # Limit inf
  plotformeIC_sexD <-rbind(plotformeIC_sexD, rep(round(max(plotformeIC_sexD[,1:num_trait])),num_trait)) 
  
  colnames(plotformeIC_sexD)<- c(trait_axes,"Sex","type")
  
  NOMpng="results/Hypervolume_Traits_by_sex_Dmorpho_radarchart.png"
  png(file = NOMpng)
  par(mar=rep(0.8,4))
  par(mfrow=c(1,1)) #Adapt
  coul <- c("grey", "grey", viridis_col2) 
  coul_leg <- c("white", "white", viridis_col2)
  coul_in <- c("grey","white",rep(alpha("white",0.1),num_cat_sexD))
  radarchart(plotformeIC_sexD[c(6,5,4,3, 1:2),1:num_trait], axistype=6, title="", maxmin = T,
             pcol=coul, plwd=2, plty=1, pfcol = coul_in, cglcol="grey", cglty=2, cglwd=1, axislabcol="black", vlcex=1)
  legend(x="bottomleft", legend = c("","",cat_list_sexD), 
         bty = "n", pch=20 , col= coul_leg , text.col = "black", cex=1, pt.cex=3, horiz = F)
  
  dev.off()
  

  save.image("results/insect_hypervolumes_sex_Dmorpho.RData")
  
}

