

# Name : Hypervolume analysis
# Authors : Pauline Eymard-Dauphin & Diane Espel
# Objectives ## creating hypervolumes and computing metrics to compare them
#----------------------------------------------------------------------------

# Hypervolume analyses can take a long time to run (hours - days). 

# This script is set up for 3 hypervolume analyses, depending on the way data
# is pooled : by altitude, by morphotype or by sex.
# You can find the beginning of each one by searching for "wantAnalysisByAltitude"
# (or -ByMorphotype or -BySex).
# The first part of the script prepares the data and should be run in any case.




# Required packages -------------------------------------------------------

library(ggplot2)
library(GGally)
library(viridis)
library(hypervolume)
library(cati)
library(mice)
library(e1071)
library(tidyverse)
library(mgcv)
library(fmsb)
library(cowplot)
library(reshape2)


# defining colors manually
viridis_col=c("#39568CFF","#1F968BFF","#29AF7FFF","#73D055FF","#B8DE29FF")
viridis_col2=c("#1F968BFF","#73D055FF")
viridis_col3=c("#B8DE29FF","#39568CFF") #light morph in green, dark morph in blue


# Z_test function, p_value adjusted ---------------------------------------------
ztest_function <- function(x,y){
  z_stat <- (mean(x) - mean(y)) / sqrt(sd(x)^2 / ind + sd(y)^2 / ind)
  p_value <- 2*pnorm(-abs(z_stat)) 
  p_adjust <- p.adjust(p_value, method = "bonferroni", n=num_comp)  
  return(p_adjust)
}

# Open data  -------------------------------------------------------------

DATA <- read.csv2("data/Altitudinal transects Amblystogenium pacificum.csv")

DATA<- na.omit(DATA)
data=DATA[,c(3,5:7,15,16,18,20)] #remove useless columns


# Exploring data to remove highly correlated variables --------------------

# PairPlot by altitude
ggpairs(data[,c(4:dim(data)[2])],aes(color=data$Altitude))
# PairPlot by morphotype
ggpairs(data[,c(4:dim(data)[2])],aes(color=data$Morphotype))
# PairPlot by sex
ggpairs(data[,c(4:dim(data)[2])],aes(color=data$Sex))


# Data preparation --------------------

# Select traits
trait_axes <- c("Mass","Body_size_index_corr","Proteines.Mass","Lipids.Mass","Sugars.Mass")
num_trait = length(trait_axes)


# Center and scale numerical data
for (i in 4:dim(data)[2]){ #scaling by columns [2]
  data[,i]<-scale(data[,i], center=TRUE, scale=TRUE)
}



########################################################################################
# Hypervolume analysis 
########################################################################################

wantAnalysisByAltitude=F
if(wantAnalysisByAltitude==TRUE){
  
  # Check the number of individuals per hypervolume (i.e. per altitude)
  table(data$Altitude)
  
  
  # Building hypervolumes ---------------------------------------------------
  
  # list the hypervolumes to be calculated
  cat_list_altitude = sort(as.character(unique(data$Altitude)))
  num_cat_altitude = length(cat_list_altitude)
  print(paste0("The number of altitude categories is : ", num_cat_altitude))
  
  # Calculate hypervolumes (1 per altitude) and join all of them 
  # (takes about 1 min to run the loop)
  hv_list_altitude = new("HypervolumeList") # list that will contain all the hypervolumes
  hv_list_altitude@HVList = vector(mode="list",length=num_cat_altitude)
  for (i in 1:num_cat_altitude){ 
    # select data
    data_ce_cat_altitude = data[data$Altitude==cat_list_altitude[i],trait_axes]
    # Calculate hypervolume for each altitude
    hv_list_altitude@HVList[[i]] <- hypervolume_gaussian(data_ce_cat_altitude,name=as.character(cat_list_altitude[i]),
                                                         verbose=FALSE,kde.bandwidth = estimate_bandwidth(data_ce_cat_altitude))
  }
  
  NOMpng="results/Hypervolume_Traits_by_altitude.png"
  png(file = NOMpng,width=400,height = 400)
  plot(hypervolume_join(hv_list_altitude@HVList[[1]], hv_list_altitude@HVList[[2]], hv_list_altitude[[3]], 
                        hv_list_altitude@HVList[[4]], hv_list_altitude@HVList[[5]]),
       show.density = T, show.data = T, show.contour = F, cex.legend = 1.6, cex.names = 1.6, col= viridis_col, 
       names = trait_axes)
  dev.off()
  
  
  
  
  
  # Assessing pairwise similatity between hypervolumes ----------------------------------
  
  # Sorensen similarity index
  # this loop takes about 4 mins to run
  overlap_altitude = matrix(NA, nrow=num_cat_altitude, ncol=num_cat_altitude, dimnames = list(cat_list_altitude,cat_list_altitude)) #empty matrix
  for (i in 1:num_cat_altitude){
    for (j in 1:i){
      if (i!=j){
        # Assemble the hypervolumes two by two
        this_set = hypervolume_set(hv_list_altitude@HVList[[i]], hv_list_altitude@HVList[[j]], check.memory=FALSE)
        # Calculate Sorensen
        overlap_altitude[i,j] = hypervolume_overlap_statistics(this_set)["sorensen"]
      }  }  } 
  

  
  # Centroid distances
  distance_altitude = matrix(NA, nrow = num_cat_altitude,ncol = num_cat_altitude, dimnames = list(cat_list_altitude,cat_list_altitude))
  for(i in 1:num_cat_altitude){
    for(j in 1:i){
      if (i!=j){
        distance_altitude[i,j]<-hypervolume_distance(hv_list_altitude@HVList[[i]], hv_list_altitude@HVList[[j]],check.memory = FALSE)
      } } }
  
  
  OverlapFile="results/overlap_altitude.csv"
  write.table(overlap_altitude,OverlapFile,sep=";",dec = ".")
  distanceFile=paste0(save_path,"distance_altitude.csv")
  write.table(distance_altitude,distanceFile,sep=";",dec = ".")
  

  
    # Volumes, contributions, centroids observed #### ---------------------------------------------------
  
  # Get volume 
  #Volume is a measure of hypervolume size and represents the width of the multidimensional,
  # trait space of the altitudes (i.e. the variability of all traits shaping the hypervolume simultaneously).
  volume_altitude<-as.data.frame(get_volume(hv_list_altitude))
  
  # Get profile trait contribution and plot them 
  # this works through a take-one-out method : consider trait T, calculate hypervolume without
  # T, compare volumes with and without. Repeat for all traits.
  # this loop takes about 7 mins to run
  forme_altitude=matrix(0,num_cat_altitude,num_trait, dimnames = list(cat_list_altitude,trait_axes))
  for (i in 1:num_cat_altitude){
    forme_altitude[i,]<-hypervolume_variable_importance(hv_list_altitude@HVList[[i]])
  }
  
  plotforme_altitude<-as.data.frame(forme_altitude)
  plotforme_altitude<- rbind(rep(round(max(plotforme_altitude)),num_trait) , rep(0,num_trait) ,plotforme_altitude) # Set plot limits
  
  # trait contribution radar charts
  NOMpng="results/Hypervolume_Traits_by_altitude_profile.png" 
  png(file = NOMpng)
  par(mar=rep(0.8,4))
  par(mfrow=c(2,3)) # To adapt 
  for(i in 1:num_cat_altitude){ 
    radarchart(plotforme_altitude[c(1,2,i+2),], axistype=2, title=cat_list_altitude[i], maxmin = T,
               plwd=3, plty=1, cglcol="grey", cglty=1, cglwd=1, axislabcol="black", vlcex=0.9)
  }
  dev.off()

  
  # Get centroids 
  centroides_altitude<-as.data.frame(get_centroid(hv_list_altitude))
  centroides_altitude$Altitude <- rownames(centroides_altitude)
  
  
  
#-------------------------------------------------------------------------------------------- 
  #  Comparing volumes (i.e. volume of each cloud) between altitudes ---------------------------------------------------
#-------------------------------------------------------------------------------------------- 
  
  # Bootstrapping volumes within each altitude level to compare them
  # this takes a long time to run (hours)
  
  # Because we could only calculate 1 hypervolume/altitude,
  # (hypervolume calculation requires at least one individual per each trait involved in hypervolume shaping), 
  # we lacked variability for comparing volumes between levels.
  # To solve this problem, we simulated 100 hypervolumes per altitude using bootstrapping
  
  table(data$Altitude) # see number of individuals by altitude
  sim=100 # Number of simulation
  ind=21 # Number of individual = Number of individuals for the case where there are the least
  
  
  # "Because the volume depends on the number of values used to build the hypervolume (Lamanna et al.,2014), 
  # we set the number of values used to build the hypervolume (i.e. number of individuals) at 21 for each simulation, 
  # which was the number of values obtained in the considered factor with the least sampled level
  # (i.e., no-competition)."
  
  volumeBS_altitude=matrix(0,num_cat_altitude,sim) # Row = Altitude Col = resampling
  rownames(volumeBS_altitude)=cat_list_altitude
  for (i in 1:num_cat_altitude){  
    data_ce_cat_altitude = data[data$Altitude==cat_list_altitude[i],trait_axes] # Data for this altitude
    for(j in 1:sim){
      
      dataBS_altitude = data_ce_cat_altitude[sample(nrow(data_ce_cat_altitude),ind,replace=T),] # sampling within the altitude
      hyperBS_altitude = hypervolume_gaussian(dataBS_altitude,verbose=FALSE, kde.bandwidth = estimate_bandwidth(dataBS_altitude))
      volumeBS_altitude[i,j]<-get_volume(hyperBS_altitude)
    }
    print(Sys.time())}
  
  ## Boxplot of resampled (i.e. boostraped) volumes 
  # We compared the volume of different sites, with all the hypervolumes built from 
  # the same number of values. 
  # We compared the variability of the hypervolumes built for each altitude against 
  # the variability of a hypervolume built using individuals from all sites at the same 
  # time
  
  plotvolumeBS_altitude<-as.data.frame(volumeBS_altitude)
  plotvolumeBS_altitude$Altitude <- rownames(plotvolumeBS_altitude)
  plotvolumeBS_altitude<-melt(plotvolumeBS_altitude,  id="Altitude")
  
  ggplot(plotvolumeBS_altitude, aes(Altitude, value, fill=Altitude))+
    geom_boxplot() + scale_fill_viridis_d(begin=0.3,end=0.9)+ theme_bw()


  ## Ztest to test significant difference between volumes
  ztest_volume_altitude <- matrix(NA,num_cat_altitude,num_cat_altitude, dimnames = list(cat_list_altitude,cat_list_altitude))
  num_comp = num_cat_altitude*(num_cat_altitude-1)/2 # Number of comparisons 
  for(i in 1:num_cat_altitude){ 
    for(j in 1:i){
      if (i!=j){
        ztest_volume_altitude[i,j]<- ztest_function(plotvolumeBS_altitude[plotvolumeBS_altitude$Altitude==cat_list_altitude[i],"value"],
                                                    plotvolumeBS_altitude[plotvolumeBS_altitude$Altitude==cat_list_altitude[j],"value"])
      } } }
  
  # Plot resampled volumes (i.e. boostraped volumes) : mean, sd, letters for differences
  volumeztest_altitude<-as.data.frame(tapply(plotvolumeBS_altitude$value, plotvolumeBS_altitude$Altitude, mean))
  colnames(volumeztest_altitude)="mean"
  volumeztest_altitude$sd<-tapply(plotvolumeBS_altitude$value, plotvolumeBS_altitude$Altitude, sd)
  volumeztest_altitude$Altitude<-rownames(volumeztest_altitude)
  volumeztest_altitude$lowersd <- volumeztest_altitude$mean-volumeztest_altitude$sd
  volumeztest_altitude$uppersd <- volumeztest_altitude$mean+volumeztest_altitude$sd
  volumeztest_altitude$letter <- c("a","ab","b","ab","ab")## Put the letters manually according to the significant differences in ztest_volume_altitude
  
  NOMpng="results/Hypervolume_Traits_by_altitude_volumeZtest.png" 
  png(file = NOMpng,width = 400, height =300)
  ggplot(volumeztest_altitude,aes(x= Altitude, y=mean, color=Altitude,label = letter)) +
    geom_point(shape  = 16, size   = 6) + 
    geom_errorbar(aes(ymin  = lowersd, ymax  =  uppersd), width =  0.1,size  =  1.2) +
    theme_bw() + ylab(expression(Volume~~(SD^5))) + ggtitle ("") +
    scale_color_viridis_d(begin=0.3,end=0.9) + 
    ylim(0,500)+
    geom_text(aes(y = 150), color  = "black") +
    theme(text= element_text(size=18), axis.title.x = element_blank(), 
          plot.title = element_text(size=18, face="italic"),legend.position = "none")
  
  dev.off()
  
  #-------------------------------------------------------------------------------------------- 
  # Comparison of trait contributions between altitudes and the null model over all 
  # altitudes
  #-------------------------------------------------------------------------------------------- 
  
  # Create a matrix containing the lower bounds, a matrix containing the upper bounds 
  # of the CI and then compare with the observed values
  
  Forme_bornes_inf_altitude=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Forme_bornes_inf_altitude[1,i]<-quantile(formeBS_altitude[,i], 0.025)
  }
  
  Forme_bornes_sup_altitude=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Forme_bornes_sup_altitude[1,i]<-quantile(formeBS_altitude[,i], 0.975)
  }
  
  Forme_IC_altitude=matrix(NA,num_cat_altitude,num_trait, dimnames = list(cat_list_altitude,trait_axes))
  for (i in 1:num_cat_altitude){
    for (j in 1:num_trait){
      if (forme_altitude[i,j]<Forme_bornes_inf_altitude[1,j]){
        Forme_IC_altitude[i,j]<- "-"
        } else if (forme_altitude[i,j]>Forme_bornes_sup_altitude[1,j]){
        Forme_IC_altitude[i,j]<- "+"
        } else {
        Forme_IC_altitude[i,j]<- "ns"
        }
      }
  }
  
  # radarplot : CI inferior and posterior need to be in the same table as the observed data;
  # and need to be the 1st 2 lines (see ?radarchart):

  forme_altitude$Altitude <- rownames(as.data.frame(forme_altitude))
  radplot_tab <- rbind(Forme_bornes_sup_altitude, Forme_bornes_inf_altitude, forme_altitude)
  
  plotformeIC_altitude<-as.data.frame(forme_altitude)
  plotformeIC_altitude$Altitude<-rownames(plotformeIC_altitude)
  plotformeIC_altitude$type<-rep("obs",num_cat_altitude)
  Forme_IC_altitude=as.data.frame(Forme_IC_altitude)
  Forme_IC_altitude$Altitude<-rownames(Forme_IC_altitude)
  Forme_bornes_inf_altitude=as.data.frame(Forme_bornes_inf_altitude)
  Forme_bornes_inf_altitude$Altitude<-rownames(Forme_bornes_inf_altitude)
  Forme_bornes_inf_altitude$type<-"inf"
  Forme_bornes_sup_altitude=as.data.frame(Forme_bornes_sup_altitude)
  Forme_bornes_sup_altitude$Altitude<-rownames(Forme_bornes_sup_altitude)
  Forme_bornes_sup_altitude$type<-"sup"
  plotformeIC_altitude<-bind_rows(plotformeIC_altitude,Forme_bornes_inf_altitude,Forme_bornes_sup_altitude)
  plotformeIC_altitude <-rbind(plotformeIC_altitude, rep(0,num_trait)) # Limit inf
  plotformeIC_altitude <-rbind(plotformeIC_altitude, rep(round(max(plotformeIC_altitude[,1:num_trait])),num_trait)) 
  
  colnames(plotformeIC_altitude)<- c(trait_axes,"Altitude","type")
  
  NOMpng=paste0("results/Hypervolume_Traits_by_altitude_radarchart.png") # Dendrogramm of Sorensen index
  png(file = NOMpng)
  par(mar=rep(0.8,4))
  par(mfrow=c(1,1)) #Adapt
  coul <- c("grey", "grey", viridis_col) 
  coul_leg <- c("white", "white", viridis_col)
  coul_in <- c("grey","white",rep(alpha("white",0.1),num_cat_altitude))
  radarchart(plotformeIC_altitude[c(9,8,7,6,1:5),1:num_trait], axistype=6, title="", maxmin = T,
             pcol=coul, plwd=2, plty=1, pfcol = coul_in, cglcol="grey", cglty=2, cglwd=1, axislabcol="black", vlcex=1)
  legend(x="bottomleft", legend = c("","",cat_list_altitude), 
         bty = "n", pch=20 , col= coul_leg , text.col = "black", cex=1, pt.cex=3, horiz = F)
  
  dev.off()
  
  # to avoid re-running : store data
  save.image("results/insect_hypervolumes_altitude.RData")
}




# ---------------------------------------------------------------------------------
# MORPHOTYPES
# ---------------------------------------------------------------------------------

wantAnalysisByMorphotype=F
if(wantAnalysisByMorphotype==TRUE){
  
  # Check the number of individual per hypervolume (i.e. per morphotype)
  table(data$Morphotype)
  
  
  # Building hypervolumes ---------------------------------------------------
  
  # list the hypervolumes to be calculated
  cat_list_morphotype = sort(as.character(unique(data$Morphotype)))
  num_cat_morphotype = length(cat_list_morphotype)
  print(paste0("The number of altitude categories is : ", num_cat_morphotype))
  
  # Calculate hypervolumes (1 per morphotype) and join all of them
  hv_list_morphotype = new("HypervolumeList") # list that will contain all the hypervolumes
  hv_list_morphotype@HVList = vector(mode="list",length=num_cat_morphotype)
  for (i in 1:num_cat_morphotype){ 
    # select data
    data_ce_cat_morphotype = data[data$Morphotype==cat_list_morphotype[i],trait_axes]
    # Calculate hypervolume for each altitude
    hv_list_morphotype@HVList[[i]] <- hypervolume_gaussian(data_ce_cat_morphotype,name=as.character(cat_list_morphotype[i]),
                                                           verbose=FALSE,kde.bandwidth = estimate_bandwidth(data_ce_cat_morphotype))
  }
  
  NOMpng="results/Hypervolume_Traits_by_morphotype.png"
  png(file = NOMpng)
  plot(hypervolume_join(hv_list_morphotype@HVList[[1]], hv_list_morphotype@HVList[[2]]),
       show.density = T, show.data = T, show.contour = F, cex.legend = 1.6, cex.names = 1.6, col= viridis_col3, 
       names = trait_axes)
  dev.off()

  
  # Assessing similatity between hypervolumes ---------------------------------------------------
  
  # Sorensen similarity index
  overlap_morphotype = matrix(NA, nrow=num_cat_morphotype, ncol=num_cat_morphotype, dimnames = list(cat_list_morphotype,cat_list_morphotype)) #empty matrix
  #overlap_morphotype_CI = matrix(NA, nrow=num_cat_morphotype, ncol=num_cat_morphotype, dimnames = list(cat_list_morphotype,cat_list_morphotype)) #empty matrix
  for (i in 1:num_cat_morphotype){
    for (j in 1:i){
      if (i!=j){
        # Assemble the hypervolumes two by two
        this_set = hypervolume_set(hv_list_morphotype@HVList[[i]], hv_list_morphotype@HVList[[j]], check.memory=FALSE)
        # Calculate Sorensen
        overlap_morphotype[i,j] = hypervolume_overlap_statistics(this_set)["sorensen"]
      }  }  } 
  
  
  # Centroid distances
  distance_morphotype = matrix(NA, nrow = num_cat_morphotype,ncol = num_cat_morphotype, dimnames = list(cat_list_morphotype,cat_list_morphotype))
  for(i in 1:num_cat_morphotype){
    for(j in 1:i){
      if (i!=j){
        distance_morphotype[i,j]<-hypervolume_distance(hv_list_morphotype@HVList[[i]], hv_list_morphotype@HVList[[j]],check.memory = FALSE)
      } } }
  
  
  
  # Volumes, contributions, centro?ds observed #### ---------------------------------------------------
  
  # Get volume 
  volume_morphotype<-as.data.frame(get_volume(hv_list_morphotype))
  
  # Get profile trait contribution and plot them 
  forme_morphotype=matrix(0,num_cat_morphotype,num_trait, dimnames = list(cat_list_morphotype,trait_axes))
  for (i in 1:num_cat_morphotype){
    forme_morphotype[i,]<-hypervolume_variable_importance(hv_list_morphotype@HVList[[i]])
  }
  
  plotforme_morphotype<-as.data.frame(forme_morphotype)
  plotforme_morphotype<- rbind(rep(round(max(plotforme_morphotype)),num_trait) , rep(0,num_trait) ,plotforme_morphotype) # Set plot limits
  
  NOMpng="results/Hypervolume_Traits_by_morphotype_profile.png" 
  png(file = NOMpng)
  par(mar=rep(0.8,4))
  par(mfrow=c(1,2)) # To adapt 
  for(i in 1:num_cat_morphotype){ 
    radarchart(plotforme_morphotype[c(1,2,i+2),], axistype=2, title=cat_list_morphotype[i], maxmin = T,
               plwd=3, plty=1, cglcol="grey", cglty=1, cglwd=1, axislabcol="black", vlcex=0.9)
  }
  dev.off()

  # Get centroids 
  centroides_morphotype<-as.data.frame(get_centroid(hv_list_morphotype))
  centroides_morphotype$Morphotype <- rownames(centroides_morphotype)
  
  
  #  Comparing volumes (i.e. volume of each cloud) between morphotypes ---------------------------------------------------
  
  # Bootstraping volumes within each altitude to compare them
  table(data$Morphotype) # see number of individuals by altitude
  sim=100 # Number of simulation
  ind=62 # Number of individual = Number of individuals for the case where there are the least

  volumeBS_morphotype=matrix(0,num_cat_morphotype,sim) # Row = Altitude Col = resampling
  rownames(volumeBS_morphotype)=cat_list_morphotype
  for (i in 1:num_cat_morphotype){  
    data_ce_cat_morphotype = data[data$Morphotype==cat_list_morphotype[i],trait_axes] # Data for this altitude
    for(j in 1:sim){
      
      dataBS_morphotype = data_ce_cat_morphotype[sample(nrow(data_ce_cat_morphotype),ind,replace=T),] # sampling within the altitude
      hyperBS_morphotype = hypervolume_gaussian(dataBS_morphotype,verbose=FALSE, kde.bandwidth = estimate_bandwidth(dataBS_morphotype))
      volumeBS_morphotype[i,j]<-get_volume(hyperBS_morphotype)
    }}
  
  ## Boxplot of resampled (i.e. boostraped) volumes 
  plotvolumeBS_morphotype<-as.data.frame(volumeBS_morphotype)
  plotvolumeBS_morphotype$Morphotype <- rownames(plotvolumeBS_morphotype)
  plotvolumeBS_morphotype<-melt(plotvolumeBS_morphotype,  id="Morphotype")
  
  NOMpng="results/Hypervolume_Traits_by_morphotype_volumeBS.png"
  png(file = NOMpng,width = 400, height =300)
  ggplot(plotvolumeBS_morphotype, aes(Morphotype, value, fill=Morphotype))+
    geom_boxplot() + scale_fill_viridis_d(begin=0.9,end=0.3)+ theme_bw()
  dev.off()
 
  ## Ztest to test significant difference between volumes
  ztest_volume_morphotype <- matrix(NA,num_cat_morphotype,num_cat_morphotype, dimnames = list(cat_list_morphotype,cat_list_morphotype))
  num_comp = num_cat_morphotype*(num_cat_morphotype-1)/2 # Number of comparisons 
  for(i in 1:num_cat_morphotype){ 
    for(j in 1:i){
      if (i!=j){
        ztest_volume_morphotype[i,j]<- ztest_function(plotvolumeBS_morphotype[plotvolumeBS_morphotype$Morphotype==cat_list_morphotype[i],"value"],
                                                      plotvolumeBS_morphotype[plotvolumeBS_morphotype$Morphotype==cat_list_morphotype[j],"value"])
      } } }
  
  # Plot resampled volumes (i.e. boostraped volumes) : mean, sd, letters for differences
  volumeztest_morphotype<-as.data.frame(tapply(plotvolumeBS_morphotype$value, plotvolumeBS_morphotype$Morphotype, mean))
  colnames(volumeztest_morphotype)="mean"
  volumeztest_morphotype$sd<-tapply(plotvolumeBS_morphotype$value, plotvolumeBS_morphotype$Morphotype, sd)
  volumeztest_morphotype$Morphotype<-rownames(volumeztest_morphotype)
  volumeztest_morphotype$lowersd <- volumeztest_morphotype$mean-volumeztest_morphotype$sd
  volumeztest_morphotype$uppersd <- volumeztest_morphotype$mean+volumeztest_morphotype$sd
  volumeztest_morphotype$letter <- c("a","b")## Put the letters manually according to the significant differences in ztest_volume_morphotype
  
  NOMpng="results/Hypervolume_Traits_by_morphotype_volumeZtest.png"
  png(file = NOMpng,width = 400, height =300)
  ggplot(volumeztest_morphotype,aes(x=Morphotype, y=mean, color=Morphotype,label = letter)) +
    geom_point(shape  = 16, size   = 6) + 
    geom_errorbar(aes(ymin  = lowersd, ymax  =  uppersd), width =  0.1,size  =  1.2) +
    theme_bw() + ylab(expression(Volume~~(SD^5))) + ggtitle ("") +
    scale_color_viridis_d(begin=0.9,end=0.3) + 
    geom_text(aes(y = 150), color  = "black") +
    ylim(250,1500)+
    theme(text= element_text(size=18), axis.title.x = element_blank(), 
          plot.title = element_text(size=18, face="italic"),legend.position = "none")
  
  dev.off()
  

  # Comparing centroids between altitude ####---------------------------------------------------
  
  # Boostraping centroids
  table(data$Morphotype) 
  sim=100
  indi=62
  
  centroBS_morphotype=matrix(0,sim,num_trait) # Resampled hypervolume centroids Row = simulation, Col = traits
  formeBS_morphotype=matrix(0,sim,num_trait)  # Resampled hypervolume trait contributions Row = simulation, Col = traits
  for(j in 1:sim){
    
    dataBS_morphotype = data[sample(nrow(data),indi,replace=F),trait_axes] 
    hyperBS_morphotype = hypervolume_gaussian(dataBS_morphotype,verbose=FALSE, kde.bandwidth = estimate_bandwidth(dataBS_morphotype))
    centroBS_morphotype[j,(1:num_trait)]<-get_centroid(hyperBS_morphotype)
    formeBS_morphotype[j,(1:num_trait)] <- hypervolume_variable_importance(hyperBS_morphotype) 
  }
  
  #### Comparison of centroids between morphotypes and the null model over all morphotypes ####
  
  #Create a matrix containing the lower bounds, a matrix containing the upper bounds of the CI and 
  # then compare with the observed values
  Centroides_bornes_inf_morphotype=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Centroides_bornes_inf_morphotype[1,i]<-quantile(centroBS_morphotype[,i], 0.025)
  }
  
  Centroides_bornes_sup_morphotype=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Centroides_bornes_sup_morphotype[1,i]<-quantile(centroBS_morphotype[,i], 0.975)
  }
  
  
  # Plot the results
  library(cowplot)
  NOMpng="results/Hypervolume_Traits_by_morphotype_plotgrid_centroides.png"
  png(file = NOMpng)
  plot_grid(
    ggplot(centroides_morphotype, aes(x = Morphotype, y = Mass, color=Morphotype))+
      ggtitle("") + ylab("Mass")  + xlab("Morphotype") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_morphotype[1,1], 
                    ymax = Centroides_bornes_sup_morphotype[1,1]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.9,end=0.3) + theme_bw() +
      theme(axis.title.x = element_blank(), legend.position = "none", plot.title = element_text(size=12, face="bold.italic")),
    
    ggplot(centroides_morphotype, aes(x = Morphotype, y = Body_size_index_corr, color=Morphotype))+
      ggtitle("") + ylab("Body_size_index")  + xlab("Morphotype") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_morphotype[1,2], 
                    ymax = Centroides_bornes_sup_morphotype[1,2]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.9,end=0.3) + theme_bw() +
      theme(axis.title.x = element_blank(), legend.position = "none", plot.title = element_text(size=12, face="bold.italic")),
    
    ggplot(centroides_morphotype, aes(x = Morphotype, y = Proteines.Mass, color=Morphotype))+
      ggtitle("")+ ylab("Proteins:Mass") + xlab("Morphotype") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_morphotype[1,3], 
                    ymax = Centroides_bornes_sup_morphotype[1,3]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.9,end=0.3) + theme_bw() +
      theme(axis.title.x = element_blank(),legend.position = "none"),
    
    ggplot(centroides_morphotype, aes(x = Morphotype, y = Lipids.Mass, color=Morphotype))+
      ggtitle("")+ ylab("Lipids:Mass") + xlab("Morphotype") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_morphotype[1,4], 
                    ymax = Centroides_bornes_sup_morphotype[1,4]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.9,end=0.3)+ theme_bw() +
      theme(axis.title.x = element_blank(),legend.position = "none"),
    
    ggplot(centroides_morphotype, aes(x = Morphotype, y = Sugars.Mass, color=Morphotype))+
      ggtitle("")+ ylab("Sugars:Mass") + xlab("Morphotype") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_morphotype[1,5], 
                    ymax = Centroides_bornes_sup_morphotype[1,5]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.9,end=0.3)+ theme_bw() +
      theme(axis.title.x = element_blank(),legend.position = "none"),
    
    ncol = 3, nrow = 2, byrow = F)
  dev.off()
  

  ### Comparison of trait contributions between morphotypes and the null model over all morphotypes #### ---------------------------------------------------------------------
  
  ##Create a matrix containing the lower bounds, a matrix containing the upper bounds of the CI and then compare with the observed values
  
  Forme_bornes_inf_morphotype=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Forme_bornes_inf_morphotype[1,i]<-quantile(formeBS_morphotype[,i], 0.025)
  }
  
  Forme_bornes_sup_morphotype=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Forme_bornes_sup_morphotype[1,i]<-quantile(formeBS_morphotype[,i], 0.975)
  }
  
  Forme_IC_morphotype=matrix(NA,num_cat_morphotype,num_trait, dimnames = list(cat_list_morphotype,trait_axes))
  for (i in 1:num_cat_morphotype){
    for (j in 1:num_trait){
      Forme_IC_morphotype[i,j]<- if (forme_morphotype[i,j]<Forme_bornes_inf_morphotype[1,j]){"-"} 
      else if (forme_morphotype[i,j]>Forme_bornes_sup_morphotype[1,j]){"+"} else {"ns"}
    }}
  
  # There are certainly simpler than all these lines but it works...
  plotformeIC_morphotype<-as.data.frame(forme_morphotype)
  plotformeIC_morphotype$Morphotype<-rownames(plotformeIC_morphotype)
  plotformeIC_morphotype$type<-rep("obs",num_cat_morphotype)
  Forme_IC_morphotype=as.data.frame(Forme_IC_morphotype)
  Forme_IC_morphotype$Morphotype<-rownames(Forme_IC_morphotype)
  Forme_bornes_inf_morphotype=as.data.frame(Forme_bornes_inf_morphotype)
  Forme_bornes_inf_morphotype$Morphotype<-rownames(Forme_bornes_inf_morphotype)
  Forme_bornes_inf_morphotype$type<-"inf"
  Forme_bornes_sup_morphotype=as.data.frame(Forme_bornes_sup_morphotype)
  Forme_bornes_sup_morphotype$Morphotype<-rownames(Forme_bornes_sup_morphotype)
  Forme_bornes_sup_morphotype$type<-"sup"
  plotformeIC_morphotype<-bind_rows(plotformeIC_morphotype,Forme_bornes_inf_morphotype,Forme_bornes_sup_morphotype)
  plotformeIC_morphotype <-rbind(plotformeIC_morphotype, rep(0,num_trait)) # Limit inf
  plotformeIC_morphotype <-rbind(plotformeIC_morphotype, rep(round(max(plotformeIC_morphotype[,1:num_trait])),num_trait)) 
  
  colnames(plotformeIC_morphotype)<- c(trait_axes,"Morphotype","type")
  
  NOMpng="results/Hypervolume_Traits_by_morphotype_radarchart.png"
  png(file = NOMpng)
  par(mar=rep(0.8,4))
  par(mfrow=c(1,1)) #Adapt
  coul <- c("grey", "grey", viridis_col3) 
  coul_leg <- c("white", "white", viridis_col3)
  coul_in <- c("grey","white",rep(alpha("white",0.1),num_cat_morphotype))
  radarchart(plotformeIC_morphotype[c(6,5,4,3, 1:2),1:num_trait], axistype=6, title="", maxmin = T,
             pcol=coul, plwd=2, plty=1, pfcol = coul_in, cglcol="grey", cglty=2, cglwd=1, axislabcol="black", vlcex=1)
  legend(x="bottomleft", legend = c("","",cat_list_morphotype), 
         bty = "n", pch=20 , col= coul_leg , text.col = "black", cex=1, pt.cex=3, horiz = F)
  
  dev.off()
  
  
  save.image("results/insect_hypervolumes_morphotype.RData")
  
}


wantAnalysisBySex=T
if(wantAnalysisBySex==TRUE){
  
  # Check the number of individual per hypervolume (i.e. per morphotype)
  table(data$Sex)
  
  
  # Building hypervolumes ---------------------------------------------------
  
  # list the hypervolumes to be calculated
  cat_list_sex = sort(as.character(unique(data$Sex)))
  num_cat_sex = length(cat_list_sex)
  print(paste0("The number of altitude categories is : ", num_cat_sex))
  
  # Calculate hypervolumes (1 per morphotype) and join all of them
  hv_list_sex = new("HypervolumeList") # list that will contain all the hypervolumes
  hv_list_sex@HVList = vector(mode="list",length=num_cat_sex)
  for (i in 1:num_cat_sex){ 
    # select data
    data_ce_cat_sex = data[data$Sex==cat_list_sex[i],trait_axes]
    # Calculate hypervolume for each altitude
    hv_list_sex@HVList[[i]] <- hypervolume_gaussian(data_ce_cat_sex,name=as.character(cat_list_sex[i]),
                                                    verbose=FALSE,kde.bandwidth = estimate_bandwidth(data_ce_cat_sex))
  }
  
  NOMpng="results/Hypervolume_Traits_by_sex.png"
  png(file = NOMpng,width=400,height=400)
  plot(hypervolume_join(hv_list_sex@HVList[[1]], hv_list_sex@HVList[[2]]),
       show.density = T, show.data = T, show.contour = F, cex.legend = 1.6, cex.names = 1.6, col= viridis_col3, 
       names = trait_axes)
  dev.off()
  
 
  # Assessing similatity between hypervolumes ---------------------------------------------------
  
  # Sorensen similarity index
  overlap_sex = matrix(NA, nrow=num_cat_sex, ncol=num_cat_sex, dimnames = list(cat_list_sex,cat_list_sex)) #empty matrix
  #overlap_sex_CI = matrix(NA, nrow=num_cat_sex, ncol=num_cat_sex, dimnames = list(cat_list_sex,cat_list_sex)) #empty matrix
  for (i in 1:num_cat_sex){
    for (j in 1:i){
      if (i!=j){
        # Assemble the hypervolumes two by two
        this_set = hypervolume_set(hv_list_sex@HVList[[i]], hv_list_sex@HVList[[j]], check.memory=FALSE)
        # Calculate Sorensen
        overlap_sex[i,j] = hypervolume_overlap_statistics(this_set)["sorensen"]
        #overlap_sex_CI[i,j] = hypervolume_overlap_confidence(this_set)["sorensen"]
      }  }  } #computing
  
  
  # Centro?d distances
  distance_sex = matrix(NA, nrow = num_cat_sex,ncol = num_cat_sex, dimnames = list(cat_list_sex,cat_list_sex))
  for(i in 1:num_cat_sex){
    for(j in 1:i){
      if (i!=j){
        distance_sex[i,j]<-hypervolume_distance(hv_list_sex@HVList[[i]], hv_list_sex@HVList[[j]],check.memory = FALSE)
      } } }
  
  
  
  # Volumes, contributions, centro?ds observed #### ---------------------------------------------------
  
  # Get volume 
  # "Volume is a measure of hypervolume size and represents the width of the multidimensional,
  # trait space of the morphotypes (i.e. the variability of all traits shaping the hypervolume simultaneously)."
  volume_sex<-as.data.frame(get_volume(hv_list_sex))
  
  # Get profile trait contribution and plot them 
  forme_sex=matrix(0,num_cat_sex,num_trait, dimnames = list(cat_list_sex,trait_axes))
  for (i in 1:num_cat_sex){
    forme_sex[i,]<-hypervolume_variable_importance(hv_list_sex@HVList[[i]])
  }
  
  plotforme_sex<-as.data.frame(forme_sex)
  plotforme_sex<- rbind(rep(round(max(plotforme_sex)),num_trait) , rep(0,num_trait) ,plotforme_sex) # Set plot limits
  
  NOMpng="results/Hypervolume_Traits_by_sex_profile.png"
  png(file = NOMpng)
  par(mar=rep(0.8,4))
  par(mfrow=c(1,2)) # To adapt 
  for(i in 1:num_cat_sex){ 
    radarchart(plotforme_sex[c(1,2,i+2),], axistype=2, title=cat_list_sex[i], maxmin = T,
               plwd=3, plty=1, cglcol="grey", cglty=1, cglwd=1, axislabcol="black", vlcex=0.9)
  }
  dev.off()
  
  
  # Get centroids 
  centroides_sex<-as.data.frame(get_centroid(hv_list_sex))
  centroides_sex$Sex <- rownames(centroides_sex)
  
  
  #  Comparing volumes (i.e. volume of each cloud) between morphotypes ---------------------------------------------------
  
  # Bootstraping volumes within each altitude to compare them
  
  # "Because we could only calculate 1 hypervolume/(Sex,
  # (hypervolume calculation requires at least one individual per each trait involved in hypervolume shaping), 
  # we lacked variability for comparing volumes between levels.
  # To solve this problem, we simulated 100 hypervolumes per morphotype using bootstrapping"
  
  table(data$Sex) # see number of individuals by altitude
  sim=100 # Number of simulation
  ind=55 # Number of individual = Number of individuals for the case where there are the least
  # "Because the volume depends on the number of values used to build the hypervolume (Lamanna et al.,2014), 
  # we set the number of values used to build the hypervolume (i.e. number of individuals) at 22 for each simulation, 
  # which was the number of values obtained in the considered factor with the least sampled level
  # (i.e., no-competition)."

  volumeBS_sex=matrix(0,num_cat_sex,sim) # Row = Altitude Col = resampling
  rownames(volumeBS_sex)=cat_list_sex
  for (i in 1:num_cat_sex){  
    data_ce_cat_sex = data[data$Sex==cat_list_sex[i],trait_axes] # Data for this altitude
    for(j in 1:sim){
      
      dataBS_sex = data_ce_cat_sex[sample(nrow(data_ce_cat_sex),ind,replace=T),] # sampling within the altitude
      hyperBS_sex = hypervolume_gaussian(dataBS_sex,verbose=FALSE, kde.bandwidth = estimate_bandwidth(dataBS_sex))
      volumeBS_sex[i,j]<-get_volume(hyperBS_sex)
    }}
  
  ## Boxplot of resampled (i.e. boostraped) volumes 
  # "we compared the volume of different sites, with all the hypervolumes built from the same number of values. 
  #  we compared the variability of the hypervolumes built for each altitude against the variability of a hypervolume built using individuals from
  # all sites at the same time"
  
  plotvolumeBS_sex<-as.data.frame(volumeBS_sex)
  plotvolumeBS_sex$Sex <- rownames(plotvolumeBS_sex)
  plotvolumeBS_sex<-melt(plotvolumeBS_sex,  id="Sex")
  
  NOMpng="results/Hypervolume_Traits_by_sex_volumeBS.png"
  png(file = NOMpng,width = 400, height =300)
  ggplot(plotvolumeBS_sex, aes(Sex, value, fill=Sex))+
    geom_boxplot() + scale_fill_viridis_d(begin=0.9,end=0.3)+ theme_bw()
  dev.off()
  
   ## Ztest to test significant difference between volumes
  ztest_volume_sex <- matrix(NA,num_cat_sex,num_cat_sex, dimnames = list(cat_list_sex,cat_list_sex))
  num_comp = num_cat_sex*(num_cat_sex-1)/2 # Number of comparisons 
  for(i in 1:num_cat_sex){ 
    for(j in 1:i){
      if (i!=j){
        ztest_volume_sex[i,j]<- ztest_function(plotvolumeBS_sex[plotvolumeBS_sex$Sex==cat_list_sex[i],"value"],
                                               plotvolumeBS_sex[plotvolumeBS_sex$Sex==cat_list_sex[j],"value"])
      } } }
  
  # Plot resampled volumes (i.e. boostraped volumes) : mean, sd, letters for differences
  volumeztest_sex<-as.data.frame(tapply(plotvolumeBS_sex$value, plotvolumeBS_sex$Sex, mean))
  colnames(volumeztest_sex)="mean"
  volumeztest_sex$sd<-tapply(plotvolumeBS_sex$value, plotvolumeBS_sex$Sex, sd)
  volumeztest_sex$Sex<-rownames(volumeztest_sex)
  volumeztest_sex$lowersd <- volumeztest_sex$mean-volumeztest_sex$sd
  volumeztest_sex$uppersd <- volumeztest_sex$mean+volumeztest_sex$sd
  volumeztest_sex$letter <- c("a","b")## Put the letters manually according to the significant differences in ztest_volume_sex
  
  NOMpng="results/Hypervolume_Traits_by_sex_volumeZtest.png"
  png(file = NOMpng,width = 400, height =300)
  ggplot(volumeztest_sex,aes(x=Sex, y=mean, color=Sex,label = letter)) +
    geom_point(shape  = 16, size   = 6) + 
    geom_errorbar(aes(ymin  = lowersd, ymax  =  uppersd), width =  0.1,size  =  1.2) +
    theme_bw() + ylab(expression(Volume~~(SD^5))) + ggtitle ("") +
    scale_color_viridis_d(begin=0.9,end=0.3) + 
    geom_text(aes(y = 150), color  = "black") +
    ylim(250,1500)+
    theme(text= element_text(size=18), axis.title.x = element_blank(), 
          plot.title = element_text(size=18, face="italic"),legend.position = "none")
  
  dev.off()
  

  # Comparing centroids between altitude ####---------------------------------------------------
  
  # Boostraping centroids
  table(data$Sex) 
  sim=100
  indi=55
  
  centroBS_sex=matrix(0,sim,num_trait) # Resampled hypervolume centro?ds Row = simulation, Col = traits
  formeBS_sex=matrix(0,sim,num_trait)  # Resampled hypervolume trait contributions Row = simulation, Col = traits
  for(j in 1:sim){
    
    dataBS_sex = data[sample(nrow(data),indi,replace=F),trait_axes] 
    hyperBS_sex = hypervolume_gaussian(dataBS_sex,verbose=FALSE, kde.bandwidth = estimate_bandwidth(dataBS_sex))
    centroBS_sex[j,(1:num_trait)]<-get_centroid(hyperBS_sex)
    formeBS_sex[j,(1:num_trait)] <- hypervolume_variable_importance(hyperBS_sex) 
  }
  
  #### Comparison of centroids between morphotypes and the null model over all morphotypes ####
  
  #Create a matrix containing the lower bounds, a matrix containing the upper bounds of the CI and 
  # then compare with the observed values
  Centroides_bornes_inf_sex=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Centroides_bornes_inf_sex[1,i]<-quantile(centroBS_sex[,i], 0.025)
  }
  
  Centroides_bornes_sup_sex=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Centroides_bornes_sup_sex[1,i]<-quantile(centroBS_sex[,i], 0.975)
  }
  
  
  # Plot the results
  library(cowplot)
  NOMpng="results/Hypervolume_Traits_by_sex_plotgrid_centroides.png"
  png(file = NOMpng)
  plot_grid(
    ggplot(centroides_sex, aes(x = Sex, y = Mass, color=Sex))+
      ggtitle("") + ylab("Mass")  + xlab("Sex") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_sex[1,1], 
                    ymax = Centroides_bornes_sup_sex[1,1]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.9,end=0.3) + theme_bw() +
      theme(axis.title.x = element_blank(), legend.position = "none", plot.title = element_text(size=12, face="bold.italic")),
    
    ggplot(centroides_sex, aes(x = Sex, y = Body_size_index_corr, color=Sex))+
      ggtitle("") + ylab("Body_size_index")  + xlab("Sex") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_sex[1,2], 
                    ymax = Centroides_bornes_sup_sex[1,2]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.9,end=0.3) + theme_bw() +
      theme(axis.title.x = element_blank(), legend.position = "none", plot.title = element_text(size=12, face="bold.italic")),
    
    ggplot(centroides_sex, aes(x = Sex, y = Proteines.Mass, color=Sex))+
      ggtitle("")+ ylab("Proteins:Mass") + xlab("Sex") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_sex[1,3], 
                    ymax = Centroides_bornes_sup_sex[1,3]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.9,end=0.3) + theme_bw() +
      theme(axis.title.x = element_blank(),legend.position = "none"),
    
    ggplot(centroides_sex, aes(x = Sex, y = Lipids.Mass, color=Sex))+
      ggtitle("")+ ylab("Lipids:Mass") + xlab("Sex") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_sex[1,4], 
                    ymax = Centroides_bornes_sup_sex[1,4]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.9,end=0.3)+ theme_bw() +
      theme(axis.title.x = element_blank(),legend.position = "none"),
    
    ggplot(centroides_sex, aes(x = Sex, y = Sugars.Mass, color=Sex))+
      ggtitle("")+ ylab("Sugars:Mass") + xlab("Sex") + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Centroides_bornes_inf_sex[1,5], 
                    ymax = Centroides_bornes_sup_sex[1,5]), alpha = 0.1, color = "Grey")+ 
      geom_point(size=3) + scale_color_viridis_d(begin=0.9,end=0.3)+ theme_bw() +
      theme(axis.title.x = element_blank(),legend.position = "none"),
    
    ncol = 3, nrow = 2, byrow = F)
  dev.off()
  
 
  ### Comparison of trait contributions between morphotypes and the null model over all morphotypes #### ---------------------------------------------------------------------
  
  ##Create a matrix containing the lower bounds, a matrix containing the upper bounds of the CI and then compare with the observed values
  
  Forme_bornes_inf_sex=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Forme_bornes_inf_sex[1,i]<-quantile(formeBS_sex[,i], 0.025)
  }
  
  Forme_bornes_sup_sex=matrix(0,1,num_trait, dimnames = list("0",trait_axes))
  for (i in 1:num_trait){
    Forme_bornes_sup_sex[1,i]<-quantile(formeBS_sex[,i], 0.975)
  }
  
  Forme_IC_sex=matrix(NA,num_cat_sex,num_trait, dimnames = list(cat_list_sex,trait_axes))
  for (i in 1:num_cat_sex){
    for (j in 1:num_trait){
      Forme_IC_sex[i,j]<- if (forme_sex[i,j]<Forme_bornes_inf_sex[1,j]){"-"} 
      else if (forme_sex[i,j]>Forme_bornes_sup_sex[1,j]){"+"} else {"ns"}
    }}
  
  # There are certainly simpler than all these lines but it works...
  plotformeIC_sex<-as.data.frame(forme_sex)
  plotformeIC_sex$Sex<-rownames(plotformeIC_sex)
  plotformeIC_sex$type<-rep("obs",num_cat_sex)
  Forme_IC_sex=as.data.frame(Forme_IC_sex)
  Forme_IC_sex$Sex<-rownames(Forme_IC_sex)
  Forme_bornes_inf_sex=as.data.frame(Forme_bornes_inf_sex)
  Forme_bornes_inf_sex$Sex<-rownames(Forme_bornes_inf_sex)
  Forme_bornes_inf_sex$type<-"inf"
  Forme_bornes_sup_sex=as.data.frame(Forme_bornes_sup_sex)
  Forme_bornes_sup_sex$Sex<-rownames(Forme_bornes_sup_sex)
  Forme_bornes_sup_sex$type<-"sup"
  plotformeIC_sex<-bind_rows(plotformeIC_sex,Forme_bornes_inf_sex,Forme_bornes_sup_sex)
  plotformeIC_sex <-rbind(plotformeIC_sex, rep(0,num_trait)) # Limit inf
  plotformeIC_sex <-rbind(plotformeIC_sex, rep(round(max(plotformeIC_sex[,1:num_trait])),num_trait)) 
  
  colnames(plotformeIC_sex)<- c(trait_axes,"Sex","type")
  
  NOMpng="results/Hypervolume_Traits_by_sex_radarchart.png"
  png(file = NOMpng)
  par(mar=rep(0.8,4))
  par(mfrow=c(1,1)) #Adapt
  coul <- c("grey", "grey", viridis_col3) 
  coul_leg <- c("white", "white", viridis_col3)
  coul_in <- c("grey","white",rep(alpha("white",0.1),num_cat_sex))
  radarchart(plotformeIC_sex[c(6,5,4,3, 1:2),1:num_trait], axistype=6, title="", maxmin = T,
             pcol=coul, plwd=2, plty=1, pfcol = coul_in, cglcol="grey", cglty=2, cglwd=1, axislabcol="black", vlcex=1)
  legend(x="bottomleft", legend = c("","",cat_list_sex), 
         bty = "n", pch=20 , col= coul_leg , text.col = "black", cex=1, pt.cex=3, horiz = F)
  
  dev.off()
}


OverlapFile="results/overlap_altitude.csv"
write.table(overlap_altitude,OverlapFile,sep=";",dec = ".")

OverlapFile="results/overlap_morphotype.csv"
write.table(overlap_morphotype,OverlapFile,sep=";",dec = ".")

OverlapFile="results/overlap_sex.csv"
write.table(overlap_sex,OverlapFile,sep=";",dec = ".")

distanceFile="results/distance_altitude.csv"
write.table(distance_altitude,distanceFile,sep=";",dec = ".")

distanceFile="results/distance_morphotype.csv"
write.table(distance_morphotype,distanceFile,sep=";",dec = ".")

distanceFile="results/distance_sex.csv"
write.table(distance_sex,distanceFile,sep=";",dec = ".")
