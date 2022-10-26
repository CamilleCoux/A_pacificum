
# Name : FAMD analysis
# Authors : Diane ESPEL
# Objectives ## 
#----------------------------------------------------------------------------

# R version : 2.15.2


# Required packages -------------------------------------------------------

library(FactoMineR)
library(readxl)
library(factoextra)
library(corrplot)
library(viridis)


# Open data  -------------------------------------------------------------

setwd(open_path)
DATA <- read.csv2("data/Altitudinal transects Amblystogenium pacificum.csv")


# Check missing rows
info=skimr::skim(DATA)
print(info)


DATA<- na.omit(DATA)
tab=DATA[,c(3,5:7,15,16,18,20)] # keep only meaningful columns
tab$Altitude=gsub("m", "", tab$Altitude) # remove m
tab$Altitude=as.numeric(tab$Altitude)
str(tab$Altitude)

# rename some levels
tab$Morphotype=gsub("B","L-morphotype",tab$Morphotype)
tab$Morphotype=gsub("N","D-morphotype",tab$Morphotype)
tab$Sex=gsub("M","Male",tab$Sex)
tab$Sex=gsub("F","Female",tab$Sex)

# Apply AFDM   -------------------------------------------------------------

res_famd=FAMD(tab,  #dateframe with n rows and p columns
              ncp = 5, # maximum number of dimensions
              sup.var = NULL, # vector indicating positions of supplementary variables 
              ind = NULL,  # vector indicating positions of supplementary variables 
              graph = T# display the graph
)



# Variances  -------------------------------------------------------------

# get eigen values
eig.val=res_famd$eig 
eig.val=as.data.frame(eig.val)
print(eig.val)


# Plot of Eigen values
fviz_eig(res_famd,title="Eigenvalues")+ 
  theme(text = element_text(size = 17),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 15))


# Graphs of all variables  -------------------------------------------------------------

#  get results for variables
var=get_famd_var(res_famd) # for both quantitative and qualitative variables
var$coord # All variables coordinates

var$cos2 # Quality of the representation
Quality=var$cos2
Quality=as.data.frame(Quality)


print(var$contrib) # All variables Contributions to dimensions
Contrib=var$contrib 
Contrib=as.data.frame(Contrib)


# Variable graphs
NOMpng="results/AllVariables_FAMD.png"
png(file = NOMpng, width = 700, height = 700)
fviz_famd_var(res_famd, repel = T, #repel to unlabel data points with overlaps
              labelsize=5,pointsize=4)+ 
  theme(text = element_text(size = 17),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 17))
dev.off()



NOMpng="results/AllVariables2_FAMD.png"
png(file = NOMpng, width = 700, height = 700)
plot(res_famd, choix="var",cex=1.5,cex.axis=1.8,font.axis=1.5)
dev.off()


# All Variables Contributions to dimensions 
# dotted line show average value if contributions were uniform 
# axis 1
fviz_contrib(res_famd, "var", axes = 1)+ # Contribution to 1st dimension (DIM1) 
  theme(text = element_text(size = 17),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 10))

# axis 2
fviz_contrib(res_famd, "var", axes = 2)+ # Contribution to 2nd dimension (DIM2)
  theme(text = element_text(size = 17),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 10))


Top10=fviz_contrib(res_famd, choice = "var", axes = 1:2, top =10,sort.val="desc")
Top10=as.data.frame(Top10[["layers"]][[1]][["data"]][["name"]]) # extract the names of the most contributives variables
names(Top10)[1] <- "variables" # Rename first column



# Graphs of quantitative variables -------------------------------------------------------------

# Get results for quantitative variables 
quanti.var <- get_famd_var(res_famd, "quanti.var") # for quantitative variables

print(quanti.var$cos2) #cos? value for each dimension
Quality=quanti.var$cos2
Quality=as.data.frame(Quality)


print(quanti.var$contrib) # variables Contributions to dimensions
Contrib=quanti.var$contrib 
Contrib=as.data.frame(Contrib)




# Graphs
NOMpng="results/Quantitative_Variables_FAMD.png"
png(file = NOMpng, width = 700, height = 700)
fviz_famd_var(res_famd, "quanti.var", repel=T, col.var = "contrib",axes=c(1,2),gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),labelsize=4,pointsize=4)+ 
  theme(text = element_text(size = 17),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 15))
dev.off()


# Correlations
corrplot(na.omit(quanti.var$cos2), is.corr=FALSE,tl.cex=1.4,tl.col = "black",cl.cex=1.4,cl.align.text="l") # show the link between variables , representation quality and correlations with variables and dimensions


# Graphique des variables qualitatives -------------------------------------------------------------

# Get results for qualitative variables 
quali.var <- get_famd_var(res_famd, "quali.var") # for qualitative variables

print(quali.var$cos2) #cos? value for each dimension
Quality=quali.var$cos2
Quality=as.data.frame(Quality)

print(quali.var$contrib) # variables Contributions to dimensions
Contrib=quali.var$contrib 
Contrib=as.data.frame(Contrib)


#Graphs 

fviz_famd_var(res_famd, "quali.var", repel=T,col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),labelsize=5,pointsize=0.7)+ 
  theme(text = element_text(size = 17),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 15))


plot(res_famd, choix="quali",cex=1,cex.axis=1.8,font.axis=1.5)




# Individual graphs-------------------------------------------------------------


# Individual graphs cos2

fviz_famd_ind(res_famd,repel=TRUE, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),labelsize=5,pointsize=4,ggrepel = TRUE, geom = c("text","point"))+ 
  theme(text = element_text(size = 17),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 15))


plot(res_famd, choix = "ind")



# Color individuals according to morphotype , altitude and sex

L=list("Altitude","Morphotype","Sex")
Nl=length(L)
library(viridis)

for (i in seq (1:Nl)){
  
  print(paste0("Individuals are coloured according to : ",L[[i]]))
  
  NOMpng="results/Individuals_MorphoSex_FAMD.png"
  png(file = NOMpng, width = 700, height = 700)
  p=fviz_famd_ind(res_famd, geom="point",habillage=c("Morphotype","Sex"),pointsize = 4,
                  # colored by group
                  addEllipses = T, ellipse.type = "confidence",repel = TRUE,ellipse.level = 0.95)+  
    scale_color_viridis(discrete = TRUE,begin=0.3,end=0.9)+
    scale_fill_viridis(discrete = TRUE,begin=0.3,end=0.9)+
    theme(text = element_text(size = 17),axis.title = element_text(size = 17),axis.text = element_text(size = 15),legend.position = "bottom")
  
  print(p)
  dev.off()
  
  
  
  NOMpng=paste0("Individuals_",L[[i]],"_FAMD.png")
  png(file = NOMpng, width = 700, height = 700)
  p=fviz_famd_ind(res_famd, geom= "point",habillage=L[[i]],pointsize = 4,
                  # colored by group
                  addEllipses = T, ellipse.type = "confidence",repel = TRUE,ellipse.level = 0.95)+  
    scale_color_viridis(discrete = TRUE,begin=0.3,end=0.9)+
    scale_fill_viridis(discrete = TRUE,begin=0.3,end=0.9)+
    theme(text = element_text(size = 17),axis.title = element_text(size = 17),axis.text = element_text(size = 15))
  print(p)
  dev.off()
  
}


