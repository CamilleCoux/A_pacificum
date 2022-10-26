
# Name : PCA analysis
# Authors : Diane ESPEL
# Objectives :  computing body size index for the study : Espel et al. (in prep)
#----------------------------------------------------------------------------

# R version : 2.15.2

# Required packages -------------------------------------------------------

library(FactoMineR)
library(factoextra)
library(corrplot)
library(naniar)

# The body size index columns are already in the datasheet, this is to show
# how we calculated them.


# Open data  -------------------------------------------------------------

DATA <- read.csv2("data/Altitudinal transects Amblystogenium pacificum.csv")

# Check missing rows
info=skimr::skim(DATA)
print(info)
print(paste0("% missing data : ",round(pct_miss(DATA)*10^2,2),"%"))


#Remove outliers
str(DATA)
newDATA=subset(DATA,Individual!=3)
newDATA=subset(newDATA,Individual!=28)
newDATA=subset(newDATA,Individual!=40)
newDATA=subset(newDATA,Individual!=45)
tab=newDATA[,c(8:13)] # Keep only morphometrics traits (length and width)


# Apply PCA  -------------------------------------------------------------

res_pca=PCA(tab,  #dateframe with n rows and p columns
            ncp = 5, # maximum number of dimensions
            graph = T# display the graph
)



# Variances  -------------------------------------------------------------

# get eigen values
eig.val=res_pca$eig 
eig.val=as.data.frame(eig.val)
print(eig.val)


# Plot of Eigen values
NOMpng="results/Eigenvalues_pca.png"
png(file = NOMpng, width = 400, height = 400)
fviz_eig(res_pca,title="Eigenvalues")+ 
  theme(text = element_text(size = 17),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 15))
dev.off()



# Graphs of all variables  -------------------------------------------------------------

#  get results for variables
var=get_pca_var(res_pca) # for both quantitative and qualitative variables

print(var$coord) # All variables coordinates

print(var$cos2) # Quality of the representation
Quality=var$cos2
Quality=as.data.frame(Quality)


print(var$contrib) # All variables Contributions to dimensions
Contrib=var$contrib 
Contrib=as.data.frame(Contrib)



# Variable graphs
NOMpng="results/AllVariables_pca.png"
png(file = NOMpng, width = 700, height = 700)
fviz_pca_var(res_pca, repel=T, col.var = "contrib",axes=c(1,2),gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),labelsize=4,pointsize=4)+ 
  theme(text = element_text(size = 17),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 15))
dev.off()

plot(res_pca,cex=1,cex.axis=1.8,font.axis=1.5)



# All Variables Contributions to dimensions 
# 1st axis:
fviz_contrib(res_pca, "var", axes = 1)+ # Contribution to 1st dimension (DIM1) 
  theme(text = element_text(size = 17),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 10))


# 2nd axis : 
fviz_contrib(res_pca, "var", axes = 2)+ # Contribution to 2nd dimension (DIM2)
  theme(text = element_text(size = 17),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 10))


Top10=fviz_contrib(res_pca, choice = "var", axes = 1:2, top =10,sort.val="desc")
Top10=as.data.frame(Top10[["layers"]][[1]][["data"]][["name"]]) # extract the names of the most contributives variables
names(Top10)[1] <- "variables" # Rename first column
# FILE8=paste0(save_path,District,"_",Island,"_10_most_contributives_variables.csv")
# write.csv(Top10,FILE8)


# Computing correlations -------------------------------------------------------------

corrplot(na.omit(var$cos2), is.corr=FALSE,tl.cex=1.4,tl.col = "black",cl.cex=1.4,cl.align.text="l") # show the link between variables , representation quality and correlations with variables and dimensions



# Individual graphs-------------------------------------------------------------

# Individual graphs cos2

fviz_pca_ind(res_pca, repel=TRUE, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),labelsize=5,pointsize=4,ggrepel = TRUE, geom = c("text","point"))+ 
  theme(text = element_text(size = 17),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 15))



# getting coordinates-------------------------------------------------------------

coordinates=as.data.frame(res_pca[["ind"]][["coord"]])
plot(coordinates$Dim.1,newDATA$Body_size_index)

newDATA$Body_size_index_corr=coordinates$Dim.1


# Save new data  -------------------------------------------------------------
library(xlsx)
save_file="data/Altitudinal transects Amblystogenium pacificum_new.csv"
write.csv2(newDATA, save_file)

