#----------------------------------------------------------------------------
# Name : Mixed model analysis 
# Authors : Camille Coux
# Objectives : fitting mixed models with each functional trait as a response 
# variable, and altitude, morphotype and sex as covariables. Model selection 
# is applied to each model
#----------------------------------------------------------------------------

# Required packages -------------------------------------------------------

library(magrittr)
library(tidyverse)
library(nlme)
library(cowplot)
library(MuMIn)
library(viridis)
library(gridExtra)

# Theme preferences for ggplots
theme_set(theme_minimal())

# read in data 

data <- read.csv2("data/Altitudinal transects Amblystogenium pacificum.csv")
data <- data[,-20]
data <- data[-c(3, 27, 44, 39),] # outliers
data$Sex %<>% as.factor
data$Site %<>% as.factor

# rename levels 
data$Altitude <- gsub("m", "", data$Altitude) %>% as.numeric
data$Morphotype <- gsub("N", "Dark", data$Morphotype)
data$Morphotype <- gsub("B", "Light", data$Morphotype)
data$Morphotype %<>% as.factor

# remove correlated traits:
data <- data[,-grep("_l|_L|Proteines$|Lipids$|Interoc", colnames(data))]

#check normality visually
par(mfrow=c(2, 2))
apply(data[,7:11], 2, hist) # lipids and sugars are not normally distributed.
par(mfrow=c(1,1))


# try log transformation : 
data$log.Lipids.Mass <- data$Lipids.Mass %>% log
hist(data$log.Lipids.Mass) # better
data$log.Sugars.Mass <- data$Sugars.Mass %>% log
hist(data$log.Sugars.Mass) # better

# remove the previous versions
#data <- data[, -grep("^Lipids.Mass|^Sugars.Mass$", colnames(data))]

# 
# # make a joint variable for morphotype and sex to facilitate ggplotting
# data$MorphoSex<-paste(data$Morphotype,data$Sex,sep="_")
# 
# # put everyone back to a factor
# data$Site %<>% as.factor
# data$Morphotype %<>% as.factor
# data$Sex %<>% as.factor
# data$MorphoSex %<>% as.factor


#--- Area plot of relative abundances across sites
# morphortypes as proportions

m <- data %>%
  group_by(Altitude) %>%
  count(Morphotype) %>%
  mutate(total = sum(n),
         Proportions = n/total) %>%
  ggplot()+
  geom_area(aes(x= Altitude, y=Proportions, fill=Morphotype), alpha = 0.4, position = 'identity') +
  ylim(0, 1) +
  #scale_fill_viridis(discrete = TRUE,  begin = 0.9, end = 0.3, alpha = 0.6)
  scale_fill_manual(values=c("#69b3a2", "#D8AE5A") )+
  geom_line(aes(x=Altitude, y=Proportions, color = Morphotype), size=1 )+
  scale_color_manual(values=c("#69b3a2", "#D8AE5A") ) +
  theme_minimal() 


# sexes as proportions
s <- data %>%
  group_by(Altitude) %>%
  count(Sex) %>%
  mutate(total = sum(n),
         Proportions = n/total) %>%
  ggplot()+
  geom_area(aes(x= Altitude, y=Proportions, fill=Sex), alpha = 0.4, position = 'identity') +
  ylim(0, 1) +
  #scale_fill_viridis(discrete = TRUE,  begin = 0.9, end = 0.3, alpha = 0.6)
  scale_fill_manual(values=c("#69b3a2", "#D8AE5A") )+
  geom_line(aes(x=Altitude, y=Proportions, color = Sex), size=1 )+
  scale_color_manual(values=c("#69b3a2", "#D8AE5A") )


# Linear Mixed Effect models for each of the traits that we selected ------------------

# MASS
# full model:
mass <- lme(Mass ~ I(scale(Altitude)^2) + scale(Altitude) + Sex + Morphotype, 
            random = ~ 1|Site, 
            data = data, method ="REML") 

summary(mass)
# nothing significant

# model selection : 
# re-specify the model and set method = ML (model selection based on AIC isn't reliable
# for REstricted Max Likelihood method)
mass.mod <- lme(Mass ~ I(scale(Altitude)^2) + scale(Altitude) * Sex * Morphotype, 
                random = ~ 1|Site, 
                data = data, method ="ML") 

# make list of all model combinations : 
mass.msl <- dredge(mass.mod) # dredge function 

# average over the set of best fitting models that have a difference in AIC score < 2 
summary(model.avg(mass.msl, subset = delta < 2))
# nothing significant



# PROTEINS
# full model : 
proteins <- lme(Proteines.Mass ~ I(scale(Altitude)^2) + scale(Altitude) + Sex + Morphotype, 
                random = ~ 1|Site, 
                data = data, method ="REML") 
summary(proteins)


# make list of all model combinations : 
prot.mod <- lme(Proteines.Mass ~ I(scale(Altitude)^2) + scale(Altitude) * Sex * Morphotype, 
                random = ~ 1|Site, 
                data = data, method ="ML") 

# model averaging : 
prot.msl <- dredge(prot.mod)
summary(model.avg(prot.msl, subset = delta < 2))



# LIPIDS
# full model : 
lipids <- lme(log.Lipids.Mass ~ I(scale(Altitude)^2) + scale(Altitude) + Sex + Morphotype, 
              random = ~ 1|Site, 
              data = data, method ="REML") 
summary(lipids)
# nothing significant

# make list of all model combinations : 
lip.mod <- lme(log.Lipids.Mass ~ I(scale(Altitude)^2) + scale(Altitude) * Sex * Morphotype, 
               random = ~ 1|Site, 
               data = data, method ="ML") 

# model averaging : 
lip.msl <- dredge(lip.mod)

# nothing significant here
summary(model.avg(lip.msl, subset = delta < 2)) 





# SUGARS
# full model :
sugar <- lme(log.Sugars.Mass ~ I(scale(Altitude)^2) + scale(Altitude) + Sex + Morphotype, 
             random = ~ 1|Site, 
             data = data, method ="REML") 
summary(sugar)

# make list of all model combinations : 
sugar.mod <- lme(log.Sugars.Mass ~ I(scale(Altitude)^2) + scale(Altitude) * Sex * Morphotype, 
                 random = ~ 1|Site, 
                 data = data, method ="ML") 
summary(sugar.mod)

# model averageing : 
sugar.msl <- dredge(sugar.mod)
summary(model.avg(sugar.msl, subset = delta < 2)) 



# BODY SIZE INDEX
# full model : 
BSI <- lme(Body_size_index ~ I(scale(Altitude)^2) + scale(Altitude) + Sex + Morphotype, 
           random = ~ 1|Site, 
           data = data, method ="REML") 
summary(BSI)

# make list of all model combinations : 
BSI.mod <- lme(Body_size_index ~ I(scale(Altitude)^2) + scale(Altitude) * Sex * Morphotype, 
               random = ~ 1|Site, 
               data = data, method ="ML") 
summary(BSI.mod)

# model averageing : 
BSI.msl <- dredge(BSI.mod)
summary(model.avg(BSI.msl, subset = delta < 2)) 


# plot all of the traits as a function of Altitude, morphotype and sex
p0 <- ggplot(data, aes(Altitude, Mass, color = Morphotype, linetype = Sex)) +
  geom_point(size = 1) +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1,alpha=0.1) +
  scale_color_viridis(discrete = TRUE, begin = 0.3, end = 0.9) 


p1<-ggplot(data, aes(Altitude, Mass, color = Morphotype, linetype = Sex)) +
  geom_point(size = 1) +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1,alpha=0.1) +
  scale_color_viridis(discrete = TRUE, begin = 0.3, end = 0.9) +
  theme(legend.position="none")

p1<-ggplot(data, aes(Altitude, Mass, color = Morphotype, linetype = Sex)) +
  geom_point(size = 1) +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1,alpha=0.1) +
  scale_color_viridis(discrete = TRUE, begin = 0.3, end = 0.9) +
  theme(legend.position="none")

p2<-ggplot(data, aes(Altitude, Proteines.Mass,color = Morphotype, linetype = Sex)) +
  geom_point(size = 1) +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1,alpha=0.1)+
  scale_color_viridis(discrete = TRUE, begin = 0.3, end = 0.9) +
  labs(y="Proteins") +
  theme(legend.position="none")

p3<-ggplot(data, aes(Altitude, log.Lipids.Mass, color = Morphotype, linetype = Sex)) +
  geom_point(size = 1) +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1,alpha=0.1)+
  scale_color_viridis(discrete = TRUE, begin = 0.3, end = 0.9) +
  labs(y="log(Lipids)") +
  theme(legend.position="none")

p4<-ggplot(data, aes(Altitude, log.Sugars.Mass, color = Morphotype, linetype = Sex)) +
  geom_point(size = 1) +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1,alpha=0.1)+
  scale_color_viridis(discrete = TRUE, begin = 0.3, end = 0.9) +
  labs(y="log(Sugars)") +
  theme(legend.position="none")

p5<-ggplot(data, aes(Altitude, Body_size_index, color = Morphotype, linetype = Sex)) +
  geom_point(size = 1) +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1,alpha=0.1)+
  scale_color_viridis(discrete = TRUE, begin = 0.3, end = 0.9) +
  labs(y="Body Size Index") +
  theme(legend.position="none")

leg <- get_legend(p0) 

grid<-grid.arrange(p5,p1,p2,p4, p3, leg, ncol=2)


ggsave("results/LMEs.png",grid,width=6,height=8)




