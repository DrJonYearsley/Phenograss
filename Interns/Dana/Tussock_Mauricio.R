### Mauricio Tussock Diameter 09 July 2020 
# Update Chambers Effect on 11 July 2020
# Update for descriptive statistics on 13 July 2020
# Update for linear model and first statistics on 14 July 2020
# Update for corrections of small typos/mistakes/legends and Division of Cultivars, Semi-naturals and Wilds on 20 July 2020

#Remove older files and clear WD
rm(list=ls())

#Set Working Directory
setwd('C:/Users/Mauricio Mantoani/Desktop/Postdoc PhenoGrass/Data Google Drive')

setwd('E:') ### Use this for the working station

#load necessary librarys to read new data in R
library(readxl) #library to read .xlsx files (otherwise .csv files are needed)
library(ggplot2) #library for graphing data
library(tidyverse)#library to organise data

#other packages and commands that are useful
library(lme4)
library(nlme)
library(arm)
library(car)
library(MASS)
library(psych)
library(agricolae)
library(PMCMRplus)
library(PMCMR)
library(emmeans)
options(max.print=1000000)

# download data and replace spaces in names with '.'
Tussock=read_excel("Tussock Diameter Data.xlsx", .name_repair = 'universal')
str(Tussock) #check if data was read in correctly

# Visualise the data
format = theme(axis.title = element_text(size=14),
            axis.text = element_text(size=12))


#Subset for Ambient and eCO2+2C
Tussock_ambient=Tussock[Tussock$Treatment=="Ambient",] #Subset to Ambient treatment only

Tussock_eCO2=Tussock[Tussock$Treatment=="eCO2 +2°C",] #Subset to eCO2 treatment only


#Boxplot with Treatments
ggplot(data=Tussock, 
       aes(x=Variety, y=Tussock.Diameter, fill=Treatment))+ #specifies the data used for x and y
  geom_boxplot(width=0.5)+ #adds tiny boxplot in violin 
  scale_fill_brewer(name='Treatment:', palette='Set2')+
  xlab("Variety")+ #adds title for x axis
  ylab("Tussock Diameter (cm)")+ #adds title for y axis
  ggtitle("Tussock Diameter (cm)")+ #adds plot title
  theme_bw()+
  format+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #rotates the axix labels on x axis


#Boxplot with Treatments and Genetic Status
ggplot(data=Tussock, 
       aes(x=Treatment, y=Tussock.Diameter, fill=Genetic))+ #specifies the data used for x and y
  geom_boxplot(width=0.5)+ #adds tiny boxplot in violin 
  scale_fill_brewer(name='Genetic Status:', palette='Set2')+
  xlab("Treatment")+ #adds title for x axis
  ylab("Tussock Diameter (cm)")+ #adds title for y axis
  ggtitle("Tussock Diameter (cm) for all Varieties by Genetic Status")+ #adds plot title
  theme_bw()+
  format


#Violin plot with Genetic Status
ggplot(data=Tussock, 
       aes(x=Genetic, y=Tussock.Diameter, fill=Treatment))+ #specifies the data used for x and y
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+ #adds violin shape without trimmed edges
  scale_fill_brewer(name='Treatment', palette='Set2')+
  xlab("Genetic Status")+ #adds title for x axis
  ylab("Tussock Diameter (cm)")+ #adds title for y axis
  ggtitle("Tussock Diameter (cm) for all Varieties by Genetic Status")+ #adds plot title
  theme_bw()+
  format


#Violin plot by Variety with boxplot 
dodge <- position_dodge(width = 0.9) ###OBS: Manual input to make boxplot inside the violin
ggplot(data=Tussock,
       aes(x=Variety, y=Tussock.Diameter, fill=Treatment))+ #specifies the data used for x and y
  geom_violin(trim=F)+ #adds violin shape without trimmed edges
  geom_boxplot(width=0.1, inherit.aes=TRUE, position=dodge)+ #adds boxplot to violin
  scale_fill_brewer(name='Treatment:', palette='Set2')+ # legend and colour of violins/boxplots
  xlab("Variety")+ #adds title for x axis
  ylab("Tussock Diameter (cm)")+ #adds title for y axis
  ggtitle("Tussock Diameter (cm) per Variety")+ #adds plot title
  theme_bw()+
  format+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #rotates the axix labels on x axis


#Violin plot by Variety with boxplot displaying Genetic Status in different colours
ggplot(data=Tussock,
       aes(x=Variety, y=Tussock.Diameter, fill=Genetic))+ #specifies the data used for x and y
  geom_violin(trim=F)+ #adds violin shape without trimmed edges
  geom_boxplot(width=0.1)+ #adds boxplot to violin
  scale_fill_brewer(name='Genetic Status:', palette='Set2')+ # legend and colour of violins/boxplots
  xlab("Variety")+ #adds title for x axis
  ylab("Tussock Diameter (cm)")+ #adds title for y axis
  ggtitle("Tussock Diameter (cm) per Variety")+ #adds plot title
  theme_bw()+
  format+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #rotates the axix labels on x axis


#Histogram plot with Frequency
hist(Tussock$Tussock.Diameter, #with frequency
     main = "Histogram of Tussock Diameter (cm) (Frequency)") #plot title

#Histogram plot with Density
hist(Tussock$Tussock.Diameter, prob=TRUE, 
     main = 'Histogram of Tussock  Diameter (cm) (Density)', col='green') #histogram 
lines(density(Tussock$Tussock.Diameter), col="black") #density as red line

####### Plots for specific subset of Ambient only


#Box plot with Treatment Ambient Only
ggplot(data=Tussock_ambient, 
       aes(x=Variety, y=Tussock.Diameter, fill=Treatment))+ #specifies the data used for x and y
  geom_boxplot(width=0.2)+ #adds tiny boxplot in violin 
  scale_fill_brewer(name='Treatment:', palette='Set2')+
  xlab("Variety")+ #adds title for x axis
  ylab("Tussock Diameter (cm) Ambient")+ #adds title for y axis
  ggtitle("Tussock Diameter (cm) for Ambient Treatment")+ #adds plot title
  theme_bw()+
  format+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #rotates the axix labels on x axis


#Violin plot with Genetic Status for Ambient Only
ggplot(data=Tussock_ambient, 
       aes(x=Genetic, y=Tussock.Diameter, fill=Treatment))+ #specifies the data used for x and y
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+ #adds violin shape without trimmed edges
  scale_fill_brewer(name='Treatment', palette='Set2')+
  xlab("Genetic Status")+ #adds title for x axis
  ylab("Tussock Diameter (cm) Ambient")+ #adds title for y axis
  ggtitle("Tussock Diameter (cm) for Ambient all Varieties")+ #adds plot title
  theme_bw()+
  format


#Violin plot by Variety with boxplot for Ambient only
dodge <- position_dodge(width = 0.9) ###OBS: Manual input to make boxplot inside the violin
ggplot(data=Tussock_ambient,
       aes(x=Variety, y=Tussock.Diameter, fill=Genetic))+ #specifies the data used for x and y
  geom_violin(trim=F)+ #adds violin shape without trimmed edges
  geom_boxplot(width=0.1, inherit.aes=TRUE, position=dodge)+ #adds boxplot to violin
  scale_fill_brewer(name='Treatment:', palette='Set2')+ # legend and colour of violins/boxplots
  xlab("Variety")+ #adds title for x axis
  ylab("Tussock Diameter (cm) Ambient")+ #adds title for y axis
  ggtitle("Tussock Diameter (cm) Ambient")+ #adds plot title
  theme_bw()+
  format+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #rotates the axix labels on x axis


#Histogram plot with Frequency for Ambient Only
hist(Tussock_ambient$Tussock.Diameter, #with frequency
     main = "Histogram of Tussock Diameter (cm) (Frequency) Ambient  all Varieties") #plot title

#Histogram plot with Density
hist(Tussock_ambient$Tussock.Diameter, prob=TRUE, 
     main = 'Histogram of Tussock Diameter (cm) (Density) Ambient  all Varieties', col='green') #histogram 
lines(density(Tussock_ambient$Tussock.Diameter), col="black") #density as red line


####### Plots for specific subset of eCO2 only


#Box plot with Treatment eCO2 Only
ggplot(data=Tussock_eCO2, 
       aes(x=Variety, y=Tussock.Diameter, fill=Treatment))+ #specifies the data used for x and y
  geom_boxplot(width=0.2)+ #adds tiny boxplot in violin 
  scale_fill_brewer(name='Treatment:', palette='Set1')+
  xlab("Variety")+ #adds title for x axis
  ylab(" Tussock Diameter (cm) eCO2")+ #adds title for y axis
  ggtitle("Tussock Diameter (cm) for eCO2 Treatment")+ #adds plot title
  theme_bw()+
  format+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #rotates the axix labels on x axis


#Violin plot with Genetic Status for eCO2 Only
ggplot(data=Tussock_eCO2, 
       aes(x=Genetic, y=Tussock.Diameter, fill=Treatment))+ #specifies the data used for x and y
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+ #adds violin shape without trimmed edges
  scale_fill_brewer(name='Treatment:', palette='Set1')+
  xlab("Genetic Status")+ #adds title for x axis
  ylab("Tussock Diameter (cm) eCO2")+ #adds title for y axis
  ggtitle("Tussock Diameter (cm) for eCO2")+ #adds plot title
  theme_bw()+
  format


#Violin plot by Variety with boxplot for eCO2 only
dodge <- position_dodge(width = 0.9) ###OBS: Manual input to make boxplot inside the violin
ggplot(data=Tussock_eCO2,
       aes(x=Variety, y=Tussock.Diameter, fill=Genetic))+ #specifies the data used for x and y
  geom_violin(trim=F)+ #adds violin shape without trimmed edges
  geom_boxplot(width=0.1, inherit.aes=TRUE, position=dodge)+ #adds boxplot to violin
  scale_fill_brewer(name='Treatment:', palette='Set2')+ # legend and colour of violins/boxplots
  xlab("Variety")+ #adds title for x axis
  ylab("Tussock Diameter (cm) eCO2")+ #adds title for y axis
  ggtitle("Tussock Diameter (cm) eCO2")+ #adds plot title
  theme_bw()+
  format+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #rotates the axix labels on x axis


#Histogram plot with Frequency for eCO2 Only
hist(Tussock_eCO2$Tussock.Diameter, #with frequency
     main = "Histogram of Tussock Diameter (cm) (Frequency) eCO2  all Varieties") #plot title

#Histogram plot with Density
hist(Tussock_eCO2$Tussock.Diameter, prob=TRUE, 
     main = 'Histogram of Tussock Diameter (cm) (Density) eCO2  all Varieties', col='red') #histogram 
lines(density(Tussock_eCO2$Tussock.Diameter), col="black") #density as red line


###=======### Update on 11 July 2020 and Division of Cultivars, Semi-natural and Wilds on 20 July 2020

####### Checking for Chamber effect on Tiller Number 


#To factor chamber for the graphs
Tussock$Chamber=factor(Tussock$Chamber, 
                          levels=unique(Tussock$Chamber))

#Violin plot of Tussock with Chamber division
dodge <- position_dodge(width = 0.9) ###OBS: Manual input to make boxplot inside the violin
ggplot(data=Tussock, 
       aes(x=Treatment, y=Tussock.Diameter, fill=Chamber)) +
  geom_violin() + 
  geom_boxplot(width=0.1, inherit.aes=TRUE, position=dodge)+ #adds boxplot to violin
  scale_fill_brewer(name='Chamber:', 
                    palette='Set2', 
                    breaks=c(1,2,7,8),
                    labels=c('1: Ambient', 
                             '2: eCO2', 
                             '7: Ambient',
                             '8: eCO2')) +
  xlab("Treatment")+ #adds title for x axis
  ylab("Tussock Diameter (cm)")+ #adds title for y axis
  ggtitle("Tussock Diameter (cm) by Chambers for all Varieties")+ #adds plot title
  theme_bw() + 
  format

#Subset to Cultivars only
Tussock_Cultivar=Tussock[Tussock$Genetic=="Cultivar",] #Subset to Cultivars only

Tussock_Cultivar$Chamber=factor(Tussock_Cultivar$Chamber, 
                                   levels=unique(Tussock_Cultivar$Chamber))

#Boxplot of Tussock with Chamber division for CULTIVARS ONLY
ggplot(data=Tussock_Cultivar, 
       aes(x=Chamber, y=Tussock.Diameter, fill=Variety)) +
  geom_boxplot(width=0.3, inherit.aes=TRUE, position=dodge)+ #adds boxplot to violin
  scale_fill_brewer(name='Variety:', palette='Set2')+
  xlab("Chamber")+ #adds title for x axis
  ylab("Tussock Diameter (cm)")+ #adds title for y axis
  ggtitle("Tussock Diameter (cm) by Chambers for Cultivars only")+ #adds plot title
  theme_bw() + 
  format

#Boxplot of Tussock with Chamber division for CULTIVARS only
ggplot(data=Tussock_Cultivar, 
       aes(x=Variety, y=Tussock.Diameter, fill=Chamber)) +
  geom_boxplot(width=0.3, inherit.aes=TRUE, position=dodge)+ #adds boxplot to violin
  scale_fill_brewer(name='Chamber:', 
                    palette='Set2', 
                    breaks=c(1,2,7,8),
                    labels=c('1: Ambient', 
                             '2: eCO2', 
                             '7: Ambient',
                             '8: eCO2')) +
  xlab("Variety")+ #adds title for x axis
  ylab("Tussock Diameter (cm)")+ #adds title for y axis
  ggtitle("Tussock Diameter (cm) by Chambers for Cultivars only")+ #adds plot title
  theme_bw() + 
  format

#Subset to WILD only
Tussock_Wild=Tussock[Tussock$Genetic=="Wild",] #Subset to Wild only

Tussock_Wild$Chamber=factor(Tussock_Wild$Chamber, 
                               levels=unique(Tussock_Wild$Chamber))

#Violin plot of Tussock with Chamber division for WILD ONLY
ggplot(data=Tussock_Wild, 
       aes(x=Chamber, y=Tussock.Diameter, fill=Variety)) +
  geom_boxplot(width=0.3, inherit.aes=TRUE, position=dodge)+ #adds boxplot to violin
  scale_fill_brewer(name='Variety:', palette='Set2')+
  xlab("Chamber")+ #adds title for x axis
  ylab("Tussock Diameter (cm)")+ #adds title for y axis
  ggtitle("Tussock Diameter (cm) by Chambers for Wild only")+ #adds plot title
  theme_bw() + 
  format

#Violin plot of Tussock with Chamber division for WILD only
ggplot(data=Tussock_Wild, 
       aes(x=Variety, y=Tussock.Diameter, fill=Chamber)) +
  geom_boxplot(width=0.3, inherit.aes=TRUE, position=dodge)+ #adds boxplot to violin
  scale_fill_brewer(name='Chamber:', 
                    palette='Set2', 
                    breaks=c(1,2,7,8),
                    labels=c('1: Ambient', 
                             '2: eCO2', 
                             '7: Ambient',
                             '8: eCO2')) +
  xlab("Variety")+ #adds title for x axis
  ylab("Tussock Diameter (cm)")+ #adds title for y axis
  ggtitle("Tussock Diameter (cm) by Chambers for Wild only")+ #adds plot title
  theme_bw() + 
  format

#Subset to Semi-natural only
Tussock_Semi=Tussock[Tussock$Genetic=="Semi-natural",] #Subset to Semi-natural only

Tussock_Semi$Chamber=factor(Tussock_Semi$Chamber, 
                            levels=unique(Tussock_Semi$Chamber))

#Violin plot of Tussock with Chamber division for Semi-natural ONLY
ggplot(data=Tussock_Semi, 
       aes(x=Chamber, y=Tussock.Diameter, fill=Variety)) +
  geom_boxplot(width=0.3, inherit.aes=TRUE, position=dodge)+ #adds boxplot to violin
  scale_fill_brewer(name='Variety:', palette='Set2')+
  xlab("Chamber")+ #adds title for x axis
  ylab("Tussock Diameter (cm)")+ #adds title for y axis
  ggtitle("Tussock Diameter (cm) by Chambers for Semi-natural only")+ #adds plot title
  theme_bw() + 
  format

#Violin plot of Tussock with Chamber division for Semi-natural only
ggplot(data=Tussock_Semi, 
       aes(x=Variety, y=Tussock.Diameter, fill=Chamber)) +
  geom_boxplot(width=0.3, inherit.aes=TRUE, position=dodge)+ #adds boxplot to violin
  scale_fill_brewer(name='Chamber:', 
                    palette='Set2', 
                    breaks=c(1,2,7,8),
                    labels=c('1: Ambient', 
                             '2: eCO2', 
                             '7: Ambient',
                             '8: eCO2')) +
  xlab("Variety")+ #adds title for x axis
  ylab("Tussock Diameter (cm)")+ #adds title for y axis
  ggtitle("Tussock Diameter (cm) by Chambers for Semi-natural only")+ #adds plot title
  theme_bw() + 
  format


###=======### Update on 13 July 2020 

# Summary to get all relevant statistical parameters
#Summary to get all relevant statistical parameters
summary(Tussock)
summary(Tussock_ambient)
summary(Tussock_eCO2)

# Or using package Psych to descrive statistics by specific group
describe.by(Tussock_eCO2, group="Variety")
describe.by(Tussock_ambient, group="Variety")

###=======### Update on 14 July 2020

# Construct Linear Models for different data for Tussock Diameter
m0 = lm(formula = Tussock.Diameter ~ Variety, data = Tussock) 
m1 = lm(formula = Tussock.Diameter ~ Treatment, data = Tussock)
m2 = lm(formula = Tussock.Diameter ~ Variety + Treatment, data = Tussock)
m3 = lm(formula = Tussock.Diameter ~ Variety + Treatment + Variety:Treatment, data = Tussock)

# Check how model is structurated  for Tussock Diameter
summary(m0)
summary(m1)
summary(m2)
summary(m3)

# Model Validation and qq plots and others to check for residuals distribution and normality  for Tussock Diameter
par(mfrow = c(3, 2)) # 3 x 2 pictures at a time
plot(m0, which=c(1:5)) # plot 5 model validation graphs
hist(m0$residuals,20)
par(mfrow=c(1, 1)) # Return to default 1x1 plot

par(mfrow = c(3, 2)) # 3 x 2 pictures at a time
plot(m1, which=c(1:5)) # plot 5 model validation graphs
hist(m1$residuals,20)
par(mfrow=c(1, 1)) # Return to default 1x1 plot

par(mfrow = c(3, 2)) # 3 x 2 pictures at a time
plot(m2, which=c(1:5)) # plot 5 model validation graphs
hist(m2$residuals,20)
par(mfrow=c(1, 1)) # Return to default 1x1 plot

par(mfrow = c(3, 2)) # 3 x 2 pictures at a time
plot(m3, which=c(1:5)) # plot 5 model validation graphs
hist(m3$residuals,20)
par(mfrow=c(1, 1)) # Return to default 1x1 plot

# Anova for the models  for Tussock Diameter
anova(m0) # Anova acounting for Variety only
anova(m1) # Anova acounting for Treatment only
anova(m2) # Anova acounting for Variety and Treatment
anova(m3) # Anova acounting for Variety and Treatment and their interaction

# Posthoc test to find differences between Varieties and Treatments since there is no interaction between Variety:Treatment  for Tussock Diameter
emmeans(m0, pairwise ~ Variety, adjust="tukey")
emmeans(m1, pairwise ~ Treatment, adjust="tukey")
emmeans(m2, pairwise ~ Variety+Treatment, adjust="tukey") #which one should we refer to?
emmeans(m3, pairwise ~ Variety:Treatment, adjust="tukey") #which one should we refer to?


###=======###