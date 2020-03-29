#Principal component analysis 
#rather use multiple factor analysis as there are both factorial and continious variables
#ressources used:
#https://www.youtube.com/watch?v=g5_hM93e8HM (as suggested by morane)
#http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/116-mfa-multiple-factor-analysis-in-r-essentials/#r-code
#by Dana Looschelders
#created: 28.03.2020
#changelog: 

library(ggplot2)
library(car)
library(vegan)
library(FactoInvestigate)
library(FactoMineR)
library(factoextra)
library(dplyr)
setwd("C:/00 Dana/Uni/Internship/Work")
data=read.table("data_for_pca.csv", sep=";", dec=",", header=T)
str(data)

#prepare the data for multiple factor analysis (mfa)
#set all values <0 for senescence as NA
data$Phyllochron1S[data$Phyllochron1S<0]=NA
str(data) #check data
data=data[,-4] #remove column with chamber 
#normalize data-> get all values between 0 and 1
for(x in 5:14){
  data[,x]=(data[,x]-min(data[,x], na.rm=T))/(max(data[,x], na.rm=T)-min(data[,x], na.rm=T))
}
str(data)  
#transpose?
#data=as.data.frame(t(data))

#set plant ID as rownames
rownames(data)=data$ï..Plant_ID
data=data[,-1] #drop column with ID as it is now the rowname
data=data[,c(13,1:12)] #rearrange type and variety together as the first two columns
data=data[,-4] #drop column with first biomass values
str(data)
?MFA
#base: data
#groups:
#1. group: kind (2 columns) Variety and Type (categorical) "n"
#2. group: Treatment (1 column) (categorical) "n"
#3. group: phenotype (2 columns) Biomass and Height (continious) "c"
#**** include both biomass samples? rn only the second one is included*****
#4. group: phyllochron: (7 columns) all phyllo and senescence (continous) "c"
#type:
#c / s : quantitative variables
#s: scaled to unit variance
#n: categorical variables
#f: frequencies
#use tab.comp for missing values?
result=MFA(data, group=c(2,1,2,7), #define how many columns belong to each group
           name.group=c("kind", "treatment", "phenotype", "phyllochron"), #set names of groups
           type=c("n", "n", "c", "c"),
           graph=TRUE) 
#try with more groups (differ height and biomass as well as variety and type)
result_2=MFA(data, group=c(1,1,1,1,1,7), #define how many columns belong to each group
             name.group=c("type", "variety", "treatment", "biomass", "height", "phyllochron"), #set names of groups
             type=c("n","n", "n","c", "c", "c"),
             graph=TRUE)

#try with each variable as one group
result_3=MFA(data, group=c(1,1,1,1,1,1,1,1,1,1,1,1),#define how many columns belong to each group
             name.group=colnames(data), #set names of groups
             type=c("n","n", "n","c", "c", "c", "c", "c", "c", "c", "c", "c"),
             graph=TRUE) 

print(result)

#analyse results

#extract eigenvalues/variances retained by each dimension
eig_value=get_eigenvalue(result_3)
#visualise eigenvalues
fviz_eig(result_3)

#results for individuals
res.indi=get_mfa_ind(result_3)
fviz_mfa_ind(result_3) #visualise results as plot

#results for quantitative and qualitative variables
res.var=get_mfa_var(result_3, "group")
fviz_mfa_var(result_3) #visualse results as plot
quanti.var=get_mfa_var(result_3, "quanti.var") #results for quantitative variables
quanti.var$contrib #display contributions

res.var$correlation #display correlation between groups and principal dimensions
fviz_mfa_var(result_3, palette="jco", repel=T) #plot correlation of quantitative variables and dimensions

res.var$contrib #display contributions to dimension
fviz_contrib(result_3, "group", axes=1) #barplot with contribution to first dimension (change values from "axes" to display other dimensions)

res.var$cos2 #display quality of representation on the factor map