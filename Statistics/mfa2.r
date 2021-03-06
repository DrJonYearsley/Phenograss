library(ggplot2)
library(FactoInvestigate)
library(FactoMineR)
library(factoextra)
library(dplyr)
library(missMDA)
library(corrplot)
library(RColorBrewer)
setwd("C:/00 Dana/Uni/Internship/Work")
data=read.table("data_for_pca2.csv", sep=";", dec=",", header=T)
str(data)

#prepare the data for multiple factor analysis (mfa)

str(data) #check data

#set plant ID as rownames
rownames(data)=data$ï..Plant_ID
data=data[,-1] #drop column with ID as it is now the rowname

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
fviz_mfa_var(result_3, "quanti.var", col.var = "contrib", 
             gradient.cols = c("red", "green", "blue"), 
             repel = TRUE,
             geom = c("point", "text")) #plot with colors in relation to contribution

# Color by cos2 values: quality on the factor map
fviz_mfa_var(result_3, col.var = "cos2",
             gradient.cols = c("red", "green", "blue"), 
             col.var.sup = "violet", repel = TRUE) #plot with colors in relation of quality of representation

fviz_cos2(result_3, choice = "quanti.var", axes = 1) #plot quality of contribution as bar plot
res.var$contrib #display contributions to dimension
fviz_contrib(result_3, "group", axes=1) #barplot with contribution to first dimension (change values from "axes" to display other dimensions)

fviz_mfa_ind(result_3, col.ind = "cos2", 
             gradient.cols = c("red", "green", "blue"),
             repel = TRUE) #plot individuals

fviz_mfa_ind(result_3, 
             habillage = "Treatment", 
             palette = c("red", "green", "blue", "brown"),
             addEllipses = TRUE, ellipse.type = "confidence", 
             repel = TRUE) # plot individuals with colour by treatment

fviz_mfa_ind(result_3, 
             habillage = "Type", 
             palette = c("red", "green", "blue"),
             addEllipses = TRUE, ellipse.type = "confidence", 
             repel = TRUE) # plot individuals with colour by type

fviz_mfa_ind(result_3, 
             habillage = "Variety", 
             palette = c("red", "green", "blue", 
                         "brown", "pink", "yellow", 
                         "black", "orange", "lightblue", 
                         "darkblue", "darkgreen", "lightgreen",
                         "turquoise", "grey"),
             addEllipses = TRUE, ellipse.type = "confidence", 
             repel = TRUE) # plot individuals with colour by variety

res.var$cos2 #display quality of representation on the factor map

fviz_ellipses(result_3, c("Treatment", "Variety"), repel = TRUE) #plot Treatment and Variety together as Eliipses

fviz_mfa_axes(result_3, repel=T) #plot relationship between principal axes and results

library(ggplot2)
library(FactoInvestigate)
library(FactoMineR)
library(factoextra)
library(dplyr)
library(missMDA)
library(corrplot)
library(RColorBrewer)
setwd("C:/00 Dana/Uni/Internship/Work")
data=read.table("data_for_mfa.csv", sep=";", dec=",", header=T)
str(data)

#prepare the data for multiple factor analysis (mfa)
#set all values <0 for senescence as NA
data$Phyllochron1S[data$Phyllochron1S<0]=NA
str(data) #check data
data=data[,-4] #remove column with chamber 
#normalize data-> get all values between 0 and 1
for(x in 5:12){
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
fviz_mfa_var(result_3, "quanti.var", col.var = "contrib", 
             gradient.cols = c("red", "green", "blue"), 
             repel = TRUE,
             geom = c("point", "text")) #plot with colors in relation to contribution

# Color by cos2 values: quality on the factor map
fviz_mfa_var(result_3, col.var = "cos2",
             gradient.cols = c("red", "green", "blue"), 
             col.var.sup = "violet", repel = TRUE) #plot with colors in relation of quality of representation

fviz_cos2(result_3, choice = "quanti.var", axes = 1) #plot quality of contribution as bar plot
res.var$contrib #display contributions to dimension
fviz_contrib(result_3, "group", axes=1) #barplot with contribution to first dimension (change values from "axes" to display other dimensions)

fviz_mfa_ind(result_3, col.ind = "cos2", 
             gradient.cols = c("red", "green", "blue"),
             repel = TRUE) #plot individuals

fviz_mfa_ind(result_3, 
             habillage = "Treatment", 
             palette = c("red", "green", "blue", "brown"),
             addEllipses = TRUE, ellipse.type = "confidence", 
             repel = TRUE) # plot individuals with colour by treatment

fviz_mfa_ind(result_3, 
             habillage = "Type", 
             palette = c("red", "green", "blue"),
             addEllipses = TRUE, ellipse.type = "confidence", 
             repel = TRUE) # plot individuals with colour by type

fviz_mfa_ind(result_3, 
             habillage = "Variety", 
             palette = c("red", "green", "blue", 
                         "brown", "pink", "yellow", 
                         "black", "orange", "lightblue", 
                         "darkblue", "darkgreen", "lightgreen",
                         "turquoise", "grey"),
             addEllipses = TRUE, ellipse.type = "confidence", 
             repel = TRUE) # plot individuals with colour by variety

res.var$cos2 #display quality of representation on the factor map

fviz_ellipses(result_3, c("Treatment", "Variety"), repel = TRUE) #plot Treatment and Variety together as Eliipses

fviz_mfa_axes(result_3, repel=T) #plot relationship between principal axes and results

