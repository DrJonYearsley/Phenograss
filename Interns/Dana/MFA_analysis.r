#Principal component analysis 
#rather use multiple factor analysis as there are both factorial and continious variables
#ressources used:
#https://www.youtube.com/watch?v=g5_hM93e8HM (as suggested by morane)
#http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/116-mfa-multiple-factor-analysis-in-r-essentials/#r-code
#by Dana Looschelders
#created: 28.03.2020
#changelog: 

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

#*************************************************************************
#try as PCA with Treatment and Variety as supplementary variables
#ressources used: 
  #http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/
dat=read.table("data_for_pca.csv", sep=";", dec=",", header=T)
str(dat)

#prepare the data for multiple factor analysis (mfa)
#set plant ID as rownames
rownames(dat)=dat$ï..Plant_ID
dat=dat[,-1] #drop column with ID as it is now the rowname
dat=dat[,-5] #drop column with first biomass values
dat$Chamber=rep(c("control", "elevated"), each= 168) #add treatment without waterlogging
colnames(dat)[4]="Treatment2"
dat[,4]=as.factor(dat[,4]) #coerce to factor
  #drop column with Chamber values? would it make sense to test for Chamber effects this way?
  #does it make sense to include both treatment categories? or drop the one with waterlogging?
  #then use chambers too?
dat=dat[,-2]
str(dat)

?PCA
res=PCA(dat, quali.sup=1:3, graph=T, scale.unit = T)
#remember: eigenvalue -> amount of variance retained by each principal component

get_eigenvalue(res) #extract eigenvalues
fviz_eig(res, addLabels=T) #visualize eigenvalues with scree plot

get_pca_ind(res) #extract results for indivduals
fviz_pca_ind(res) #visualize results for individuals

var=get_pca_var(res) #extract results for variables
fviz_pca_var(res, col.var="black", repel=T) #visualze results for variables
var$coord #coordinates of variables
var$cor #correlation variables and dimensions
var$cos2 #quality of representation for variables
  #high cos2 -> good representation of variable on PC
corrplot(var$cos2, is.corr=F) #correlation plot
fviz_cos2(res, choice="var") #barplot
fviz_pca_var(res, col.var="cos2", repel=T) #plot with color according to co2 value

var$contrib #contributions of the variables in %
fviz_pca_var(res, col.var="contrib", repel=T) #plot with color according to contribution to PC
corrplot(var$contrib, is.corr=F)
fviz_contrib(res, choice="var")

fviz_pca_biplot(res) #visualize both (add axes= e.g. 1 to see contribution to single dimension or 1:2 for multiple)

ind=get_pca_ind(res)
fviz_pca_ind(res, col.ind="cos2", repel=T)

res.dim=dimdesc(res, axes=c(1,2,3), proba = 0.05) #identify most significantly associated variables with given PC
res.dim$Dim.1
res.dim$Dim.2
res.dim$Dim.3

fviz_pca_ind(res, col.ind="cos2", repel=T) #plot individuals colored in relation to cos2 value
fviz_pca_ind(res, pointsize="cos2", repel=T) #plot individuals colored in relation to cos2 value

fviz_cos2(res, choice="ind", axes=1:2) #bar plot for dim 1 and 2 with contribution of individuals

#color by groups
#for Treatment
fviz_pca_ind(res, geom.ind="point", 
             col.ind=dat$Treatment2, addEllipses = T, 
             legend.title="Treatment")
#for Variety
fviz_pca_ind(res, geom.ind="point", 
             col.ind=dat$Variety, addEllipses = T,
             legend.title="Variety")
#for type
fviz_pca_ind(res, geom.ind="point", 
             col.ind=dat$Type, addEllipses = T,
             legend.title="Type")

#biplot with dimensions (grouped)
#by treatment
fviz_pca_biplot(res, col.ind=dat$Treatment, addEllipses = T,
                label="var", col.var="black", repel=T, legend.title="Treatment")
#by variety
fviz_pca_biplot(res, col.ind=dat$Variety, addEllipses = T,
                label="var", col.var="black", repel=T, legend.title="Variety")
#by treatment
fviz_pca_biplot(res, col.ind=dat$Type, addEllipses = T,
                label="var", col.var="black", repel=T, legend.title="Type")

res$quali.sup #predicited results for the supplementary qualitative variables
fviz_pca_ind(res, habillage = 3, #graph according to supplementary variables
             addEllipses =TRUE, ellipse.type = "confidence",
             repel = TRUE) 
#habillage= index of supplementary variable to show
    #1 = Variety
    #2 = Type
    #3 = Treatment

#plot selected variables and individuals
fviz_pca_var(res, select.var=list(cos2=0.6)) #visualize only variables with cos2>0.6

fviz_pca_biplot(res, select.ind=list(contrib=20)) #display the 20 most contributing individuals

#export the results to csv
write.infile(res, "pca.csv", sep=";")
save.image()
