#Ryegrass Experiment - Biomass
#Dana Looschelders
#Changelog: 
  #03.02.2020: created
  #04.02.2020: tidy up script, add standard deviation to plots
#The script performs the following tasks: graph a singe measurement
  #aggregate data per variety/chamber
  #graphs for one measurement: barplot with linegraph for mean

setwd("C:/00 Dana/Uni/Internship/Work") #set working directory (location of files)
library(ggplot2) #load package 
#import data
Biomass=read.table("Dried Biomass Data Collection.csv", sep=";", dec=".", header=T)
str(Biomass) #check data
#delete empty cloumns off data.frame
Biomass=Filter(function(x)!all(is.na(x)), Biomass)
str(Biomass)#check again

#add factor for average values for chambers 1/2 and 7/8
Biomass$Chamber_avg=c(rep("control", 168), rep("elevated", 168))

#calculate mean Biomass for every variety in every chamber (4 groups) 
stats_Biomass_Cut_1=aggregate(Biomass$Biomass.Cut.1..23.01.2020., by=list(Varietys=Biomass$Variety, Chamber=Biomass$Chamber), mean)
names(stats_Biomass_Cut_1)[3]="avg_Biomass_cut_1"
#calculate standard deviation
sd=aggregate(Biomass$Biomass.Cut.1..23.01.2020., by=list(Varietys=Biomass$Variety, Chamber=Biomass$Chamber), sd)
stats_Biomass_Cut_1$sd=sd$x #add sd to dataframe
str(stats_Biomass_Cut_1) #check data

#calculate for 2 groups (1/2 and 7/8)
stats2_Biomass_Cut_1=aggregate(Biomass$Biomass.Cut.1..23.01.2020., by=list(Varietys=Biomass$Variety, Chamber=Biomass$Chamber_avg), mean)
names(stats2_Biomass_Cut_1)[3]="avg_Biomass_cut_1"
#calculate standard deviation 
sd2=aggregate(Biomass$Biomass.Cut.1..23.01.2020., by=list(Varietys=Biomass$Variety, Chamber=Biomass$Chamber_avg), sd)
stats2_Biomass_Cut_1$sd=sd2$x #add sd to dataframe
str(stats2_Biomass_Cut_1) #check data

#calculate mean per chamber (4 groups)
mean_Biomass_chamber=aggregate(stats_Biomass_Cut_1$avg_Biomass_cut_1, 
                               by = list(Chamber=stats_Biomass_Cut_1$Chamber), mean) #calculate mean per chamber
str(mean_Biomass_chamber) #check data

#add mean per chamber (2 groups)
mean_Biomass_chamber_avg2=aggregate(stats2_Biomass_Cut_1$avg_Biomass_cut_1, 
                                    by = list(Chamber=stats2_Biomass_Cut_1$Chamber), mean) #calculate mean per chamber
str(mean_Biomass_chamber_avg2) #check data


#plot four groups (every Chamber)
str(Biomass) #check data
Biomass$Chamber=as.factor(Biomass$Chamber) #Chambers as factors
H_bars=ggplot(Biomass, aes(Variety, Biomass.Cut.1..23.01.2020.)) #create bar plot layer
H_lines=geom_hline(data=mean_Biomass_chamber, aes(yintercept=x), color="red") #create line layer for mean
H_bars + stat_summary(fun.y=mean, geom="bar")+ #use stat_summaery to calculate mean
  theme_minimal()+ #select design for plot
  theme(axis.text.x = element_text(angle = 90))+ #rotate x axis labels (variety names)
  H_lines+ylab(label="Average Biomass [g]")+ #label for y axis
  xlab(label="Ryegrass Varieties")+ #label for x axis
  labs(title="Lolium perenne biomass 1 Cut Jan 2020", #choose main title
       subtitle = "Chambers: 1,2,7,8 \n red bar: mean per chamber") + #choose sub title
  stat_summary(fun.data=mean_se, geom="errorbar")+ #add standard error as errorbars
  facet_grid(.~Chamber) #split by defined groups

#plot for average for 1/2 and 7/8
str(Biomass) #check data
Biomass$Chamber=as.factor(Biomass$Chamber_avg) #groups as factor
colnames(mean_Biomass_chamber_avg2)=c("Chamber_avg", "x") #rename columns for mean (the column name it is split by has to be the same as specified in facet.grid!)
H_bars2=ggplot(Biomass, aes(Variety, Biomass.Cut.1..23.01.2020.)) #create bar plot layer
H_lines2=geom_hline(data=mean_Biomass_chamber_avg2, aes(yintercept=x), color="red") #create line layer for mean
H_bars2 + stat_summary(fun.y=mean, geom="bar")+ #use stat_summaery to calculate mean
  theme_minimal()+ #select design for plot
  theme(axis.text.x = element_text(angle = 90))+ #rotate x axis labels (variety names)
  H_lines2+ylab(label="Average Biomass [g]")+ #label for y axis
  xlab(label="Ryegrass Varieties")+ #label for x axis
  labs(title="Lolium perenne biomass 1 Cut Jan 2020", #main title
       subtitle = "Chambers: control: 1,2  elevated: 7,8 \n red bar: mean per chambers")+ #subtitle
  stat_summary(fun.data=mean_se, geom="errorbar")+ #add standard error as errorbars
  facet_grid(.~Chamber_avg) #split by defined groups


#Export stats
write.csv(stats_Biomass_Cut_1, "Biomass_statistics.csv")

