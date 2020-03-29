#Ryegrass exeriment - Tiller number
#Dana Looschelders
#Changelog:
#30.01.2020: update data
Tiller=read.table("Tiller_number.csv", sep=";", dec=".", header=T)
str(Tiller)
library(tidyverse) #load package 

###Average Tiller number
#for the 17th
stats_Tiller_17=aggregate(Tiller$Tiller..N..17.01.2020, by=list(Varietys=Height$Variety, Chamber=Height$Chamber), mean)
names(stats_Tiller_17)[3]="avg_Tiller_17"
#plot with ggplot2
T_plot=ggplot(stats_Tiller_17, aes(y=avg_Tiller_17 ,x=Varietys)) #create plot
T_bars=geom_bar(stat="identity") #create bar plot layer
#add mean per chamber
mean_tiller_chamber=aggregate(stats_Tiller_17$avg_Tiller_17, 
                              by = list(Chamber=stats_Tiller_17$Chamber), mean) #calculate mean per chamber
T_lines=geom_hline(data=mean_tiller_chamber, aes(yintercept=x), color="red") #create line layer for mean
T_plot+T_bars+facet_grid(.~Chamber)+labs(title="Lolium perenne tiller number 17th Jan 2020", 
                                         subtitle = "Chambers: 1,2,7,8 \n red bar: mean per chamber") +theme_minimal()+theme(axis.text.x = element_text(angle = 90))+T_lines+ylab(label="Average Tiller number")+xlab(label="Ryegrass Varieties")
#for the 23rd
stats_Tiller_23=aggregate(Tiller$Tiller..N..23.01.2020, by=list(Varietys=Height$Variety, Chamber=Height$Chamber), mean)
names(stats_Tiller_23)[3]="avg_Tiller_23"

#for the 30th
stats_Tiller_30=aggregate(Tiller$Tiller.N.Cut.3.30.01.2020, by=list(Varietys=Height$Variety, Chamber=Height$Chamber), mean)
names(stats_Tiller_30)[3]="avg_Tiller_30"


#plot with ggplot2
T_plot2=ggplot(stats_Tiller_30, aes(y=avg_Tiller_30 ,x=Varietys)) #create plot
T_bars2=geom_bar(stat="identity") #create bar plot layer
#add mean per chamber
mean_tiller_chamber=aggregate(stats_Tiller_30$avg_Tiller_30, 
                              by = list(Chamber=stats_Tiller_30$Chamber), mean) #calculate mean per chamber
T_lines2=geom_hline(data=mean_tiller_chamber, aes(yintercept=x), color="red") #create line layer for mean
T_plot2+T_bars2+facet_grid(.~Chamber)+labs(title="Lolium perenne tiller number 30rd Jan 2020", 
                                           subtitle = "Chambers: 1,2,7,8 \n red bar: mean per chamber") +theme_minimal()+theme(axis.text.x = element_text(angle = 90))+T_lines2+ylab(label="Average Tiller number")+xlab(label="Ryegrass Varieties")

