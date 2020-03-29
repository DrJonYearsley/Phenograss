#Ryegrass exeriment - Height 
#Dana Looschelders
#Changelog:
#30.01.2020: update data
#install.packages("ggplot2")
#install.packages("tidyr")
#install.packages("tidyverse")
#install.packages("dplyr")
#install.packages("installr")
library(ggplot2) #load packages
library(tidyr)
library(tidyverse)

#*********************************************************************
#choose your working directory (should contain a csv.file with data)
#**********************************************************************
setwd("C:/00 Dana/Uni/Internship/Work") #set working directory (place of files)
#import data
Height=read.table("Plant Height Data Collection.csv", sep=";", dec=".", header=T)
str(Height) #check data
#delete empty cloumns off data.frame
Height=Filter(function(x)!all(is.na(x)), Height)
str(Height)#check again

#jitter plots
Height_jitter=cbind(Height$Variety, Height$Chamber, Height$Height.Week.6..20.02.2020.)

plot(Height_jitter)
plot(Height)

#calculate mean height for every variety in every Treatment
#function for statistics for every week (just add data and run function)
#names: names for columns in output dataframe (calles stats_Height)
#no_col: number of columns to process (DON'T count first column as dataframe is created with it)

#*************************************************************************************************
#names: names for columns in output dataframe (calles stats_Tiller)
#no_col: number of columns to process (DON'T count first column as dataframe is created with it)
#***************************************************************************************************

stats.height=function(no_col=5, names=c("week 2", "week 3", "week 4", "week 5", "week 6")){
  cols=seq(from=6, by=1, to=6+no_col-1)
  stats_Height=aggregate(Height[,5], by=list(Varietys=Height$Variety, Treatment=Height$Treatment), mean)
  names(stats_Height)[3]="week 1"
  for (i in cols) {
    dummy_stats=aggregate(Height[,i], by=list(Varietys=Height$Variety, Treatment=Height$Treatment), mean)
    col=i-2
    stats_Height[,col]=dummy_stats$x
    names(stats_Height)[col]=names[i-5]
  }
  print(stats_Height)
}

stats_Height=stats.height()
write.csv(x = stats_Height, file="stats_height.csv")
#create time series
#***********************************************************************************************************
#to modify for new data: add Treatment_x$'weex x' to function rbind and and the number for the week to cbind
#************************************************************************************************************
unique(stats_Height$Treatment)
get.height.time.series=function(Treatment="CON"){
  Treatment=stats_Height[stats_Height$Treatment==Treatment,]
  Treatment_dummy=rbind(Treatment$`week 1`, Treatment$`week 2`, Treatment$`week 3`, Treatment$`week 4`, Treatment$`week 5`, Treatment$`week 6`)
  colnames(Treatment_dummy)=Treatment$Varietys
  Treatment_dummy=cbind(Treatment_dummy,week=c(1,2,3,4,5,6))
  Treatment_dummy=as.data.frame(Treatment_dummy)
  #the semi-natural varieties need to be renamed because the - messes with R (an underscore is fine)
  names(Treatment_dummy)[8]="semi_natural11" 
  names(Treatment_dummy)[9]="semi_natural16"
  names(Treatment_dummy)[10]="semi_natural17"
  #collapse the measures for the different varieties into one column in order to plot the time series
  df <- Treatment_dummy %>% #creates new dataframe
    select(Aberchoice, Abergain, Aspect, Carraig, Dunluce, Lilora, Moy, semi_natural11, semi_natural16, semi_natural17,Solomon, Wild4, Wild6, Wild7, week) %>%
    gather(key = "variable", value = "value", -week) 
  #plot as time series
  ggplot(df, aes(x = week, y = value)) + #basic plot 
    geom_line(aes(color = variable), size = 1) + #adds lines
    geom_point(aes(shape=variable), size=3)+ #adds points
    theme_minimal()+ #defines design of plot
    scale_shape_manual(values=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14))+ #changes shape of points
    ylab(label="Plant Height [cm]")+ #labels y axes
    xlab(label="Week")+ #labels x axes
    labs(title=paste("Treatment",Treatment,"Lolium perenne height 07. Feb 2020"), #title of plot
         subtitle = "bar: biomass cut")+ #subtitle of plot
    geom_vline(xintercept=2)+ #adds bar to indicate biomass cut (intercept determines where bar is set)
    geom_vline(xintercept=5)
}
get.height.time.series(Treatment="CON")

Treatments=c("CON", "WAT", "eCO2", "eCO2W")




#function to save plots for all four Treatments
save.height.plot=function(){
  for (i in Treatments){
    title=paste("Height",i,".pdf") 
    get.height.time.series(Treatment=i)
    ggsave(title, width=12, height=8)
  }
}
save.height.plot() #call function to execute

#for a single date plot with ggplot2
#example with 30th Jan
H_plot3=ggplot(stats_Height_30, aes(y=avg_Height_30 ,x=Varietys)) #create plot
H_bars3=geom_bar(stat="identity") #create bar plot layer
#add mean per Treatment
mean_height_Treatment3=aggregate(stats_Height_30$avg_Height_30, 
                               by = list(Treatment=stats_Height_30$Treatment), mean) #calculate mean per Treatment
str(mean_height_Treatment) #check data
H_lines3=geom_hline(data=mean_height_Treatment3, aes(yintercept=x), color="red") #create line layer for mean
H_plot3+H_bars3+facet_grid(.~Treatment)+
  labs(title="Lolium perenne Height 30. Jan 2020", 
       subtitle = "Treatments: 1,2,7,8 \n red bar: mean per Treatment") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90))+
  H_lines3+ylab(label="Average Height [cm]")+
  xlab(label="Ryegrass Varieties")


#plot with growth
#re-arrange dataframe in order to display growth as stacked bars
stats_Height_dummy=data.frame(Dummy=rep(1,112), "Height"=stats_Height$avg_Height_17, "Varieties"=stats_Height$Varietys, "Treatment"=stats_Height$Treatment)
stats_Height_dummy[57:112,2]=stats_Height$Growth
stats_Height_dummy[57:112,3]=stats_Height$Varietys
stats_Height_dummy[57:112,4]=stats_Height$Treatment
stats_Height_dummy[57:112,1]=rep(2,56)
str(stats_Height_dummy)
#plot with stacked bars
H_plot3=ggplot(stats_Height_dummy, aes(fill=Dummy,y=Height ,x=Varieties)) #create plot
H_bars3=geom_bar(stat="identity", position="stack") #create bar plot layer
H_lines=geom_hline(data=mean_height_Treatment, aes(yintercept=x), color="grey") #create line layer for mean
H_lines2=geom_hline(data=mean_height_Treatment2, aes(yintercept=x), color="black") #create line layer for mean
H_plot3+H_bars3+facet_grid(.~Treatment)+
  labs(title="Lolium perenne growth 23rd Jan 2020",
       subtitle = "Treatments: 1,2,7,8\n darkblue: growth until 17th \n lightblue: growth until 23d \n grey bar: mean per Treatment on 17th \n black bar: mean per Treatment on 23rd") +
  theme_minimal()+theme(axis.text.x = element_text(angle = 90))+
  H_lines+H_lines2+ylab(label="Average Height [cm]")+
  xlab(label="Ryegrass Varieties")+guides(fill=FALSE)

#Export stats
write.csv(stats_Height, "Height_statistics.csv")
