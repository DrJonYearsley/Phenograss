#Data analysis ryegrass
#Dana Looschelders
#Change log:
  #28.01.2020: Calculate mean days for leaf apperance after germination, percentage of plants without 1/2/3 leaf
  #30.01.2020 Calculate and graph difference in days for leaf apperance
  #04.02.2020 tidy up script

require(ggplot2) #load package
setwd("C:/00 Dana/Uni/Internship/Work") #set working directory (place of files)
Leaf=read.table("Leaf Phenology Data Collection.csv", header=TRUE, dec=".", sep=";")
str(Leaf) #check data
#delete empty cloumns off data.frame
Leaf=Filter(function(x)!all(is.na(x)), Leaf)
str(Leaf)#check again

#Seeds sown on 16.12.2019
#express dates of leaf appereance as days passed since germination
  #convert dates from factor to posixct
Leaf$X1st.Leaf=strptime(Leaf$X1st.Leaf, "%d.%m.%Y")
Leaf$X2nd.Leaf=strptime(Leaf$X2nd.Leaf, "%d.%m.%Y")
Leaf$X3rd.Leaf=strptime(Leaf$X3rd.Leaf, "%d.%m.%Y")

start_date=as.POSIXlt("2019-12-16") #create start date (date of sowing)
str(start_date) #check data
str(Leaf$X1st.Leaf) #check data

#calculate days until 1 leaf
Leaf$days_1st.Leaf=difftime(Leaf$X1st.Leaf, start_date, unit="days") 
Leaf$days_1st.Leaf=as.numeric(Leaf$days_1st.Leaf) #change class to numeric
#calculate days until 2 leaf
Leaf$days_2nd.Leaf=difftime(Leaf$X2nd.Leaf, start_date, unit="days") 
Leaf$days_2nd.Leaf=as.numeric(Leaf$days_2nd.Leaf) #change class to numeric
#calculate days until 3 leaf
Leaf$days_3rd.Leaf=difftime(Leaf$X3rd.Leaf, start_date, unit="days")
Leaf$days_3rd.Leaf=as.numeric(Leaf$days_3rd.Leaf) #change class to numeric

#aggregate mean leaf apperance by variety and chamber nr
#for first leaf
stats_1_leaf=aggregate(Leaf$days_1st.Leaf, by=list(Varieties=Leaf$Variety, Chamber=Leaf$Chamber), mean, na.rm=TRUE)
names(stats_1_leaf)[names(stats_1_leaf) == "x"]="mean_Leaf" #rename column
#calculate means for chamber
mean_1_Leaf_chamber=aggregate(stats_1_leaf$mean_Leaf, 
                              by = list(Chamber=stats_1_leaf$Chamber), mean) 
#for 2 leaf
stats_2_leaf=aggregate(Leaf$days_2nd.Leaf,by=list(Varieties=Leaf$Variety, Chamber=Leaf$Chamber), mean, na.rm=TRUE)
names(stats_2_leaf)[names(stats_2_leaf) == "x"]="mean_Leaf" #rename column
#calculate means for chamber
mean_2_Leaf_chamber=aggregate(stats_2_leaf$mean_Leaf, 
                              by = list(Chamber=stats_2_leaf$Chamber), mean) 

#for 3 leaf
stats_3_leaf=aggregate(Leaf$days_3rd.Leaf,by=list(Varieties=Leaf$Variety, Chamber=Leaf$Chamber), mean, na.rm=TRUE)
names(stats_3_leaf)[names(stats_3_leaf) == "x"]="mean_Leaf" #rename 
mean_3_Leaf_chamber=aggregate(stats_3_leaf$mean_Leaf, 
                              by = list(Chamber=stats_3_leaf$Chamber), mean) 


#plot
L_plot=ggplot(stats_1_leaf, aes(y=mean_Leaf, x=Varieties)) #create plot
L1_points=geom_point(aes(color="First"))
L1_lines=geom_hline(data=mean_1_Leaf_chamber, aes(yintercept=x), color="lightgreen") #create line layer for mean
L2_points=geom_point(data=stats_2_leaf, aes(color="Second"))
L2_lines=geom_hline(data=mean_2_Leaf_chamber, aes(yintercept=x), color="darkgreen") #create line layer for mean
L3_points=geom_point(data=stats_3_leaf, aes(color="Third"))
L3_lines=geom_hline(data=mean_3_Leaf_chamber, aes(yintercept=x), color="turquoise") #create line layer for mean
L_plot+L1_points+L1_lines+L2_points+L2_lines+L3_points+L3_lines+
  facet_grid(.~Chamber)+
  labs(title="Lolium perenne Phenology 30th Jan 2020", subtitle = " Chambers: 1,2,7,8 \n Seeds sown on 16th Dec 2019") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90))+
  ylab(label="Duration since seeds sown [days]")+
  xlab(label="Ryegrass Varieties")+
  scale_colour_manual(name=" Legend: \n Appearance of leaves", values=c(First="lightgreen", Second="darkgreen", Third="turquoise"))

 #calculate number of plants without second/third leaf 
varieties=unique(Leaf$Variety)
leafs=seq(1:3)
str(varieties)

#Function to count NAs per chamber per variety
count.nas <- function(Chamber=One, Chambernr=1){
Chamber=data.frame(na_values_leaf1=rep(0,14),na_values_leaf2=rep(0,14), 
                  na_values_leaf3=rep(0,14), total_plants=rep(0,14), 
                  na_percent_Leaf1=rep(0,14),na_percent_Leaf2=rep(0,14), 
                  na_percent_Leaf3=rep(0,14), "varieties"=varieties,row.names = varieties) #create data.frame to store result
for (i in varieties) {
  Chamber[i,1]=sum(is.na(Leaf$X1st.Leaf[Leaf$Variety==i&Leaf$Chamber==Chambernr]))
  Chamber[i,2]=sum(is.na(Leaf$X2nd.Leaf[Leaf$Variety==i&Leaf$Chamber==Chambernr]))
  Chamber[i,3]=sum(is.na(Leaf$X3rd.Leaf[Leaf$Variety==i&Leaf$Chamber==Chambernr]))
  Chamber[i,4]=sum(Leaf$Variety==i&Leaf$Chamber==Chambernr)
}
for(i in leafs) {
  Chamber[,4+i]=Chamber[,i]/Chamber[,4]
}
return(Chamber)
}

#call function to calculate
Chamber_1=count.nas()
Chamber_2=count.nas(Chamber=Two, Chambernr=2)
Chamber_7=count.nas(Chamber=Seven, Chambernr=7)
Chamber_8=count.nas(Chamber=Eight, Chambernr=8)

#plot percentage of plants without leaf
old.par <- par(mfrow=c(2, 1))
plot(Chamber_1$varieties, Chamber_1$na_values_leaf3, ylab="Plants without leaf [%]", xlab= "Variety", main="Chamber 1")
plot(Chamber_2$varieties, Chamber_2$na_values_leaf3, ylab="Plants without leaf [%]", xlab= "Variety", main="Chamber 2")
par(old.par)

old.par <- par(mfrow=c(2, 1))
plot(Chamber_7$varieties, Chamber_7$na_values_leaf3, ylab="Plants without leaf [%]", xlab= "Variety", main="Chamber 7")
plot(Chamber_8$varieties, Chamber_8$na_values_leaf3, ylab="Plants without leaf [%]", xlab= "Variety", main="Chamber 8")
par(old.par)

#Calculate difference in days for leaf one
diff_days=data.frame(Leaf, "Dif_1_2"=rep(0, length(Leaf[,1])), "Dif_2_3"=rep(0, length(Leaf[,1])))
diff_days$Dif_1_2=as.numeric(difftime(Leaf$X2nd.Leaf, Leaf$X1st.Leaf, units="days"))
diff_days$Dif_2_3=as.numeric(difftime(Leaf$X3rd.Leaf, Leaf$X2nd.Leaf, units="days"))

stats_0_to_1=aggregate(diff_days$days_1st.Leaf, by=list(varieties=diff_days$Variety, Chamber=diff_days$Chamber), mean, na.rm=TRUE)
names(stats_0_to_1)[names(stats_0_to_1) == "x"]="mean_days_to_1_leaf"
stats_1_to_2=aggregate(diff_days$Dif_1_2, by=list(Varieties=diff_days$Variety, Chamber=diff_days$Chamber), mean, na.rm=TRUE)
names(stats_1_to_2)[names(stats_1_to_2) == "x"]="mean_days_to_2_leaf" #rename column
stats_2_to_3=aggregate(diff_days$Dif_2_3, by=list(Varieties=diff_days$Variety, Chamber=diff_days$Chamber), mean, na.rm=TRUE)
names(stats_2_to_3)[names(stats_2_to_3) == "x"]="mean_days_to_3_leaf" #rename column

#regroup the percentage data to fit in plot
percent_leaves=data.frame(Leaf, "percent_leaves"=rep(0, length(Leaf[,1])))
percent_leaves$percent_leaves[1:length(Chamber_1[,1])]=Chamber_1$na_percent_Leaf3


#reorganise data.frame to plot philochrom and percent of plant with leaves
percent_leaves=data.frame(stats_1_leaf)
str(percent_leaves)
percent_leaves$na_percent_third_leaf[1:length(Chamber_1[,1])]=Chamber_1$na_percent_Leaf3 
percent_leaves$na_percent_third_leaf[15:28]=Chamber_2$na_percent_Leaf3
percent_leaves$na_percent_third_leaf[29:42]=Chamber_7$na_percent_Leaf3
percent_leaves$na_percent_third_leaf[43:56]=Chamber_8$na_values_leaf3

percent_leaves$na_percent_second_leaf[1:length(Chamber_1[,1])]=Chamber_1$na_percent_Leaf2 
percent_leaves$na_percent_second_leaf[15:28]=Chamber_2$na_percent_Leaf2
percent_leaves$na_percent_second_leaf[29:42]=Chamber_7$na_percent_Leaf2
percent_leaves$na_percent_second_leaf[43:56]=Chamber_8$na_values_leaf2

percent_leaves$na_percent_first_leaf[1:length(Chamber_1[,1])]=Chamber_1$na_percent_Leaf1 
percent_leaves$na_percent_first_leaf[15:28]=Chamber_2$na_percent_Leaf1
percent_leaves$na_percent_first_leaf[29:42]=Chamber_7$na_percent_Leaf1
percent_leaves$na_percent_first_leaf[43:56]=Chamber_8$na_values_leaf1

#add values to leaf specific data.frames
stats_0_to_1$percent_na=percent_leaves$na_percent_first_leaf
stats_1_to_2$percent_na=percent_leaves$na_percent_second_leaf
stats_2_to_3$percent_na=percent_leaves$na_percent_third_leaf
#plot phyllocron and percent of plants with no leaves together

#first leaf
Dif_plot_3=ggplot(stats_0_to_1, aes(y=percent_na, x=varieties))
Dif_plot_3+geom_bar(stat="identity", fill="blue")+
  geom_point(data=stats_0_to_1, aes(x=varieties, y=mean_days_to_1_leaf))+
  facet_grid(.~Chamber)+
  labs(title="Lolium perenne Phyllochron 30th Jan 2020", subtitle = " Chambers: 1,2,7,8 \n Seeds sown on 16th Dec 2019") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90), axis.text.y.right=element_text(color = "blue"), axis.title.y.right = element_text(color="blue"))+
  scale_y_continuous(name="Duration between germination and 1st leaf [days]", sec.axis=sec_axis(~.*1, name="Plants without 2nd leaf [%]"))+
  xlab(label="Ryegrass Varieties")

#second leaf
Dif_plot_3=ggplot(stats_1_to_2, aes(y=percent_na, x=Varieties))
Dif_plot_3+geom_bar(stat="identity", fill="blue")+
  geom_point(data=stats_1_to_2, aes(x=Varieties, y=mean_days_to_2_leaf))+
  facet_grid(.~Chamber)+
  labs(title="Lolium perenne Phyllochron 30th Jan 2020", subtitle = " Chambers: 1,2,7,8 \n Seeds sown on 16th Dec 2019") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90), axis.text.y.right=element_text(color = "blue"), axis.title.y.right = element_text(color="blue"))+
  scale_y_continuous(name="Duration between 1st and 2nd leaf [days]", sec.axis=sec_axis(~.*1, name="Plants without 2nd leaf [%]"))+
  xlab(label="Ryegrass Varieties")

#third leaf
Dif_plot_3=ggplot(stats_2_to_3, aes(y=percent_na, x=Varieties))
Dif_plot_3+geom_bar(stat="identity", fill="blue")+
  geom_point(data=stats_2_to_3, aes(x=Varieties, y=mean_days_to_3_leaf))+
  facet_grid(.~Chamber)+
  labs(title="Lolium perenne Phyllochron 30th Jan 2020", subtitle = " Chambers: 1,2,7,8 \n Seeds sown on 16th Dec 2019") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90), axis.text.y.right=element_text(color = "blue"), axis.title.y.right = element_text(color="blue"))+
  scale_y_continuous(name="Duration between 2nd and 3rd leaf [days]", sec.axis=sec_axis(~.*1, name="Plants without 3rd leaf [%]"))+
  xlab(label="Ryegrass Varieties")

#in percent plants with leaves
stats_0_to_1$reverse_percent_na=(100-stats_0_to_1$percent_na)*0.3 #add 0.3 as scaling factor just for test display

Dif_plot_3=ggplot(stats_0_to_1, aes(y=reverse_percent_na, x=varieties))
Dif_plot_3+geom_bar(stat="identity", fill="lightblue")+
  geom_point(data=stats_0_to_1, aes(x=varieties, y=mean_days_to_1_leaf))+
  facet_grid(.~Chamber)+
  labs(title="Lolium perenne Phyllochron 30th Jan 2020", subtitle = " Chambers: 1,2,7,8 \n Seeds sown on 16th Dec 2019") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90), axis.text.y.right=element_text(color = "blue"), axis.title.y.right = element_text(color="blue"))+
  scale_y_continuous(name="Duration between germination and 1st leaf [days]", sec.axis=sec_axis(~.*1, name="Plants without 2nd leaf [%]"))+
  xlab(label="Ryegrass Varieties")

stats_2_to_3$reverse_percent_na=(100-stats_2_to_3$percent_na)*0.3 #add 0.3 as scaling factor just for test display

Dif_plot_3=ggplot(stats_2_to_3, aes(y=reverse_percent_na, x=Varieties))
Dif_plot_3+geom_bar(stat="identity", fill="lightblue")+
  geom_point(data=stats_2_to_3, aes(x=Varieties, y=mean_days_to_3_leaf))+
  facet_grid(.~Chamber)+
  labs(title="Lolium perenne Phyllochron 30th Jan 2020", subtitle = " Chambers: 1,2,7,8 \n Seeds sown on 16th Dec 2019") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90), axis.text.y.right=element_text(color = "blue"), axis.title.y.right = element_text(color="blue"))+
  scale_y_continuous(name="Duration between 2nd and 3rd leaf [days]", sec.axis=sec_axis(~.*1, name="Plants without 3rd leaf [%]"))+
  xlab(label="Ryegrass Varieties")

