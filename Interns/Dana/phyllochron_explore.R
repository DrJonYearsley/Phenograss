# Trial analysis/visualisation of phyllochron data
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# 16 May 2020
# ******************************************

rm(list=ls())

library(readxl)
library(tidyr)
library(ggplot2)

# download data and replace spaces in names with '.'
phyllo12 = read_excel('./data_rosemount/Leaf Phyllochron 1-2.xlsx', .name_repair = 'universal')
phyllo23 = read_excel('./data_rosemount/Leaf Phyllochron 2-3.xlsx', .name_repair = 'universal')
phyllo34 = read_excel('./data_rosemount/Leaf Phyllochron 3-4.xlsx', .name_repair = 'universal')


# Put data in long format
phyllo12_long = pivot_longer(phyllo12,
                             cols = -c(1:4),
                             names_to = 'month',
                             values_to = 'days',
                             values_drop_na = TRUE)
phyllo23_long = pivot_longer(phyllo23,
                             cols = -c(1:4),
                             names_to = 'month',
                             values_to = 'days',
                             values_drop_na = TRUE)
phyllo34_long = pivot_longer(phyllo34,
                             cols = -c(1:4),
                             names_to = 'month',
                             values_to = 'days',
                             values_drop_na = TRUE)

# Combine all the data into one data frame
phyllo12_long$phylloID = '1-2'
phyllo23_long$phylloID = '2-3'
phyllo34_long$phylloID = '3-4'

phyllo = rbind(phyllo12_long, phyllo23_long, phyllo34_long)

# Clean up unused data frames
rm(list=c('phyllo12', 'phyllo23','phyllo34'))

phyllo$Variety = as.factor(phyllo$Variety)
phyllo$Treatment = as.factor(phyllo$Treatment)
phyllo$Chamber = as.factor(phyllo$Chamber)
phyllo$month = factor(phyllo$month, levels=c('May','June','July','August','September'))
phyllo$phylloID = as.factor(phyllo$phylloID)
phyllo$Plant.ID = as.factor(phyllo$Plant.ID)


# Visualise the data
format = theme(axis.title = element_text(size=14),
               axis.text = element_text(size=12))



ggplot(data=phyllo, 
       aes(x=month, y=days, fill=phylloID)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75)) + 
  scale_fill_brewer(name='Phyllochron', palette='Set1') +
  labs(title='Both Treatments') +
  theme_bw() + 
  format

ggplot(data=subset(phyllo,phylloID='3-4'), 
       aes(x=month, y=days, fill=Treatment)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75)) + 
  scale_fill_brewer(name='Treatment', palette='Set1') +
  labs(title='Phyllochron 3-4') +
  theme_bw() + 
  format

ggplot(data=subset(phyllo,phylloID='3-4'), 
       aes(x=month, y=days, fill=Treatment)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), adjust=2) + 
  scale_fill_brewer(name='Treatment', palette='Set1') +
  labs(title='Phyllochron 3-4: Broader bandwidth') +
  theme_bw() + 
  format

ggplot(data=subset(phyllo,phylloID='3-4'), 
       aes(x=month, y=days, fill=Chamber)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75)) + 
  scale_fill_brewer(name='Chamber', 
                    palette='Set1', 
                    breaks=c(1,2,7,8),
                    labels=c('1: Ambient', 
                             '2: Treatment', 
                             '7:Ambient',
                             '8: Treatment')) +
  labs(title='Phyllochron 3-4') +
  theme_bw() + 
  format


ggplot(data=subset(phyllo,phylloID='1-2'), 
       aes(x=month, y=days, fill=Chamber)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75)) + 
  scale_fill_brewer(name='Chamber', 
                    palette='Set1', 
                    breaks=c(1,2,7,8),
                    labels=c('1: Ambient', 
                             '2: Treatment', 
                             '7:Ambient',
                             '8: Treatment')) +
  labs(title='Phyllochron 1-2') +
  theme_bw() + 
  format

