# Latent variable model for growth chamber data
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# Feb 2021
# ******************************************

rm(list=ls())

library(readxl)
library(tidyr)
library(ggplot2)
library(lavaan)

# download data and replace spaces in names with '.'
biomass = read_excel('./data_rosemount/Mixed Growth Dry Weight.xlsx', 
                     .name_repair = 'universal')

height = read_excel('./data_rosemount/Mixed Maximum Plant Height.xlsx', 
                    .name_repair = 'universal')



d = merge(biomass, height)

dNames = names(d)
dNames[12] = 'Biomass_GR'
dNames[13] = 'Height_GR'
names(d) = dNames

d$Harvest = as.ordered(d$Harvest)
d$Chamber = as.factor(d$Chamber)
d$Treatment = as.factor(d$Treatment)
d$Variety = as.factor(d$Variety)


# Some data visualisation
ggplot(d, aes(x=Height_GR,
              y=Biomass_GR,
              colour=Harvest)) +
  geom_point() + 
  scale_colour_brewer(palette='Set1') +
  theme_bw()


ggplot(d, aes(x=Harvest,
              y=Biomass_GR,
              fill=Treatment)) + 
  geom_boxplot() +
  scale_fill_brewer(palette='Dark2') + 
  theme_bw()


ggplot(d, aes(x=Harvest,
              y=Height_GR,
              fill=Treatment)) + 
  geom_boxplot() +
  scale_fill_brewer(palette='Dark2') + 
  theme_bw()

ggplot(subset(d, Harvest!=1 & Variety%in%levels(Variety)[4]), 
       aes(y=Biomass_GR,
              x=Height_GR,
              fill=Harvest,
           colour=Harvest)) + 
  geom_point(shape=21, size=2) +
  geom_smooth(method=lm,formula=y~x) +
  scale_fill_brewer(palette='Dark2') + 
  scale_colour_brewer(palette='Dark2') + 
  theme_bw()


ggplot(subset(d, Harvest!=1 & Variety%in%levels(Variety)[4]), 
       aes(y=Biomass_GR,
           x=Temperature)) + 
  geom_point(shape=21, size=2) +
  geom_smooth(method=lm,formula=y~x) +
  scale_fill_brewer(palette='Dark2') + 
  scale_colour_brewer(palette='Dark2') + 
  theme_bw()


ggplot(subset(d, Harvest!=1 & Variety%in%levels(Variety)[4]), 
       aes(x=CO2,
           y=Temperature,
         z=Height_GR)) + 
  geom_contour() +
  scale_fill_brewer(palette='Dark2') + 
  scale_colour_brewer(palette='Dark2') + 
  theme_bw()



aggregate(cbind(Height_GR,Biomass_GR)~Temperature+CO2, data=d, FUN=mean)

# Try developing a latent variable model

dsub = subset(d, Variety%in%levels(Variety)[4] & Harvest!=1)


TreatmentID = model.matrix(~0+Treatment, data=dsub)
HarvestID = model.matrix(~0+Harvest, data=dsub)

dsub = cbind(dsub, HarvestID, TreatmentID)

summary(lm(Height_GR~1 + Harvest2 + Harvest3 + Harvest4, data=dsub))

mod = '# Regressions
       Biomass_GR ~ 1 + v1 
       Height_GR ~ 1 + v1
       v1 ~ Harvest2 + Harvest3 + Harvest4 

       v1 =~ Biomass_GR + Height_GR
'

m_sem1 = sem(model=mod, 
                data=dsub)

summary(m_sem1)






# Specify with lava package

library(lava)


m1 = lvm()
regression(m1) = Biomass_GR~Temperature + CO2 
regression(m1) = Height_GR~Temperature + CO2 
#covariance(m1) = Biomass_GR~Height_GR
#regression(m1) = v1~Harvest2 + Harvest3 + Harvest4
#latent(m1) = ~v1

e = estimate(m1, dsub)

summary(e)
