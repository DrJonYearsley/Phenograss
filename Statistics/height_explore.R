# Trial analysis/visualisation of height data
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# 16 May 2020
# ******************************************

rm(list=ls())

library(readxl)
library(tidyr)
library(ggplot2)

# download data and replace spaces in names with '.'
height = read_excel('./data_rosemount/Plant Height Data.xlsx', .name_repair = 'universal')


# Put data in long format
height_long = pivot_longer(height,
                           cols = -c(1:4),
                           names_to = 'date',
                           values_to = 'height',
                           values_drop_na = TRUE,
                           names_pattern = "([0-9]{2}[[:punct:]]{1}[0-9]{2}[[:punct:]][0-9]{4})")

# Calculate days from 16th Dec 2019
height_long$date = as.Date(height_long$date, format="%d.%m.%Y")
height_long$days = julian(height_long$date, origin=as.Date('2019-12-16'))

# Biomass was cut week 5, 9 and 13
height_long$cut = 'Cut_1'

ind = difftime(height_long$date, as.Date('2020-02-13'), units='days')<=0
height_long$cut[ind] = 'Cut_0'

# ind = difftime(height_long$date, as.Date('2020-02-13'), units='days')>0 &
#   difftime(height_long$date, as.Date('2020-03-12'), units='days')<=0
# height_long$cut[ind] = 'Cut_1'

ind = difftime(height_long$date, as.Date('2020-03-12'), units='days')>0 &
  difftime(height_long$date, as.Date('2020-04-09'), units='days')<=0
height_long$cut[ind] = 'Cut_2'

ind = difftime(height_long$date, as.Date('2020-04-09'), units='days')>0
height_long$cut[ind] = 'Cut_3'

# Recode variables
height_long$cut = as.factor(height_long$cut)
height_long$Variety = as.factor(height_long$Variety)
height_long$Treatment = as.factor(height_long$Treatment)
height_long$Chamber = as.factor(height_long$Chamber)

# Calculate growth rates
tmp = aggregate(cbind(height,days)~Plant.ID+Variety+Treatment+cut, data=height_long, FUN=diff)
tmp2 = aggregate(days~Plant.ID+Variety+Treatment+cut, data=height_long, FUN=min)


# Construct final dataframe
for (r in 1:nrow(tmp)) {
  if (r==1) {
    gr = data.frame(Plant.ID=tmp$Plant.ID[r],
                    Variety = tmp$Variety[r],
                    Treatment = tmp$Treatment[r],
                    cut = tmp$cut[r],
                    time = tmp2$days[r]+cumsum(tmp$days[[r]]),
                    growth.rate =  tmp$height[[r]]/tmp$days[[r]])
  } else {
    gr = rbind(gr,
               data.frame(Plant.ID=tmp$Plant.ID[r],
                    Variety = tmp$Variety[r],
                    Treatment = tmp$Treatment[r],
                    cut = tmp$cut[r],
                    time = tmp2$days[r]+cumsum(tmp$days[[r]]),
                    growth.rate =  tmp$height[[r]]/tmp$days[[r]]))
  }
}


head(gr)
tail(gr)

# ******************************************************
# Visualise the data
format = theme(axis.title = element_text(size=14),
               axis.text = element_text(size=12))



table(height_long$days)


# Average of varieties and plants
d = aggregate(height~days+Treatment+cut, data=height_long, FUN=median)

ggplot(data=d, 
       aes(x=days, y=height, colour=Treatment, fill=cut, group=Treatment)) +
  geom_path(size=1) +
  geom_point(shape=21, size=3) +
  scale_fill_brewer(name='Cuts', palette='Dark2') +
  scale_colour_brewer(name='Treatment', palette='Set1') +
  labs(title='Average over plant ID') +
  theme_bw() + 
  format


ggplot(d=gr, 
       aes(y=growth.rate, x=cut, fill=Treatment)) +
  geom_boxplot() +
  scale_fill_brewer(name='Treatment', palette='Set1') +
  labs(title='Growth rates (height gain per day)') +
  theme_bw() + 
  format

ggplot(d=subset(gr, Treatment=='Ambient'), 
       aes(y=growth.rate, x=factor(time), fill=cut)) +
  geom_boxplot() +
  scale_fill_brewer(name='Treatment', palette='Set1') +
  labs(title='Growth rates (height gain per day)') +
  theme_bw() + 
  format


# *****************************
# Look at some initial analyses

library(DHARMa)
library(tweedie)
library(lme4)
library(performance)
library(broom)
library(merTools)

m = glm(growth.rate~cut*time*Treatment, 
        data=subset(gr,cut!='Cut 3'), 
        family=gaussian)

sim = simulateResiduals(m, n=200)
plot(sim)

plot(m)  # Some over-dispersion, otherwise pretty good
summary(m)

# Remove three way interaction
m0 = update(m, .~.-cut:time:Treatment)

anova(m0,m, test='F')  # Slight 3 way interaction

# *******************************
# Include Variety as a random effect on intercept
m2 = lmer(growth.rate~time*Treatment +(1|Variety), 
          data=subset(gr,cut!='Cut 3'),
          REML=FALSE)

# Validate model (DHARMa)
sim_lme = simulateResiduals(m2, n=200)
plot(sim_lme)

#performance
check_model(m2)  # Not too bad... some over-dispersion
check_collinearity(m2, component='all')  # Multicollinearity is fine
model_performance(m2)

# Model summary
summary(m2)


m2_null = update(m2, . ~ . - time:Treatment)
m2_null2 = update(m2, . ~ . - Treatment)


anova(m2_null, m2)  # Evidence of a Treatment:time interaction
anova(m2_null2, m2)  # Intercept does not depend on Treatment

m3 = update(m2, . ~ . - Treatment)

# Plot random effects term
tmp = REsim(m3)
plotREsim(tmp, labs=TRUE)

# Plot fixed effects
tmp2 = FEsim(m3)
plotFEsim(tmp2)

