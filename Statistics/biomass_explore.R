# Trial analysis/visualisation of biomass data
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# 15 Sept 2020
# ******************************************

rm(list=ls())

library(readxl)
library(tidyr)
library(ggplot2)

# download data and replace spaces in names with '.'
biomass = read_excel('./data_rosemount/Biomass Data.xlsx', 
                     .name_repair = 'universal')


# Put data in long format
biomass_long = pivot_longer(biomass,
                           cols = -c(1:4),
                           names_to = 'date',
                           values_to = 'biomass',
                           values_drop_na = TRUE,
                           names_pattern = "([0-9]{2}[[:punct:]]{1}[0-9]{2}[[:punct:]][0-9]{4})")

# Calculate days from 16th Dec 2019
biomass_long$date = as.Date(biomass_long$date, format="%d.%m.%Y")
biomass_long$days = julian(biomass_long$date, origin=as.Date('2019-12-16'))

# Biomass was cut week 5, 9 and 13
biomass_long$cut = 'Cut_1'

ind = difftime(biomass_long$date, as.Date('2020-02-13'), units='days')<=0
biomass_long$cut[ind] = 'Cut_0'

# ind = difftime(biomass_long$date, as.Date('2020-02-13'), units='days')>0 &
#   difftime(biomass_long$date, as.Date('2020-03-12'), units='days')<=0
# biomass_long$cut[ind] = 'Cut_1'

ind = difftime(biomass_long$date, as.Date('2020-03-12'), units='days')>0 &
  difftime(biomass_long$date, as.Date('2020-04-09'), units='days')<=0
biomass_long$cut[ind] = 'Cut_2'

ind = difftime(biomass_long$date, as.Date('2020-04-09'), units='days')>0
biomass_long$cut[ind] = 'Cut_3'

# Recode variables
biomass_long$cut = as.factor(biomass_long$cut)
biomass_long$Variety = as.factor(biomass_long$Variety)
biomass_long$Treatment = as.factor(biomass_long$Treatment)
biomass_long$Chamber = as.factor(biomass_long$Chamber)

tmp2 = aggregate(days~Plant.ID+Variety+Treatment+cut, data=biomass_long, FUN=min)



# ******************************************************
# Visualise the data
format = theme(axis.title = element_text(size=14),
               axis.text = element_text(size=12))



table(biomass_long$days)


# Average of varieties and plants
d = aggregate(biomass~days+Treatment+cut, data=biomass_long, FUN=median)

ggplot(data=d, 
       aes(x=days, y=biomass, colour=Treatment, fill=cut, group=Treatment)) +
  geom_path(size=1) +
  geom_point(shape=21, size=3) +
  scale_fill_brewer(name='Cuts', palette='Dark2') +
  scale_colour_brewer(name='Treatment', palette='Set1') +
  labs(title='Average over plant ID') +
  theme_bw() + 
  format





# *****************************
# Look at some initial analyses

library(DHARMa)

library(lme4)
library(performance)
library(merTools)


# Create a subset to remove cut-3 (not much data in cut-3)
gr_sub = subset(gr,cut!='Cut 3')


# Try fitting a linear model first
m_lm = glm(growth.rate~Variety + factor(time)*Treatment,
        data=gr,
        family=gaussian)

sim = simulateResiduals(m_lm, n=200)
plot(sim)

# Remove three interaction
m_lm_null1 = update(m_lm, .~.-factor(time)*Treatment)

# Remove Variety
m_lm_null2 = update(m_lm, .~.-Variety)

anova(m_lm_null1,m_lm, test='F')  # Significant ineraction interaction
anova(m_lm_null2,m_lm, test='F')  # Significant effect of Variety

summary(m_lm)

# *******************************
# Include Variety
# as a random effect on intercept
m2 = lmer(growth.rate~factor(time)*Treatment +(1|Variety), 
          data=gr_sub,
          REML=FALSE)

# Validate model (DHARMa)
sim_lme = simulateResiduals(m2, n=200)
plot(sim_lme)

#performance
check_model(m2)  # Not too bad... some mild over-dispersion perhaps
check_collinearity(m2, component='all')  # Multicollinearity is fine
model_performance(m2)

# Model summary
summary(m2)


m2_null = update(m2, . ~ . - factor(time):Treatment)
m2_null2 = update(m2, . ~ . - Treatment)


anova(m2_null, m2)  # Evidence of a Treatment:time interaction
anova(m2_null2, m2)  

# Refit model 
m3 = lmer(growth.rate~1 + factor(time)*Treatment  + (1|Variety), 
          data=gr_sub,
          REML=FALSE)

summary(m3)

# Plot random effects term
tmp = REsim(m3)
plotREsim(tmp, labs=TRUE)

# Plot fixed effects
tmp2 = FEsim(m3)
plotFEsim(tmp2)


# Function to check for overdispersion
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(m2)  # This test shows overdispersion is not an issue


# # Make prediction with model aggregating over Variety
gr_agg = aggregate(growth.rate~time+Treatment+cut, data=gr_sub, FUN=median)
# gr_agg$pred = predict(m3, newdata=gr_agg, re.form=NA, type='response')

# Use predictInterval for 95% CI then aggregate
tmp = predictInterval(m3, newdata=gr_sub, 
                      which='fixed', 
                      type='linear',
                      include.resid.var=0,
                      level=0.95)
gr_sub = cbind(gr_sub,tmp)
tmp_agg = aggregate(cbind(growth.rate,fit,lwr,upr)~time+Treatment+cut, data=gr_sub, FUN=median)




ggplot(data=tmp_agg, 
       aes(x=time, y=fit, ymax=upr, ymin=lwr,colour=Treatment, fill=cut, group=Treatment)) +
  geom_point(aes(y=growth.rate),shape=25, size=4) +
  geom_point(aes(y=fit),fill=NA,shape=21, size=2) +
  geom_linerange() +
  scale_fill_brewer(name='Cuts', palette='Dark2') +
  scale_colour_brewer(name='Treatment', palette='Set1') +
  labs(title='Average over plant ID & 95%CI') +
  theme_bw() + 
  format
