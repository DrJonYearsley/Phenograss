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
                           cols = c(5:9),
                           names_to = c('cut','date'),
                           names_sep = '[[:punct:]]{2}',
                           values_to = 'biomass',
                           values_drop_na = TRUE,
                           names_pattern = "([0-9]{1}[[:punct:]]{2}[0-9]{2}[[:punct:]]{1}[0-9]{2}[[:punct:]]{1}[0-9]{4})")

biomass_long = pivot_longer(biomass,
                            cols = c(5:9),
                            names_to = c('cut','date'),
                            names_sep = "[[:punct:]]{2}",
                            values_to = 'biomass',
                            values_drop_na = TRUE,
                            names_prefix = "(Biomass[[:punct:]]{1}Cut[[:punct:]]{1})")


# Calculate days from 16th Dec 2019
biomass_long$date = as.Date(biomass_long$date, format="%d.%m.%Y.")
biomass_long$days = julian(biomass_long$date, origin=as.Date('2019-12-16'))


# Recode variables
biomass_long$cut = as.factor(biomass_long$cut)
biomass_long$Variety = as.factor(biomass_long$Variety)
biomass_long$Treatment = as.factor(biomass_long$Treatment)
biomass_long$Chamber = as.factor(biomass_long$Chamber)

#tmp2 = aggregate(days~Plant.ID+Variety+Treatment+cut, data=biomass_long, FUN=min)



# ******************************************************
# Visualise the data
format = theme(axis.title = element_text(size=14),
               axis.text = element_text(size=12))


ggplot(data=subset(biomass_long, Treatment=='Ambient' & Variety=='Aberchoice'),
       aes(x=days,
           y=biomass)) +
  geom_point(aes(colour=cut)) + 
  geom_path()


table(biomass_long$days)

# Total biomass per variety
tmp = aggregate(biomass~Plant.ID+Treatment+Variety, data=biomass_long, FUN=sum)
d1 = aggregate(biomass~Variety+Treatment, data=tmp, FUN=mean)
d2 = aggregate(biomass~Variety+Treatment, data=tmp, FUN=function(x){sd(x)/sqrt(sum(is.finite(x)))})

names(d1) = c('Variety','Treatment','mean')
names(d2) = c('Variety','Treatment','se')
d = merge(d1,d2)

varietyList = sample(unique(d$Variety), size=8)

ggplot(data=subset(d, Variety%in% varietyList),
       aes(x=Treatment,
           y=mean,
           ymin=mean-2*se,
           ymax=mean+2*se,
           colour=Variety,
           group=Variety)) +
  geom_point(size=3) +
  geom_path() +
  labs(y='Mean Biomass [grams]') +
  scale_color_brewer('Variety',palette='Dark2') + 
  theme_bw() + 
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        legend.text = element_text(size=18),
        legend.title = element_text(size=20))

ggsave(filename='biomass_reaction_norm.png', width=8, height=4)

# # Average of varieties and plants
# d = aggregate(biomass~days+Treatment+cut, 
#               data=biomass_long, 
#               FUN=median)
# 
# ggplot(data=d, 
#        aes(x=days, y=biomass, colour=Treatment, fill=cut, group=Treatment)) +
#   geom_path(size=1) +
#   geom_point(shape=21, size=3) +
#   scale_fill_brewer(name='Cuts', palette='Dark2') +
#   scale_colour_brewer(name='Treatment', palette='Dark2') +
#   labs(title='Average over plant ID') +
#   theme_bw() + 
#   format
# 




# *****************************
# Look at some initial analyses

library(DHARMa)

library(lme4)
library(performance)
library(merTools)




# Try fitting a linear model first
m_lm = lm(biomass~Variety * Treatment,
        data=tmp)

sim = simulateResiduals(m_lm, n=200)
plot(sim)

# Remove three interaction
m_lm_null1 = update(m_lm, .~.-Treatment:Variety)

# Remove Variety
m_lm_null2 = update(m_lm, .~.-Variety)

anova(m_lm_null1,m_lm, test='F')  # Significant ineraction interaction
anova(m_lm_null2,m_lm, test='F')  # Significant effect of Variety

summary(m_lm)

# *******************************
# Include Variety
# as a random effect on intercept
m2 = lmer(biomass~Treatment +(1|Variety), 
          data=tmp,
          REML=FALSE)

# Validate model (DHARMa)
sim_lme = simulateResiduals(m2, n=200)
plot(sim_lme)

#performance
check_model(m2)  # Not too bad... some mild over-dispersion perhaps
model_performance(m2)

# Model summary
summary(m2)


m2_null = update(m2, . ~ . - Treatment)


anova(m2_null, m2)  # Evidence of a Treatment effect


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

