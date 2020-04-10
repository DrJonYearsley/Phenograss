Phenograss: Survival analysis of germination data
================
Jon Yearsley
10 April, 2020

## Setup and load data

``` r
rm(list=ls())

library(survival)
library(survminer)
library(timereg)
library(broom)
library(pander)
library(emmeans)
library(DHARMa)

# Load data
load('germination_data.RData')
```

Order varieties by median germination time in the control

``` r
control = subset(germination, Treatment=='Ambient')
varietyOrder = aggregate(Day~Variety, data=control, FUN=median)

# Order Varieties
varietyOrder = varietyOrder[order(varietyOrder$Day), ]
germination$Variety = reorder(germination$Variety, 
                              match(germination$Variety,
                                    varietyOrder$Variety))
```

## Overall visualisation

### Germination time distributions

![](analysis_report_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

plot the treatment and varieties in order of the median germination day
under the control

![](analysis_report_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

Plot just the control data

![](analysis_report_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

![](analysis_report_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

-----

Plot just the treatment data

![](analysis_report_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

![](analysis_report_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

Look at overall numbers that germinated

``` r
germin_tab = table(germination$Treatment, 
                   germination$germinated)
pander(germin_tab, 
       caption = "Numbers germinating & not germinating within observation period")
```

|               | FALSE | TRUE |
| :-----------: | :---: | :--: |
| **+CO2+Temp** |  433  | 827  |
|  **Ambient**  |  440  | 820  |

Numbers germinating & not germinating within observation period

Test to see if there’s an overall difference

``` r
pander(chisq.test(germin_tab), 
       caption = "Chi squared test of association on germin_tab")
```

| Test statistic | df | P value |
| :------------: | :-: | :-----: |
|     0.0631     | 1  | 0.8017  |

Chi squared test of association on germin\_tab

No difference in overall germination rate due to treatment

### Linear modelling of germination probability

Fit logistic model and test for difference between species and
treatments

``` r
m=glm(germinated~Treatment*Variety, data=germination, family='binomial')
```

Do some model validation

``` r
summary(m)
```

    ## 
    ## Call:
    ## glm(formula = germinated ~ Treatment * Variety, family = "binomial", 
    ##     data = germination)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -2.0074  -0.9833   0.6039   0.7585   1.8420  
    ## 
    ## Coefficients:
    ##                                        Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                             1.38629    0.22822   6.074 1.24e-09 ***
    ## TreatmentAmbient                        0.41522    0.34730   1.196 0.231861    
    ## VarietyAberchoice                      -0.19671    0.31411  -0.626 0.531157    
    ## VarietyAspect                          -0.28768    0.31069  -0.926 0.354473    
    ## VarietySolomon                          0.22314    0.33479   0.667 0.505078    
    ## VarietyWild4                            0.34831    0.42755   0.815 0.415273    
    ## VarietySemi-natural6                   -0.10110    0.38766  -0.261 0.794260    
    ## VarietyLilora                          -1.01543    0.29424  -3.451 0.000558 ***
    ## VarietyWild7                            0.48551    0.44307   1.096 0.273178    
    ## VarietyWild6                           -0.69315    0.35649  -1.944 0.051850 .  
    ## VarietyDunluce                         -0.98083    0.29463  -3.329 0.000871 ***
    ## VarietyAbergain                        -1.25276    0.29252  -4.283 1.85e-05 ***
    ## VarietyMoy                             -2.23359    0.36256  -6.161 7.24e-10 ***
    ## VarietySemi-natural11                  -2.88022    0.40423  -7.125 1.04e-12 ***
    ## VarietySemi-natural7                   -1.86172    0.35013  -5.317 1.05e-07 ***
    ## TreatmentAmbient:VarietyAberchoice     -0.31961    0.46508  -0.687 0.491951    
    ## TreatmentAmbient:VarietyAspect          0.15685    0.47707   0.329 0.742325    
    ## TreatmentAmbient:VarietySolomon        -0.47406    0.48820  -0.971 0.331525    
    ## TreatmentAmbient:VarietyWild4          -0.96024    0.58694  -1.636 0.101840    
    ## TreatmentAmbient:VarietySemi-natural6  -0.09098    0.58208  -0.156 0.875793    
    ## TreatmentAmbient:VarietyLilora         -0.31066    0.43630  -0.712 0.476452    
    ## TreatmentAmbient:VarietyWild7          -1.18841    0.59476  -1.998 0.045701 *  
    ## TreatmentAmbient:VarietyWild6          -0.48933    0.51853  -0.944 0.345332    
    ## TreatmentAmbient:VarietyDunluce        -0.55242    0.43506  -1.270 0.204166    
    ## TreatmentAmbient:VarietyAbergain       -0.54875    0.43293  -1.268 0.204968    
    ## TreatmentAmbient:VarietyMoy            -0.75751    0.54143  -1.399 0.161788    
    ## TreatmentAmbient:VarietySemi-natural11 -0.41522    0.58588  -0.709 0.478499    
    ## TreatmentAmbient:VarietySemi-natural7  -0.95140    0.52569  -1.810 0.070324 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 3251.8  on 2519  degrees of freedom
    ## Residual deviance: 2816.7  on 2492  degrees of freedom
    ## AIC: 2872.7
    ## 
    ## Number of Fisher Scoring iterations: 4

The residual deviance divided by degrees of freedom is about 1, so an
indication that there’s no overdispersion.

``` r
# Simulate residuals from the model
val = simulateResiduals(m, n=500, refit=FALSE)
```

Create some validation plots

``` r
plot(val)
```

![](analysis_report_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
testResiduals(val)
```

![](analysis_report_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

    ## $uniformity
    ## 
    ##  One-sample Kolmogorov-Smirnov test
    ## 
    ## data:  simulationOutput$scaledResiduals
    ## D = 0.014857, p-value = 0.6342
    ## alternative hypothesis: two-sided
    ## 
    ## 
    ## $dispersion
    ## 
    ##  DHARMa nonparametric dispersion test via sd of residuals fitted vs.
    ##  simulated
    ## 
    ## data:  simulationOutput
    ## ratioObsSim = 1.0006, p-value = 0.94
    ## alternative hypothesis: two.sided
    ## 
    ## 
    ## $outliers
    ## 
    ##  DHARMa outlier test based on exact binomial test
    ## 
    ## data:  simulationOutput
    ## outLow = 5.000e+00, outHigh = 9.000e+00, nobs = 2.520e+03, freqH0 =
    ## 1.996e-03, p-value = 0.1397
    ## alternative hypothesis: two.sided

    ## $uniformity
    ## 
    ##  One-sample Kolmogorov-Smirnov test
    ## 
    ## data:  simulationOutput$scaledResiduals
    ## D = 0.014857, p-value = 0.6342
    ## alternative hypothesis: two-sided
    ## 
    ## 
    ## $dispersion
    ## 
    ##  DHARMa nonparametric dispersion test via sd of residuals fitted vs.
    ##  simulated
    ## 
    ## data:  simulationOutput
    ## ratioObsSim = 1.0006, p-value = 0.94
    ## alternative hypothesis: two.sided
    ## 
    ## 
    ## $outliers
    ## 
    ##  DHARMa outlier test based on exact binomial test
    ## 
    ## data:  simulationOutput
    ## outLow = 5.000e+00, outHigh = 9.000e+00, nobs = 2.520e+03, freqH0 =
    ## 1.996e-03, p-value = 0.1397
    ## alternative hypothesis: two.sided

Logistic model looks valid. Fit null models and do some hypothesis
testing.

``` r
m0 = update(m, .~.-Treatment:Variety)

pander(anova(m0,m, test='Chisq'))
```

| Resid. Df | Resid. Dev | Df | Deviance | Pr(\>Chi) |
| :-------: | :--------: | :-: | :------: | :-------: |
|   2505    |    2829    | NA |    NA    |    NA     |
|   2492    |    2817    | 13 |  12.39   |   0.496   |

Analysis of Deviance Table

No evidence of an interaction for total germination proportion

Test for effect of variety

``` r
m0b = update(m, .~.-Treatment:Variety - Variety)
pander(anova(m0b,m0, test='Chisq'))
```

| Resid. Df | Resid. Dev | Df | Deviance | Pr(\>Chi) |
| :-------: | :--------: | :-: | :------: | :-------: |
|   2518    |    3252    | NA |    NA    |    NA     |
|   2505    |    2829    | 13 |  422.6   | 3.683e-82 |

Analysis of Deviance Table

A strong effect of Variety

Test for effect of Treatement

``` r
m0c = update(m, .~.-Treatment:Variety - Treatment)
pander(anova(m0c,m0, test='Chisq'))
```

| Resid. Df | Resid. Dev | Df | Deviance | Pr(\>Chi) |
| :-------: | :--------: | :-: | :------: | :-------: |
|   2506    |    2829    | NA |    NA    |    NA     |
|   2505    |    2829    | 1  |  0.1032  |  0.7481   |

Analysis of Deviance Table

No effect of Treatment

Posthoc test for the effect of Variety and compare each

``` r
m_eff = emmeans(m0, specs = 'Variety')
m_posthoc = contrast(m_eff, method='eff')
m_posthoc
```

    ##  contrast              estimate    SE  df z.ratio p.value
    ##  Carraig effect          0.9855 0.167 Inf  5.917  <.0001 
    ##  Aberchoice effect       0.6425 0.152 Inf  4.236  <.0001 
    ##  Aspect effect           0.7662 0.157 Inf  4.895  <.0001 
    ##  Solomon effect          0.9855 0.167 Inf  5.917  <.0001 
    ##  Wild4 effect            0.8450 0.221 Inf  3.832  0.0002 
    ##  Semi-natural6 effect    0.8450 0.221 Inf  3.832  0.0002 
    ##  Lilora effect          -0.1714 0.132 Inf -1.297  0.2094 
    ##  Wild7 effect            0.8450 0.221 Inf  3.832  0.0002 
    ##  Wild6 effect            0.0616 0.185 Inf  0.333  0.7393 
    ##  Dunluce effect         -0.2578 0.131 Inf -1.965  0.0577 
    ##  Abergain effect        -0.5276 0.130 Inf -4.068  0.0001 
    ##  Moy effect             -1.6060 0.198 Inf -8.127  <.0001 
    ##  Semi-natural11 effect  -2.0883 0.224 Inf -9.318  <.0001 
    ##  Semi-natural7 effect   -1.3252 0.187 Inf -7.075  <.0001 
    ## 
    ## Results are averaged over the levels of: Treatment 
    ## Results are given on the log odds ratio (not the response) scale. 
    ## P value adjustment: fdr method for 14 tests

Plot the posthoc analysis

![](analysis_report_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

## Calculate Kaplan-Meier survival curve

This is a non-parametric approach, and doesn’t have the same assumptions
as Cox PH model. We use the approach to estimate the survival curve for
germination, and to test if Treatment has an effect (averaged across
varieties)

The survival curve normally corresponds to the probability of dying. In
our case it is the probability of germinating. The hazard rate has a
similar interpretation.

### Survival curve for each Variety and Treatment

``` r
# Effect of Treatment for each Variety
germ_fit = survfit(Surv(Day, germinated) ~ Variety + Treatment, 
                             data=germination)
```

Display median survival times (i.e. germination times) and confidence
intervals

``` r
# Print estimates of time to germination
pander(germ_fit)
```

|                                                 | records | n.max | n.start | events | median |
| :---------------------------------------------: | :-----: | :---: | :-----: | :----: | :----: |
|    **Variety=Carraig, Treatment=+CO2+Temp**     |   120   |  120  |   120   |   96   |   10   |
|     **Variety=Carraig, Treatment=Ambient**      |   120   |  120  |   120   |  103   |   15   |
|   **Variety=Aberchoice, Treatment=+CO2+Temp**   |   120   |  120  |   120   |   92   |   13   |
|    **Variety=Aberchoice, Treatment=Ambient**    |   120   |  120  |   120   |   94   |   16   |
|     **Variety=Aspect, Treatment=+CO2+Temp**     |   120   |  120  |   120   |   90   |  14.5  |
|      **Variety=Aspect, Treatment=Ambient**      |   120   |  120  |   120   |  101   |   17   |
|    **Variety=Solomon, Treatment=+CO2+Temp**     |   120   |  120  |   120   |  100   |   12   |
|     **Variety=Solomon, Treatment=Ambient**      |   120   |  120  |   120   |   99   |   17   |
|     **Variety=Wild4, Treatment=+CO2+Temp**      |   60    |  60   |   60    |   51   |   14   |
|      **Variety=Wild4, Treatment=Ambient**       |   60    |  60   |   60    |   46   |   17   |
| **Variety=Semi-natural6, Treatment=+CO2+Temp**  |   60    |  60   |   60    |   47   |   14   |
|  **Variety=Semi-natural6, Treatment=Ambient**   |   60    |  60   |   60    |   50   |   18   |
|     **Variety=Lilora, Treatment=+CO2+Temp**     |   120   |  120  |   120   |   71   |  18.5  |
|      **Variety=Lilora, Treatment=Ambient**      |   120   |  120  |   120   |   74   |   19   |
|     **Variety=Wild7, Treatment=+CO2+Temp**      |   60    |  60   |   60    |   52   |   13   |
|      **Variety=Wild7, Treatment=Ambient**       |   60    |  60   |   60    |   45   |   19   |
|     **Variety=Wild6, Treatment=+CO2+Temp**      |   60    |  60   |   60    |   40   |   18   |
|      **Variety=Wild6, Treatment=Ambient**       |   60    |  60   |   60    |   39   |   22   |
|    **Variety=Dunluce, Treatment=+CO2+Temp**     |   120   |  120  |   120   |   72   |   22   |
|     **Variety=Dunluce, Treatment=Ambient**      |   120   |  120  |   120   |   68   |   27   |
|    **Variety=Abergain, Treatment=+CO2+Temp**    |   120   |  120  |   120   |   64   |  22.5  |
|     **Variety=Abergain, Treatment=Ambient**     |   120   |  120  |   120   |   60   |   30   |
|      **Variety=Moy, Treatment=+CO2+Temp**       |   60    |  60   |   60    |   18   |   NA   |
|       **Variety=Moy, Treatment=Ambient**        |   60    |  60   |   60    |   14   |   NA   |
| **Variety=Semi-natural11, Treatment=+CO2+Temp** |   60    |  60   |   60    |   11   |   NA   |
|  **Variety=Semi-natural11, Treatment=Ambient**  |   60    |  60   |   60    |   11   |   NA   |
| **Variety=Semi-natural7, Treatment=+CO2+Temp**  |   60    |  60   |   60    |   23   |   NA   |
|  **Variety=Semi-natural7, Treatment=Ambient**   |   60    |  60   |   60    |   16   |   NA   |

Table continues below

|                                                 | 0.95LCL | 0.95UCL |
| :---------------------------------------------: | :-----: | :-----: |
|    **Variety=Carraig, Treatment=+CO2+Temp**     |   10    |   13    |
|     **Variety=Carraig, Treatment=Ambient**      |   14    |   15    |
|   **Variety=Aberchoice, Treatment=+CO2+Temp**   |   12    |   15    |
|    **Variety=Aberchoice, Treatment=Ambient**    |   15    |   17    |
|     **Variety=Aspect, Treatment=+CO2+Temp**     |   13    |   17    |
|      **Variety=Aspect, Treatment=Ambient**      |   17    |   18    |
|    **Variety=Solomon, Treatment=+CO2+Temp**     |   11    |   13    |
|     **Variety=Solomon, Treatment=Ambient**      |   16    |   18    |
|     **Variety=Wild4, Treatment=+CO2+Temp**      |   13    |   15    |
|      **Variety=Wild4, Treatment=Ambient**       |   16    |   19    |
| **Variety=Semi-natural6, Treatment=+CO2+Temp**  |   13    |   17    |
|  **Variety=Semi-natural6, Treatment=Ambient**   |   17    |   19    |
|     **Variety=Lilora, Treatment=+CO2+Temp**     |   15    |   NA    |
|      **Variety=Lilora, Treatment=Ambient**      |   16    |   26    |
|     **Variety=Wild7, Treatment=+CO2+Temp**      |   13    |   16    |
|      **Variety=Wild7, Treatment=Ambient**       |   18    |   21    |
|     **Variety=Wild6, Treatment=+CO2+Temp**      |   17    |   26    |
|      **Variety=Wild6, Treatment=Ambient**       |   19    |   30    |
|    **Variety=Dunluce, Treatment=+CO2+Temp**     |   19    |   29    |
|     **Variety=Dunluce, Treatment=Ambient**      |   24    |   NA    |
|    **Variety=Abergain, Treatment=+CO2+Temp**    |   18    |   NA    |
|     **Variety=Abergain, Treatment=Ambient**     |   25    |   NA    |
|      **Variety=Moy, Treatment=+CO2+Temp**       |   NA    |   NA    |
|       **Variety=Moy, Treatment=Ambient**        |   NA    |   NA    |
| **Variety=Semi-natural11, Treatment=+CO2+Temp** |   NA    |   NA    |
|  **Variety=Semi-natural11, Treatment=Ambient**  |   NA    |   NA    |
| **Variety=Semi-natural7, Treatment=+CO2+Temp**  |   NA    |   NA    |
|  **Variety=Semi-natural7, Treatment=Ambient**   |   NA    |   NA    |

### Look at survival curves, averaging over varieties

``` r
# Effect of treatment averaging over variety
germ_fit_treatment = survfit(Surv(Day, germinated) ~ Treatment, 
                             data=germination)
```

Display

``` r
# Display results of average time to germination across varieties
pander(germ_fit_treatment)
```

|                         | records | n.max | n.start | events | median |
| :---------------------: | :-----: | :---: | :-----: | :----: | :----: |
| **Treatment=+CO2+Temp** |  1260   | 1260  |  1260   |  827   |   17   |
|  **Treatment=Ambient**  |  1260   | 1260  |  1260   |  820   |   19   |

Table continues below

|                         | 0.95LCL | 0.95UCL |
| :---------------------: | :-----: | :-----: |
| **Treatment=+CO2+Temp** |   16    |   17    |
|  **Treatment=Ambient**  |   19    |   21    |

Plot survival curves for Treatments (averaging across varieties)
![](analysis_report_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

#### Hypothesis test

Perform a log-ratio test of difference between survival curves for
control - treatment

``` r
pander(survdiff(Surv(Day, germinated) ~ Treatment, 
         data=germination))
```

|                         |  N   | Observed | Expected | (O-E)^2/E | (O-E)^2/V |
| :---------------------: | :--: | :------: | :------: | :-------: | :-------: |
| **Treatment=+CO2+Temp** | 1260 |   827    |  740.2   |   10.19   |   19.86   |
|  **Treatment=Ambient**  | 1260 |   820    |  906.8   |   8.315   |   19.86   |

Call: Surv(Day, germinated) \~ Treatment Chisq = 19.864297 on 1 degrees
of freedom, p = 0.000008

## Parametric models of germination

Try to fit a Cox proportional hazards model. This is a parametric model
that assumes a constant proportional hazard rate. This assumption must
be validated

Fit model with Treatment and Variety. There’s no evidence of an
interaction between Treatment and Variety

``` r
# Fit the full model with main effects and interaction
coxph.fit <- coxph(Surv(Day, germinated) ~ Treatment+Variety, 
                   data=germination, 
                   method="breslow")  # Could use efron
```

Same model but treat Variety as a random effect

``` r
# Fit the full model with main effects and random term for Variety
coxph.rand <- coxph(Surv(Day, germinated) ~ Treatment + frailty(Variety),
                    data=germination)
```

Validate the proportional hazard assumption

``` r
# Validate hazard function assumption
test.ph <- cox.zph(coxph.fit)
pander(test.ph$table)
```

|               | chisq | df |     p     |
| :-----------: | :---: | :-: | :-------: |
| **Treatment** |  222  | 1  | 3.368e-50 |
|  **Variety**  | 173.6 | 13 | 3.369e-30 |
|  **GLOBAL**   | 393.7 | 14 | 2.677e-75 |

Looks like the assumption is violated. Meaning that the effect of
Treatment and Variety on germination rate varies with time.

Look at some visuals

``` r
plot(test.ph)
```

![](analysis_report_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->![](analysis_report_files/figure-gfm/unnamed-chunk-31-2.png)<!-- -->

Looks like proportional hazard decreases over time for the treatment
effect. The effect of Variety seems to be less of an issue.

Look at model with Variety as a random effect

``` r
test2.ph <- cox.zph(coxph.rand)
pander(test2.ph$table)
```

|               | chisq |   df   |     p     |
| :-----------: | :---: | :----: | :-------: |
| **Treatment** | 236.9 | 0.9999 | 1.848e-53 |
|  **GLOBAL**   | 236.9 | 14.85  | 4.946e-42 |

``` r
plot(test2.ph)
```

![](analysis_report_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

Same message. Looks like proportional hazard decreases over time for the
treatment effect

Need to find time varying version of Cox model.
