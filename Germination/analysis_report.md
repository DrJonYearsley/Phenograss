Phenograss: Survival analysis of germination data
================
Jon Yearsley
26th March 2020

## Setup and load data

``` r
rm(list=ls())

library(survival)
library(survminer)
library(timereg)
library(broom)
library(pander)
library(emmeans)

# Load data
load('germination_data.RData')
```

## Overall visualisation

![](analysis_report_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

Look at overall numbers that germinated

``` r
germin_tab = table(germination$Treatment, 
                   germination$germinated)
pander(germin_tab, 
       caption = "Numbers germinating & not germinating within observation period")
```

|          | FALSE | TRUE |
| :------: | :---: | :--: |
| **CON**  |  440  | 820  |
| **eCO2** |  433  | 827  |

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
    ##  Aberchoice effect       0.6425 0.152 Inf  4.236  <.0001 
    ##  Abergain effect        -0.5276 0.130 Inf -4.068  0.0001 
    ##  Aspect effect           0.7662 0.157 Inf  4.895  <.0001 
    ##  Carraig effect          0.9855 0.167 Inf  5.917  <.0001 
    ##  Dunluce effect         -0.2578 0.131 Inf -1.965  0.0577 
    ##  Lilora effect          -0.1714 0.132 Inf -1.297  0.2094 
    ##  Moy effect             -1.6060 0.198 Inf -8.127  <.0001 
    ##  Semi-natural11 effect  -2.0883 0.224 Inf -9.318  <.0001 
    ##  Semi-natural6 effect    0.8450 0.221 Inf  3.832  0.0002 
    ##  Semi-natural7 effect   -1.3252 0.187 Inf -7.075  <.0001 
    ##  Solomon effect          0.9855 0.167 Inf  5.917  <.0001 
    ##  Wild4 effect            0.8450 0.221 Inf  3.832  0.0002 
    ##  Wild6 effect            0.0616 0.185 Inf  0.333  0.7393 
    ##  Wild7 effect            0.8450 0.221 Inf  3.832  0.0002 
    ## 
    ## Results are averaged over the levels of: Treatment 
    ## Results are given on the log odds ratio (not the response) scale. 
    ## P value adjustment: fdr method for 14 tests

Plot the posthoc analysis

![](analysis_report_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

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

|                                            | records | n.max | n.start | events | median |
| :----------------------------------------: | :-----: | :---: | :-----: | :----: | :----: |
|   **Variety=Aberchoice, Treatment=CON**    |   120   |  120  |   120   |   94   |   16   |
|   **Variety=Aberchoice, Treatment=eCO2**   |   120   |  120  |   120   |   92   |   13   |
|    **Variety=Abergain, Treatment=CON**     |   120   |  120  |   120   |   60   |   30   |
|    **Variety=Abergain, Treatment=eCO2**    |   120   |  120  |   120   |   64   |  22.5  |
|     **Variety=Aspect, Treatment=CON**      |   120   |  120  |   120   |  101   |   17   |
|     **Variety=Aspect, Treatment=eCO2**     |   120   |  120  |   120   |   90   |  14.5  |
|     **Variety=Carraig, Treatment=CON**     |   120   |  120  |   120   |  103   |   15   |
|    **Variety=Carraig, Treatment=eCO2**     |   120   |  120  |   120   |   96   |   10   |
|     **Variety=Dunluce, Treatment=CON**     |   120   |  120  |   120   |   68   |   27   |
|    **Variety=Dunluce, Treatment=eCO2**     |   120   |  120  |   120   |   72   |   22   |
|     **Variety=Lilora, Treatment=CON**      |   120   |  120  |   120   |   74   |   19   |
|     **Variety=Lilora, Treatment=eCO2**     |   120   |  120  |   120   |   71   |  18.5  |
|       **Variety=Moy, Treatment=CON**       |   60    |  60   |   60    |   14   |   NA   |
|      **Variety=Moy, Treatment=eCO2**       |   60    |  60   |   60    |   18   |   NA   |
| **Variety=Semi-natural11, Treatment=CON**  |   60    |  60   |   60    |   11   |   NA   |
| **Variety=Semi-natural11, Treatment=eCO2** |   60    |  60   |   60    |   11   |   NA   |
|  **Variety=Semi-natural6, Treatment=CON**  |   60    |  60   |   60    |   50   |   18   |
| **Variety=Semi-natural6, Treatment=eCO2**  |   60    |  60   |   60    |   47   |   14   |
|  **Variety=Semi-natural7, Treatment=CON**  |   60    |  60   |   60    |   16   |   NA   |
| **Variety=Semi-natural7, Treatment=eCO2**  |   60    |  60   |   60    |   23   |   NA   |
|     **Variety=Solomon, Treatment=CON**     |   120   |  120  |   120   |   99   |   17   |
|    **Variety=Solomon, Treatment=eCO2**     |   120   |  120  |   120   |  100   |   12   |
|      **Variety=Wild4, Treatment=CON**      |   60    |  60   |   60    |   46   |   17   |
|     **Variety=Wild4, Treatment=eCO2**      |   60    |  60   |   60    |   51   |   14   |
|      **Variety=Wild6, Treatment=CON**      |   60    |  60   |   60    |   39   |   22   |
|     **Variety=Wild6, Treatment=eCO2**      |   60    |  60   |   60    |   40   |   18   |
|      **Variety=Wild7, Treatment=CON**      |   60    |  60   |   60    |   45   |   19   |
|     **Variety=Wild7, Treatment=eCO2**      |   60    |  60   |   60    |   52   |   13   |

Table continues below

|                                            | 0.95LCL | 0.95UCL |
| :----------------------------------------: | :-----: | :-----: |
|   **Variety=Aberchoice, Treatment=CON**    |   15    |   17    |
|   **Variety=Aberchoice, Treatment=eCO2**   |   12    |   15    |
|    **Variety=Abergain, Treatment=CON**     |   25    |   NA    |
|    **Variety=Abergain, Treatment=eCO2**    |   18    |   NA    |
|     **Variety=Aspect, Treatment=CON**      |   17    |   18    |
|     **Variety=Aspect, Treatment=eCO2**     |   13    |   17    |
|     **Variety=Carraig, Treatment=CON**     |   14    |   15    |
|    **Variety=Carraig, Treatment=eCO2**     |   10    |   13    |
|     **Variety=Dunluce, Treatment=CON**     |   24    |   NA    |
|    **Variety=Dunluce, Treatment=eCO2**     |   19    |   29    |
|     **Variety=Lilora, Treatment=CON**      |   16    |   26    |
|     **Variety=Lilora, Treatment=eCO2**     |   15    |   NA    |
|       **Variety=Moy, Treatment=CON**       |   NA    |   NA    |
|      **Variety=Moy, Treatment=eCO2**       |   NA    |   NA    |
| **Variety=Semi-natural11, Treatment=CON**  |   NA    |   NA    |
| **Variety=Semi-natural11, Treatment=eCO2** |   NA    |   NA    |
|  **Variety=Semi-natural6, Treatment=CON**  |   17    |   19    |
| **Variety=Semi-natural6, Treatment=eCO2**  |   13    |   17    |
|  **Variety=Semi-natural7, Treatment=CON**  |   NA    |   NA    |
| **Variety=Semi-natural7, Treatment=eCO2**  |   NA    |   NA    |
|     **Variety=Solomon, Treatment=CON**     |   16    |   18    |
|    **Variety=Solomon, Treatment=eCO2**     |   11    |   13    |
|      **Variety=Wild4, Treatment=CON**      |   16    |   19    |
|     **Variety=Wild4, Treatment=eCO2**      |   13    |   15    |
|      **Variety=Wild6, Treatment=CON**      |   19    |   30    |
|     **Variety=Wild6, Treatment=eCO2**      |   17    |   26    |
|      **Variety=Wild7, Treatment=CON**      |   18    |   21    |
|     **Variety=Wild7, Treatment=eCO2**      |   13    |   16    |

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

|                    | records | n.max | n.start | events | median | 0.95LCL | 0.95UCL |
| :----------------: | :-----: | :---: | :-----: | :----: | :----: | :-----: | :-----: |
| **Treatment=CON**  |  1260   | 1260  |  1260   |  820   |   19   |   19    |   21    |
| **Treatment=eCO2** |  1260   | 1260  |  1260   |  827   |   17   |   16    |   17    |

Plot survival curves for Treatments (averaging across varieties)
![](analysis_report_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

#### Hypothesis test

Perform a log-ratio test of difference between survival curves for
control - treatment

``` r
pander(survdiff(Surv(Day, germinated) ~ Treatment, 
         data=germination))
```

|                    |  N   | Observed | Expected | (O-E)^2/E | (O-E)^2/V |
| :----------------: | :--: | :------: | :------: | :-------: | :-------: |
| **Treatment=CON**  | 1260 |   820    |  906.8   |   8.315   |   19.86   |
| **Treatment=eCO2** | 1260 |   827    |  740.2   |   10.19   |   19.86   |

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

![](analysis_report_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->![](analysis_report_files/figure-gfm/unnamed-chunk-19-2.png)<!-- -->

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

![](analysis_report_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

Same message. Looks like proportional hazard decreases over time for the
treatment effect

Need to find time varying version of Cox model.

### Effect of variety

Use Cox PH model on data from the control

``` r
# Fit the full model
coxph.fit_full <- coxph(Surv(Day, germinated) ~ Variety, 
                        data=subset(germination, Treatment=='CON'), 
                        method="breslow")  # Could use efron

# Fit null model testing for Variety on control data
coxph.fit_null <- coxph(Surv(Day, germinated) ~ 1, 
                         data=subset(germination, Treatment=='CON'), 
                         method="breslow")  # Could use efron

# Perform hypothesis test
anova(coxph.fit_full, coxph.fit_null, test='ChiSq')
```

    ## Analysis of Deviance Table
    ##  Cox model: response is  Surv(Day, germinated)
    ##  Model 1: ~ Variety
    ##  Model 2: ~ 1
    ##    loglik  Chisq Df P(>|Chi|)    
    ## 1 -5378.2                        
    ## 2 -5527.9 299.32 13 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Strong effect of Variety on germination. Plot the hazard ratio for the
model with main effect variety

    ## Warning: Removed 1 rows containing missing values (geom_errorbar).

![](analysis_report_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

``` r
survfit(coxph.fit_full)
```

    ## Call: survfit(formula = coxph.fit_full)
    ## 
    ##       n  events  median 0.95LCL 0.95UCL 
    ##    1260     820      21      19      22

## Time varying coefficients

### Modify Cox PH model (fit with time bins)

Bin data into 3 bins t\<14, 14\<t\<18, t\>18

``` r
# Create dataframe with time bins
germin2 <- survSplit(Surv(Day, germinated) ~ ., 
                  data= germination, 
                  cut=c(14,18),
                  episode= "tgroup", id="id")
```

``` r
# Fit model
coxph.fit_full2 <- coxph(Surv(Day, germinated) ~ Treatment:strata(tgroup), 
                        data=germin2, 
                        method="breslow")  # Could use efron
```

``` r
# Look at mean fitted values
coxph.fit_full2$means
```

    ##  TreatmentCON:strata(tgroup)tgroup=1 TreatmenteCO2:strata(tgroup)tgroup=1 
    ##                            0.2272727                            0.2272727 
    ##  TreatmentCON:strata(tgroup)tgroup=2 TreatmenteCO2:strata(tgroup)tgroup=2 
    ##                            0.1856061                            0.1351010 
    ##  TreatmentCON:strata(tgroup)tgroup=3 TreatmenteCO2:strata(tgroup)tgroup=3 
    ##                            0.1237374                            0.1010101

### Fit model with time varying coefficients

The hazard rate for Treatment will be allowed to vary with time

``` r
m = timecox(Surv(Day, germinated) ~ const(Variety)+Treatment, 
            data=germination, 
            max.time=25, 
            residuals=TRUE)
```

Plot results, and add in null expectation of no variation with time.
![](analysis_report_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->![](analysis_report_files/figure-gfm/unnamed-chunk-29-2.png)<!-- -->

``` r
summary(m)
```

    ## Multiplicative Hazard Model 
    ## 
    ## Test for nonparametric terms 
    ## 
    ## Test for non-significant effects 
    ##               Supremum-test of significance p-value H_0: B(t)=0
    ## (Intercept)                           157.0                   0
    ## TreatmenteCO2                          43.7                   0
    ## 
    ## Test for time invariant effects 
    ##                     Kolmogorov-Smirnov test p-value H_0:constant effect
    ## (Intercept)                             110                           0
    ## TreatmenteCO2                           102                           0
    ##                       Cramer von Mises test p-value H_0:constant effect
    ## (Intercept)                          131000                           0
    ## TreatmenteCO2                        110000                           0
    ## 
    ## Parametric terms :     
    ##                               Coef.    SE Robust SE       z    P-val lower2.5%
    ## const(Variety)Abergain       -1.230 0.121     0.161  -7.650 2.00e-14   -1.4700
    ## const(Variety)Aspect         -0.101 0.105     0.157  -0.644 5.20e-01   -0.3070
    ## const(Variety)Carraig         0.388 0.104     0.145   2.670 7.49e-03    0.1840
    ## const(Variety)Dunluce        -0.202 0.120     0.312  -0.649 5.16e-01   -0.4370
    ## const(Variety)Lilora         -0.567 0.114     0.181  -3.130 1.75e-03   -0.7900
    ## const(Variety)Moy            -1.100 0.218     0.731  -1.500 1.33e-01   -1.5300
    ## const(Variety)Semi-natural11 -2.600 0.241     0.259 -10.000 0.00e+00   -3.0700
    ## const(Variety)Semi-natural6  -0.103 0.128     0.189  -0.545 5.86e-01   -0.3540
    ## const(Variety)Semi-natural7  -0.755 0.188     0.582  -1.300 1.94e-01   -1.1200
    ## const(Variety)Solomon         0.176 0.104     0.150   1.170 2.40e-01   -0.0278
    ## const(Variety)Wild4          -0.062 0.128     0.182  -0.341 7.33e-01   -0.3130
    ## const(Variety)Wild6          -0.809 0.140     0.195  -4.150 3.37e-05   -1.0800
    ## const(Variety)Wild7          -0.231 0.127     0.174  -1.330 1.83e-01   -0.4800
    ##                              upper97.5%
    ## const(Variety)Abergain          -0.9930
    ## const(Variety)Aspect             0.1050
    ## const(Variety)Carraig            0.5920
    ## const(Variety)Dunluce            0.0332
    ## const(Variety)Lilora            -0.3440
    ## const(Variety)Moy               -0.6730
    ## const(Variety)Semi-natural11    -2.1300
    ## const(Variety)Semi-natural6      0.1480
    ## const(Variety)Semi-natural7     -0.3870
    ## const(Variety)Solomon            0.3800
    ## const(Variety)Wild4              0.1890
    ## const(Variety)Wild6             -0.5350
    ## const(Variety)Wild7              0.0179
    ##    
    ##   Call: 
    ## timecox(formula = Surv(Day, germinated) ~ const(Variety) + Treatment, 
    ##     data = germination, max.time = 25, residuals = TRUE)
