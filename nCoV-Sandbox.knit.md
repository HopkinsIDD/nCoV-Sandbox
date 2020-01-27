---
title: "nCoV 2019 Sandbox"
output:
  html_document:
    df_print: paged
---



# What is this?

The nCoV Sandbox is a running analytic blog we are "writing" as we try to apply some methods we had in the very early stages of development, and some old friends, to the 2019 nCoV outbreak. It also is us trying to run some analyses to get our on handle on, and keep up to date on, the epidmeiology of the emerging epidemic.

This is a bit of an excercise in radical transparency, and things are going start out very messy...but will hopefully get cleaner and more meaningful as things go. But the old stuff will (for the moment) remain at the botto for posterity.

# Analytic Blog

## Basic Epi Summary 1-25-2020

Three goals for today:

1. Make functions to automate most ofbasic epi report
2. Add a few new basic analyses.

### Summarizing the line list data
Age distribution and epicurve for cases where we have 
individual line list information. 


```
## Warning: Removed 18 rows containing non-finite values (stat_count).
```

<img src="nCoV-Sandbox_files/figure-html/unnamed-chunk-1-1.png" width="672" />

Now lets look at some basic infomration on survival by age group
and gender.


```
## Waiting for profiling to be done...
```



age_cat    alive   dead           OR  CI                       
--------  ------  -----  -----------  -------------------------
(0,10]         2      0    0.0000000  NA,2.4442929329126e+305  
(10,20]        4      0    0.0000000  NA,1.54362475336342e+123 
(20,30]       16      0    0.0000000  NA,5.68642814649887e+36  
(30,40]       25      1    0.1900000  0.01,1.41                
(40,50]       26      1    0.1826923  0.01,1.36                
(50,60]       19      4    1.0000000  -                        
(60,70]       10     16    7.6000000  2.15,32.49               
(70,80]        1      7   33.2500000  4.34,723.13              
(80,90]        1     10   47.5000000  6.56,1014.13             

```
## Waiting for profiling to be done...
```



gender    alive   dead          OR  CI        
-------  ------  -----  ----------  ----------
male         72     27   1.0000000  -         
female       35     12   0.9142857  0.58,1.99 

**Take aways from the line list data:**

- There is a huge survival effect by age.
- No apparent effect of gender.
- If the line list data is reflective of when deaths got
sick in the cumlative data, the overall  CFR may increase quite a bit in coming weeks.


### Bringing in the cumulative case data.

Now lets start to look at the aggregate cumulative case
data as that is going to be the most widely available, complete
and the basis for most of our predictive style analyses.

First we will focuse on Mainland China, Hong Kong and 
Macau.


```r
  jhucsse <- read_JHUCSSE_cases("2020-01-25 23:59", append_wiki = TRUE)
    
  ##Filter to China:
  jhucsse_china <- jhucsse %>% 
    filter(Country_Region%in%c("Mainland China", "Macau", "Hong Kong"))

 
  
  jhucsse_china %>% drop_na(Confirmed) %>% 
    filter(Update>"2020-01-14") %>% 
  ggplot(aes(x=Update, y=Confirmed, col=Province_State)) +
    geom_line() + scale_y_log10()
```

<img src="nCoV-Sandbox_files/figure-html/unnamed-chunk-3-1.png" width="672" />

Looking at all provinces, so let's narrow it to places that at some 
point experience at least 25 confimed cases and 
look vs. a straight log-linear line. 

Note that is is not quite right for real exponential growth since we 
are looking at the cumulative report rather than the 


```r
    tmp <- jhucsse_china%>%filter(Confirmed>=25)
    tmp <- unique(tmp$Province_State)
    
  
    ## Look at consitencey in exponential groqth by areas.
    analyze <-   jhucsse_china %>% drop_na(Confirmed) %>% 
      filter(Update>"2020-01-14") %>%
      filter(Province_State%in%tmp)
    
    #Get the slopes for each province. 
    slopes <- analyze %>% nest(-Province_State) %>%
      mutate(slope=map_dbl(data, ~lm(log10(.$Confirmed)~as.Date(.$Update))$coef[2])) %>%
      select(-data) %>% mutate(exp_scale=10^(slope))
```

```
## Warning: All elements of `...` must be named.
## Did you want `data = c(Country_Region, Update, Confirmed, Suspected, Recovered, Deaths, 
##     Demised)`?
```

```r
    kable(slopes, digits=2)
```



Province_State    slope   exp_scale
---------------  ------  ----------
Hubei              0.14        1.38
Zhejiang           0.33        2.12
Guangdong          0.19        1.57
Henan              0.39        2.46
Chongqing          0.38        2.38
Hunan              0.41        2.59
Anhui              0.56        3.59
Beijing            0.19        1.55
Sichuan            0.31        2.04
Shanghai           0.20        1.60
Shandong           0.43        2.67
Jiangxi            0.42        2.62
Guangxi            0.43        2.70
Jiangsu            0.47        2.97

```r
    #ggplot(slopes, aes(x=Province_State, y=slope)) +
    #         geom_bar(stat="identity") + coord_flip()
    
    ##Plot the exponential growth rate in eaach against a linear rate. 
    jhucsse_china %>% drop_na(Confirmed) %>% 
      filter(Update>"2020-01-14") %>%
      filter(Province_State%in%tmp)%>%
      ggplot(aes(x=Update, y=Confirmed, col=Province_State)) +
        geom_point() + scale_y_log10() + stat_smooth(method="lm", se=FALSE)
```

<img src="nCoV-Sandbox_files/figure-html/unnamed-chunk-4-1.png" width="672" />

Leaving it there for the moment due to lack of aggregate data. 

**Cumulative analysis preliminary so too early to say much but:**

- Growing everywhere with at least a few cases
- Rates seem very roughly similar



## Basic Epi Summary 1-24-2020

Simple snapshot as of 2020-24-1 based on snapshot of linelist data
derived from public sources from:
https://docs.google.com/spreadsheets/d/1jS24DjSPVWa4iuxuD4OAXrE3QeI8c9BC1hSlqr-NMiU/edit#gid=1449891965

(AKA the Kudos list).

This is some very basic episnapshots that should be improved 
in the coming days. 


First just take a rough look at the age distribution of cases.
Ten year increments.

```r
  source("R/DataLoadUtils.r")

  kudos <- readKudos2("data/Kudos Line List-1-24-2020.csv") %>%
   mutate(age_cat = cut(age, seq(0,100,10)))
  
  #Age distribution of cases.
  require(ggplot2)
  ggplot(drop_na(kudos, age_cat), 
         aes(x=age_cat, fill=as.factor(death))) + 
    geom_bar( color="grey") + coord_flip() + xlab("Age Catergory")
```

<img src="nCoV-Sandbox_files/figure-html/unnamed-chunk-5-1.png" width="672" />

Next, are we seeing any obvious differences in mortality
by gender or age?


```
## Waiting for profiling to be done...
```



gender    alive   dead     OR  CI        
-------  ------  -----  -----  ----------
male         54     15   1.00  -         
female       23      8   1.25  0.48,3.31 

```
## Waiting for profiling to be done...
```



age_cat    alive   dead       OR  CI                      
--------  ------  -----  -------  ------------------------
(10,20]        2      0     0.00  0,2.393735337629e+91    
(20,30]        8      0     0.00  NA,8.12332729629327e+59 
(30,40]       20      1     0.70  0.03,18.71              
(40,50]       21      1     0.67  0.02,17.8               
(50,60]       14      1     1.00  -                       
(60,70]        8      9    15.75  2.34,319.1              
(70,80]        1      3    42.00  2.83,1775.73            
(80,90]        1      8   112.00  9.58,4388.42            

Even as sparse as this data is, this is showing some clear
evidence of and age relationship. 


Epidemic curve of line list cases. Not
super informative at this point. 


```r
  ggplot(kudos, aes(x=symptom_onset, fill=as.factor(death))) +
  geom_bar()
```

```
## Warning: Removed 11 rows containing non-finite values (stat_count).
```

<img src="nCoV-Sandbox_files/figure-html/unnamed-chunk-7-1.png" width="672" />
A touch interesting that all deaths are early on. This suggests either (A) surveillance was really biased towards deaths in the early days, or (B) a lot of the later reports have not had time to die. 

[Note that there was perviously a 1-23-2020 summary 
but that was too preliminary even for this]

# Planning Notes/Ideas

- Apply basic framework to data so far, focusing on final size first
- Also do some basic epi summaries
    - here is a place for cool visulizations.
- Focus on province/state level in China (including Hong Kong/Macau SAR)
- Run this as a very open excercise in open science.
- Take snapshots of data...starting with stuff from : https://docs.google.com/spreadsheets/d/1jS24DjSPVWa4iuxuD4OAXrE3QeI8c9BC1hSlqr-NMiU/edit#gid=1449891965.
- Get total case snapshot from https://docs.google.com/spreadsheets/d/169AP3oaJZSMTquxtrkgFYMSp4gTApLTTWqo25qCpjL0/edit#gid=975486030
- Basically take snapshots and then post.
- Longer term...use approach and inference stuff for effect of actions.
- Keep a database of dates of intervention actions, major events.
    - 1/22/2020 : Wuhan Quarintine


### Some Tasks
- Download first snapshot, do some basic data cleaning and epi summaries.
- Do some data cleaning and loading
- Make some basic summaries
- Get province level covariates
    - population
    - population density
      - average/high/low
    - Average Feb temperature
    - Average/high-low Absolute humidity
    - Basic demographics (if available)
    - [OTHER WEATHER]
    - [OTHER INDEXES OF URBANIZATOIN/ECONOMY]

- Get the most basic model work
