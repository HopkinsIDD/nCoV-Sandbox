---
title: "nCoV 2019 Sandbox"
output:
  html_document:
    df_print: paged
---



# What is this?

The nCoV Sandbox is a running analytic blog we are "writing" as we try to apply some methods we had in the very early stages of development, and some old friends, to the 2019 nCoV outbreak. It also is us trying to run some analyses to get our on handle on, and keep up to date on, the epidmeiology of the emerging epidemic.

This is a bit of an excercise in radical transparency, and things are going start out very messy...but will hopefully get cleaner and more meaningful as things go. But the old stuff will (for the moment) remain at the bottom for posterity.

# Analytic Blog

## Reconstructing Past Incidence Using Cumulative Reports (1-28-2020)

Goal is to do a better job of recreating the daily case counts
in each area so we have implied epidemic curves to work with for
some of the more sophisticated stuff (hopefully) to come. 

First let's load in the data. Currently using only
confirmed cases (driven a bit by data source),
but unclear how long this will be viable.

<img src="nCoV-Sandbox_files/figure-html/unnamed-chunk-1-1.png" width="672" /><img src="nCoV-Sandbox_files/figure-html/unnamed-chunk-1-2.png" width="672" />

Things looks a little funny prior to the first, but this
does seem like it should give a rough pseudo epidemic curve for
the purpose of anlaysis.

## Basic Epi Summary 1-27-2020

First goal for the day, dig in deeper on the age specific data
and compare with the MERS-CoV data in a bit more detail.

First as always, load and sumarize the most recent Kudos line 
list (https://docs.google.com/spreadsheets/d/1jS24DjSPVWa4iuxuD4OAXrE3QeI8c9BC1hSlqr-NMiU/edit#gid=1187587451)


```
## Warning: The following named parsers don't match the column names: date
```

```
## Warning: Removed 31 rows containing non-finite values (stat_count).
```

<img src="nCoV-Sandbox_files/figure-html/unnamed-chunk-2-1.png" width="672" />
Note that we don't have any linelist information on the deaths
that occured before arou 1/15 in this line lisat. Moving forward with this data comparing with MERS-CoV data from Saudi Arabia 
through summer 2014.

<img src="nCoV-Sandbox_files/figure-html/unnamed-chunk-3-1.png" width="672" />

**Figure:** Odds ratio of death by age group for MERS=CoV and nCoV-2019. Log-scale.


**Table:** Odds ratio of death by age group for MERS=CoV and nCoV-2019

age_cat   nCoV                  MERS              
--------  --------------------  ------------------
0-9       -                     0.41 (0.11, 1.26) 
10-19     -                     0.19 (0.05, 0.52) 
20-29     -                     0.22 (0.12, 0.41) 
30-39     0.14 (0.01, 1.02)     0.20 (0.11, 0.35) 
40-49     0.16 (0.01, 1.19)     0.52 (0.31, 0.87) 
50-59     1                     1                 
60-69     5.88 (1.80, 23.39)    2.86 (1.59, 5.26) 
70+       17.71 (4.74, 82.15)   4.92 (2.79, 8.95) 

**Take aways from OR of death comparison**

- It looks like the pattern of relative mortality is similar
  in MERS-CoV and nCoV-2019 even if absolute rates are different.
- The paucity of data on age specific deaths in the current data 
  set for nCoV-2019 means uncertainty is huge
- Still assuming these are similar in a relative sense seems
  reasonable.

### Thought Experiment

What if nCoV symptomatic and death rates were identical to 
those of MERS-CoV. How many cases would the current line
list represent? How about the full data if they follow 
a similar age distribution?

Using mortality and infection rates for this paper 
in AJE on MERS-CoV symptomatic ratios and IFRs
ratios (10.1093/aje/kwv452), and a lot of assumptions:

1. Age distribution of all cases looks like line list cases.
2. Age distribution of deaths looks like line list deaths.
3. Confirmed cases (4,474) are roughly equal to symptomatic cases.
4. There are 107 deaths.
5. The symptomatic ratio is the same as MERS.
6. All line list cases that will die have died.

**Table:** Implied number of cases and needed ratio of IFR
in nCoV and MERS-CoV to reconcile deaths and implied cases.


Age        pr alive   pr dead   est. cases   est. dead   MERS symptomatic ratio   MERS IFR   Implied Infections by SR   Implied Infections by IFR   IFR Ratio to Reconcile
--------  ---------  --------  -----------  ----------  -----------------------  ---------  -------------------------  --------------------------  -----------------------
0-9            0.02      0.00        89.48        0.00                     0.11       0.10                     813.45                        0.00                     0.00
10-19          0.04      0.00       178.96        0.00                     0.11       0.05                    1626.91                        0.00                     0.00
20-29          0.10      0.00       425.03        0.00                     0.14       0.05                    3035.93                        0.00                     0.00
30-39          0.22      0.03      1006.65        2.74                     0.23       0.08                    4376.74                       34.29                     0.01
40-49          0.20      0.03       872.43        2.74                     0.39       0.17                    2237.00                       16.14                     0.01
50-59          0.14      0.10       648.73       10.97                     0.60       0.38                    1081.22                       28.88                     0.03
60-69          0.16      0.41       738.21       43.90                     0.78       0.63                     946.42                       69.68                     0.07
70+            0.12      0.44       514.51       46.64                     0.88       0.79                     584.67                       59.04                     0.10
Overall        1.00      1.00      4474.00      107.00                     0.46       0.31                   14702.34                      208.03                     0.02

So, if the symptomatic ratio for nCoV 2019 is similar to what was
implied by the confirmed cases of MERS-CoV (and other assumptions
hold) the following things are
true.:

- There are are at least 14,700 nCoV-2019 infections out there on
  27-1-2020. This is likely low as the 4,474 reported cases are
  likely a bit lower than there actually are.
- If this is the case, the IFR for nCoV is likely smaller than
  1/50th of that of MERS-CoV or lower 
  (so less than 6 deaths per 1,000 infections)
- The difference is bigger in younger individuals (1/100th or less)
  than older ones (1/10th).

**Note this is interesting note it is the result of a 
thought experiment only!!!**


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

<img src="nCoV-Sandbox_files/figure-html/unnamed-chunk-6-1.png" width="672" />

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
```

```
## Warning: All formats failed to parse. No formats found.
```

```r
  ##Filter to China:
  jhucsse_china <- jhucsse %>% 
    filter(Country_Region%in%c("Mainland China", "Macau", "Hong Kong"))

 
  
  jhucsse_china %>% drop_na(Confirmed) %>% 
    filter(Update>"2020-01-14") %>% 
  ggplot(aes(x=Update, y=Confirmed, col=Province_State)) +
    geom_line() + scale_y_log10()
```

<img src="nCoV-Sandbox_files/figure-html/unnamed-chunk-8-1.png" width="672" />

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
## Did you want `data = c(Country_Region, Update, Confirmed, Deaths, Recovered, Suspected, 
##     Demised)`?
```

```r
    kable(slopes, digits=2)
```



Province_State    slope   exp_scale
---------------  ------  ----------
Hubei              0.14        1.37
Zhejiang           0.29        1.96
Guangdong          0.24        1.75
Henan              0.61        4.07
Chongqing          0.46        2.89
Hunan              0.44        2.77
Anhui              0.41        2.58
Beijing            0.18        1.52
Sichuan            0.37        2.35
Shanghai           0.20        1.58
Shandong           0.41        2.55
Jiangxi            0.36        2.27
Guangxi            0.41        2.57
Jiangsu            0.40        2.49

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

<img src="nCoV-Sandbox_files/figure-html/unnamed-chunk-9-1.png" width="672" />

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

<img src="nCoV-Sandbox_files/figure-html/unnamed-chunk-10-1.png" width="672" />

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

<img src="nCoV-Sandbox_files/figure-html/unnamed-chunk-12-1.png" width="672" />
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
