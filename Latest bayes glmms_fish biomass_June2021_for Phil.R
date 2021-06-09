rm(list=ls())
options(max.print = 999999999)
setwd("C:/Users/Lrr19prl/OneDrive - Bangor University/COFUND/R")
fish <- read.csv("wsdbiotg.csv", header = TRUE,strip.white = T) #Uncapped data
names(fish)
head(fish)
str(fish)
summary(fish)

fish$SITE_SLOPE_100m <- as.numeric(fish$SITE_SLOPE_100m)
fish$SITE_SLOPE_200m <- as.numeric(fish$SITE_SLOPE_200m)
fish$SITE_SLOPE_400m <- as.numeric(fish$SITE_SLOPE_400m)
fish$DEPTH_BIN <- factor(fish$DEPTH_BIN, levels = c("Shallow", "Mid","Deep"))
fish$POP_STATUS <- factor(fish$POP_STATUS, levels = c("U","P"))
fish$SITE <- factor(fish$SITE)
fish$REGION <- factor(fish$REGION)
fish$REGION_NAME <- factor(fish$REGION_NAME)
fish$ISLAND <- factor(fish$ISLAND)
fish$OBS_YEAR <- factor(fish$OBS_YEAR, levels = c("2010", "2011","2012","2013","2014"))
is.factor(fish$DEPTH_BIN)
is.factor(fish$SITE)
is.factor(fish$OBS_YEAR)
str(fish)
summary(fish)

###################################################################
#Load packages and library files
library(dplyr)
library(lattice)  #Needed for multi-panel graphs
library(nlme)
library(lme4)
library(mgcv)
library(gamm4)
library(multcomp)
library(MASS)
library(MuMIn)
library(ggplot2)
library(effects)
library(car)
library(glmmTMB)
library(rstanarm)
library(rstan)
library(brms)
library(tidybayes)
library(bayesplot)
library(sjPlot)
library(loo)

library(here)
library(ggridges)
library(plyr)
library(dplyr)


#options help Stan run faster:
rstan_options(auto_write = TRUE)
#See how many cores the pc has: parallel::detectCores()
options(mc.cores = parallel::detectCores())
#Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
###################################################################

#Sort out weird things in the data:

#Remove islands: Gardner, Laysan, Maro, Midway, Necker, Nihoa, South Bank 
str(fish$ISLAND) #43 islands
fish <- fish[!(fish$ISLAND %in% c('Gardner', 'Laysan','Maro','Midway','Necker','Nihoa','South Bank')),] %>% droplevels
str(fish$ISLAND) #36 islands
levels(fish$ISLAND) 

#Remove weird sites: PAL-00107, GUA-00410, PHR-00341
str(fish$SITE) #2357 sites
fish <- fish[!(fish$SITE %in% c('PAL-00107', 'GUA-00410','PHR-00341')),] %>% droplevels
str(fish$SITE) #2354 sites
summary(fish)

#Also remove inner Maug sites because site slope data unreliable:
str(fish$SITE) #2354 sites
fish <- fish[!(fish$SITE %in% c('MAU-00077', 'MAU-00124', 'MAU-00162', 'MAU-00163', 'MAU-00174', 'MAU-00178', 'MAU-00179',
                                'MAU-00208', 'MAU-00218', 'MAU-00279','MAU-00281','MAU-00282','MAU-00284','MAU-00285','MAU-00294',
                                'MAU-00334','MAU-00339')),] %>% droplevels
str(fish$SITE) #2337 sites

#write.csv(fish,'SPC_sites.csv')

#Remove SPLE hammerheads from TotFish: (already removed from uncapped data)
#fish$TotFish1<-fish$TotFish-fish$SPLE

#Remove SPLE hammerheads from PISCIVORE: (already removed from uncapped data)
#fish$PISCIVORE1<-fish$PISCIVORE-fish$SPLE

#Subset zeros from biomass data for individual trophic group models:
PISC<-subset(fish,PISCIVORE>0) %>% droplevels
PLANK<-subset(fish,PLANKTIVORE>0) %>% droplevels
PRIM<-subset(fish,PRIMARY>0) %>% droplevels

POP<-subset(fish,POP_STATUS=="P") %>% droplevels
UNPOP<-subset(fish,POP_STATUS=="U") %>% droplevels
str(POP)

#Rescale variables:
fish$DEPTH_std <- fish$DEPTH / max(fish$DEPTH)
summary(fish$DEPTH_std)
fish$DEPTH_c <-scale(fish$DEPTH, center = TRUE, scale = FALSE)
fish$ISL_SLOPE_c<-scale(fish$ISL_SLOPE, center = TRUE, scale = FALSE)
fish$SITE_SLOPE_100m_c<-scale(fish$SITE_SLOPE_100m, center = TRUE, scale = FALSE)
fish$SITE_SLOPE_200m_c<-scale(fish$SITE_SLOPE_200m, center = TRUE, scale = FALSE)
fish$SITE_SLOPE_400m_c<-scale(fish$SITE_SLOPE_400m, center = TRUE, scale = FALSE)
summary(fish$DEPTH_c)
summary(fish$ISL_SLOPE_c)
summary(fish$SITE_SLOPE_100m_c)

PISC$DEPTH_c<-scale(PISC$DEPTH, center = TRUE, scale = FALSE)
PLANK$DEPTH_c<-scale(PLANK$DEPTH, center = TRUE, scale = FALSE)
PRIM$DEPTH_c<-scale(PRIM$DEPTH, center = TRUE, scale = FALSE)

PISC$ISL_SLOPE_c<-scale(PISC$ISL_SLOPE, center = TRUE, scale = FALSE)
PISC$SITE_SLOPE_100m_c<-scale(PISC$SITE_SLOPE_100m, center = TRUE, scale = FALSE)
PISC$SITE_SLOPE_200m_c<-scale(PISC$SITE_SLOPE_200m, center = TRUE, scale = FALSE)
PISC$SITE_SLOPE_400m_c<-scale(PISC$SITE_SLOPE_400m, center = TRUE, scale = FALSE)

PLANK$ISL_SLOPE_c<-scale(PLANK$ISL_SLOPE, center = TRUE, scale = FALSE)
PLANK$SITE_SLOPE_100m_c<-scale(PLANK$SITE_SLOPE_100m, center = TRUE, scale = FALSE)
PLANK$SITE_SLOPE_200m_c<-scale(PLANK$SITE_SLOPE_200m, center = TRUE, scale = FALSE)
PLANK$SITE_SLOPE_400m_c<-scale(PLANK$SITE_SLOPE_400m, center = TRUE, scale = FALSE)

PRIM$ISL_SLOPE_c<-scale(PRIM$ISL_SLOPE, center = TRUE, scale = FALSE)
PRIM$SITE_SLOPE_100m_c<-scale(PRIM$SITE_SLOPE_100m, center = TRUE, scale = FALSE)
PRIM$SITE_SLOPE_200m_c<-scale(PRIM$SITE_SLOPE_200m, center = TRUE, scale = FALSE)
PRIM$SITE_SLOPE_400m_c<-scale(PRIM$SITE_SLOPE_400m, center = TRUE, scale = FALSE)

###################################################################

#Some plots of the data:

ggplot(fish,aes(y=SITE_SLOPE_100m,x=SITE_SLOPE_200m))+geom_point()+geom_smooth(method = "lm")
ggplot(fish,aes(y=ISL_SLOPE,x=POP_STATUS))+geom_boxplot()

ggplot(fish,aes(y=SITE_SLOPE_100m,x=ISL_SLOPE))+geom_point()+geom_smooth(method = "lm")
ggplot(fish,aes(y=SITE_SLOPE_200m,x=ISL_SLOPE))+geom_point()+geom_smooth(method = "lm")
ggplot(fish,aes(y=SITE_SLOPE_400m,x=ISL_SLOPE))+geom_point()+geom_smooth(method = "lm")

mu100 <- ddply(fish, "ISLAND", summarise, grp.mean=mean(SITE_SLOPE_100m))
mu200 <- ddply(fish, "ISLAND", summarise, grp.mean=mean(SITE_SLOPE_200m))
mu400 <- ddply(fish, "ISLAND", summarise, grp.mean=mean(SITE_SLOPE_400m))

ggplot(fish, aes(x=SITE_SLOPE_100m)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")+
  geom_vline(data=mu100, aes(xintercept=grp.mean),linetype="dashed",colour="grey")+
  facet_wrap(~ISLAND)+
  theme_classic()
ggplot(fish, aes(x=SITE_SLOPE_200m)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")+
  geom_vline(data=mu200, aes(xintercept=grp.mean),linetype="dashed",colour="grey")+
  facet_wrap(~ISLAND)+
  theme_classic()
ggplot(fish, aes(x=SITE_SLOPE_400m)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")+
  geom_vline(data=mu400, aes(xintercept=grp.mean),linetype="dashed",colour="grey")+
  facet_wrap(~ISLAND)+
  theme_classic()

ggplot(fish,aes(y=TotFish,x=DEPTH,col=POP_STATUS))+geom_point()+geom_smooth(method='lm')+facet_grid(~OBS_YEAR)+theme_bw()
ggplot(fish,aes(y=TotFish,x=DEPTH,col=POP_STATUS))+geom_point()+geom_smooth(method='lm')+facet_grid(~REGION)+theme_bw()
ggplot(fish,aes(y=TotFish,x=DEPTH,col=POP_STATUS))+geom_point()+geom_smooth(method='lm')+facet_grid(~REGION_NAME)+theme_bw()
ggplot(fish,aes(y=TotFish,x=DEPTH,col=POP_STATUS))+geom_point()+geom_smooth(method='lm')+facet_wrap(~ISLAND)+theme_bw()

ggplot(fish,aes(y=PISCIVORE,x=DEPTH,fill=POP_STATUS))+geom_point()+geom_smooth(method='lm')+facet_grid(~REGION_NAME)+theme_bw()
ggplot(fish,aes(y=PLANKTIVORE,x=DEPTH,fill=POP_STATUS))+geom_point()+geom_smooth(method='lm')+facet_grid(~REGION_NAME)+theme_bw()
ggplot(fish,aes(y=PRIMARY,x=DEPTH,fill=POP_STATUS))+geom_point()+geom_smooth(method='lm')+facet_grid(~REGION_NAME)+theme_bw()
ggplot(fish,aes(y=SECONDARY,x=DEPTH,fill=POP_STATUS))+geom_point()+geom_smooth(method='lm')+facet_grid(~REGION_NAME)+theme_bw()

###################################################################

#Multicollinearity?

pairs( ~ TotFish + POP_STATUS + ISL_SLOPE + SITE_SLOPE_100m + SITE_SLOPE_200m + SITE_SLOPE_400m , data=fish, col = 'blue')
pairs( ~ TotFish + POP_STATUS + DEPTH_c + ISL_SLOPE_c + SITE_SLOPE_100m_c + SITE_SLOPE_200m_c + SITE_SLOPE_400m_c , 
       lower.panel=panel.cor, upper.panel=panel.smooth,
       data=fish)

#Look with vif?
library(lme4)
library(car)
m4<-glmer(TotFish ~ DEPTH_c +
            POP_STATUS +
            SITE_SLOPE_100m_c +
            SITE_SLOPE_400m_c +
            ISL_SLOPE_c +
            OBS_YEAR +
            (1|OBS_YEAR:REGION) + 
            (1|OBS_YEAR:ISLAND) + 
            (1|REGION)  + 
            (1|SITE) +
            (1+DEPTH_c||ISLAND),
          data=fish, family = Gamma(link=log),na.action = na.exclude)#,na.action=na.fail)
vif(m4) #All values < 3


###################################################################

###Total fish biomass (TotFish)

#Remove NAs in slope predictors for polys to work:
fish1<-fish[!is.na(fish$SITE_SLOPE_100m_c),]
fish1<-fish1[!is.na(fish1$SITE_SLOPE_400m_c),]
fish1<-fish1[!is.na(fish1$ISL_SLOPE_c),]
summary(fish1)

#1) Build the model:

psd <- rnorm(100000,2,2)
psd <- psd[psd>0]

# prior probability of multipliers/scale of variation between groups
mean(psd < log(100))
hist(psd, 100)

x<-seq(0,200,1) # make a sequence of the predictor 

y <- dlnorm(x, log(124.07), 0.44)

plot(x,y)

#Set priors:

#MacNeil et al. 2015 Nature: Resident fish biomass in absense of fishing averages 1,013 (963, 1469) kg ha-1 
#(posterior median (95% highest posterior density intervals)). So in g-m2: 101.3 (96.3, 146.9)
#So logged:
log(101.3) #4.618086

priors = c(set_prior('normal(0,2)', class='b'),
           set_prior('cauchy(0,5)', class='sd'),
           set_prior('cauchy(0,1)', class='sd', coef = 'DEPTH_c', group = 'ISLAND'),
           set_prior('normal(4.6,1)', class='Intercept'))

#Looking to see if priors are way out by drawing on the prior distribution with 'sample_prior = "only"':

fish1$POP_STATUS <- factor(fish1$POP_STATUS, levels = c('U', 'P'))

TotFish.Gamma.brms.prior.1 = brm(TotFish ~ DEPTH_c +
                                   POP_STATUS +
                                   poly(SITE_SLOPE_100m_c,2) +
                                   poly(SITE_SLOPE_400m_c,2) +
                                   poly(ISL_SLOPE_c,2) +
                                   DEPTH_c:POP_STATUS +
                                   DEPTH_c:poly(SITE_SLOPE_100m_c,2) +
                                   DEPTH_c:poly(SITE_SLOPE_400m_c,2) +
                                   DEPTH_c:poly(ISL_SLOPE_c,2) +
                                   OBS_YEAR +
                                   (1|OBS_YEAR:REGION) + 
                                   (1|OBS_YEAR:ISLAND) + 
                                   (1|REGION)  + 
                                   (1|SITE) +
                                   (1+DEPTH_c||ISLAND),
                                 data=fish1, 
                                 family = Gamma(link=log),
                                 prior = priors,
                                 chains = 4,
                                 cores  = 4,
                                 iter = 1000,
                                 warmup = 100,
                                 sample_prior = "only")
summary(TotFish.Gamma.brms.prior.1)
prior_summary(TotFish.Gamma.brms.prior.1)
hist(log(fitted(TotFish.Gamma.brms.prior.1)[, 1]))

#Fit model with chosen priors:

## NEW MAKE INIT LIST 
start_names <- unique(paste0('sd_',TotFish.Gamma.brms.prior.1$prior$group,'__',TotFish.Gamma.brms.prior.1$prior$coef)[TotFish.Gamma.brms.prior.1$prior$class=='sd'])
start_names <- start_names[!grepl('.__$', start_names)]

SLIST <- as.list(start_names)
names(SLIST) <- start_names
SLIST <- lapply(SLIST, function(x) 0.1)

inilist <- c(SLIST,
             list(b_Intercept=4, 
                  b_DEPTH_c=0,
                  `b_DEPTH_c:POP_STATUSP`=0,
                  `b_DEPTH_c:ISL_SLOPE_c`=0,
                  `b_DEPTH_c:SITE_SLOPE_100m_c`=0,
                  `b_DEPTH_c:SITE_SLOPE_400m_c`=0,
                  b_OBS_YEAR_2011=0,
                  b_OBS_YEAR_2012=0,
                  b_OBS_YEAR_2013=0,
                  b_OBS_YEAR_2014=0,
                  b_POP_STATUSP=0,
                  b_ISL_SLOPE_c=0,
                  b_SITE_SLOPE_100m_c=0,
                  b_SITE_SLOPE_400m_c=0,
                  shape=2.5))

TotFish.Gamma.brms.full.1 = brm(TotFish ~ DEPTH_c +
                                  POP_STATUS +
                                  poly(SITE_SLOPE_100m_c,2) +
                                  poly(SITE_SLOPE_400m_c,2) +
                                  poly(ISL_SLOPE_c,2) +
                                  DEPTH_c:POP_STATUS +
                                  DEPTH_c:poly(SITE_SLOPE_100m_c,2) +
                                  DEPTH_c:poly(SITE_SLOPE_400m_c,2) +
                                  DEPTH_c:poly(ISL_SLOPE_c,2) +
                                  OBS_YEAR +
                                  (1|OBS_YEAR:REGION) + 
                                  (1|OBS_YEAR:ISLAND) + 
                                  (1|REGION)  + 
                                  (1|SITE) +
                                  (1+DEPTH_c||ISLAND),
                                inits = function() inilist,
                                data=fish1, 
                                family = Gamma(link=log),
                                prior = priors,
                                chains = 4,
                                cores  = 4,
                                thin = 4,
                                iter = 2500,
                                warmup = 500, 
                                control = list(adapt_delta=0.95))


#2) Taking a look at the model:
mcmc_plot(TotFish.Gamma.brms.full.1, type='trace')
pp = brms::pp_check(TotFish.Gamma.brms.full.1, type = 'ecdf_overlay', nsamples=200) 
pp + theme_bw()+ xlim(c(0,500)) + xlab('Density (g.m2)') + ylab("Cumulative probability")
brms::pp_check(TotFish.Gamma.brms.full.1, type = 'loo_pit_qq', nsamples=200) 
brms::pp_check(TotFish.Gamma.brms.full.1, type = 'error_scatter_avg_vs_x', nsamples=200, x='DEPTH_c') 

summary(TotFish.Gamma.brms.full.1)
conditional_effects(TotFish.Gamma.brms.full.1)

#CHECK TO SEE IF YOU NEED RANDOM DEPTH; THEN SEE IF NEED RANDOM OBS_YEAR terms:
#1. Need random slope for depth?
TotFish.Gamma.brms.full.1 = brm(TotFish ~ DEPTH_c +
                                  POP_STATUS +
                                  poly(SITE_SLOPE_100m_c,2) +
                                  poly(SITE_SLOPE_400m_c,2) +
                                  poly(ISL_SLOPE_c,2) +
                                  DEPTH_c:POP_STATUS +
                                  DEPTH_c:poly(SITE_SLOPE_100m_c,2) +
                                  DEPTH_c:poly(SITE_SLOPE_400m_c,2) +
                                  DEPTH_c:poly(ISL_SLOPE_c,2) +
                                  OBS_YEAR +
                                  (1|OBS_YEAR:REGION) + 
                                  (1|OBS_YEAR:ISLAND) + 
                                  (1|REGION)  + 
                                  (1|SITE) +
                                  (1+DEPTH_c||ISLAND),
                                inits = function() inilist,
                                data=fish1, 
                                family = Gamma(link=log),
                                prior = priors,
                                chains = 4,
                                cores  = 4,
                                thin = 4,
                                iter = 1000,
                                warmup = 100, 
                                control = list(adapt_delta=0.95))

TotFish.Gamma.brms.full.1a = brm(TotFish ~ DEPTH_c +
                                   POP_STATUS +
                                   poly(SITE_SLOPE_100m_c,2) +
                                   poly(SITE_SLOPE_400m_c,2) +
                                   poly(ISL_SLOPE_c,2) +
                                   DEPTH_c:POP_STATUS +
                                   DEPTH_c:poly(SITE_SLOPE_100m_c,2) +
                                   DEPTH_c:poly(SITE_SLOPE_400m_c,2) +
                                   DEPTH_c:poly(ISL_SLOPE_c,2) +
                                   OBS_YEAR +
                                   (1|OBS_YEAR:REGION) + 
                                   (1|OBS_YEAR:ISLAND) + 
                                   (1|REGION)  + 
                                   (1|SITE) +
                                   (1|ISLAND),
                                 inits = function() inilist,
                                 data=fish1, 
                                 family = Gamma(link=log),
                                 prior = c(set_prior('normal(0,2)', class='b'),
                                           set_prior('cauchy(0,5)', class='sd'),
                                           set_prior('normal(4.6,1)', class='Intercept')),
                                 chains = 1,
                                 cores  = 1,
                                 thin = 4,
                                 iter = 1000,
                                 warmup = 100)#, 
#control = list(adapt_delta=0.95))

loofull<-loo(TotFish.Gamma.brms.full.1)
looA<-loo(TotFish.Gamma.brms.full.1a)

loo_compare(loofull,looA)

mcmc_plot(TotFish.Gamma.brms.full.1a, type='trace')
pp = brms::pp_check(TotFish.Gamma.brms.full.1a, type = 'ecdf_overlay', nsamples=200) 
pp + theme_bw()+ xlim(c(0,500)) + xlab('Density (g.m2)') + ylab("Cumulative probability")
brms::pp_check(TotFish.Gamma.brms.full.1a, type = 'loo_pit_qq', nsamples=200) 

#2. Need random OBS_YEAR terms?
TotFish.Gamma.brms.full.1b = brm(TotFish ~ DEPTH_c +
                                   POP_STATUS +
                                   poly(SITE_SLOPE_100m_c,2) +
                                   poly(SITE_SLOPE_400m_c,2) +
                                   poly(ISL_SLOPE_c,2) +
                                   DEPTH_c:POP_STATUS +
                                   DEPTH_c:poly(SITE_SLOPE_100m_c,2) +
                                   DEPTH_c:poly(SITE_SLOPE_400m_c,2) +
                                   DEPTH_c:poly(ISL_SLOPE_c,2) +
                                   OBS_YEAR +
                                   (1|OBS_YEAR:ISLAND) + 
                                   (1|REGION)  + 
                                   (1|SITE) +
                                   (1+DEPTH_c||ISLAND),
                                 inits = function() inilist,
                                 data=fish1, 
                                 family = Gamma(link=log),
                                 prior = priors,
                                 chains = 1,
                                 cores  = 1,
                                 thin = 4,
                                 iter = 1000,
                                 warmup = 100)#, 
#control = list(adapt_delta=0.95))

TotFish.Gamma.brms.full.1c = brm(TotFish ~ DEPTH_c +
                                   POP_STATUS +
                                   poly(SITE_SLOPE_100m_c,2) +
                                   poly(SITE_SLOPE_400m_c,2) +
                                   poly(ISL_SLOPE_c,2) +
                                   DEPTH_c:POP_STATUS +
                                   DEPTH_c:poly(SITE_SLOPE_100m_c,2) +
                                   DEPTH_c:poly(SITE_SLOPE_400m_c,2) +
                                   DEPTH_c:poly(ISL_SLOPE_c,2) +
                                   OBS_YEAR +
                                   (1|OBS_YEAR:REGION) + 
                                   (1|REGION)  + 
                                   (1|SITE) +
                                   (1+DEPTH_c||ISLAND),
                                 inits = function() inilist,
                                 data=fish1, 
                                 family = Gamma(link=log),
                                 prior = priors,
                                 chains = 1,
                                 cores  = 1,
                                 thin = 4,
                                 iter = 1000,
                                 warmup = 100)#, 
#control = list(adapt_delta=0.95))

TotFish.Gamma.brms.full.1d = brm(TotFish ~ DEPTH_c +
                                   POP_STATUS +
                                   poly(SITE_SLOPE_100m_c,2) +
                                   poly(SITE_SLOPE_400m_c,2) +
                                   poly(ISL_SLOPE_c,2) +
                                   DEPTH_c:POP_STATUS +
                                   DEPTH_c:poly(SITE_SLOPE_100m_c,2) +
                                   DEPTH_c:poly(SITE_SLOPE_400m_c,2) +
                                   DEPTH_c:poly(ISL_SLOPE_c,2) +
                                   OBS_YEAR +
                                   (1|REGION)  + 
                                   (1|SITE) +
                                   (1+DEPTH_c||ISLAND),
                                 inits = function() inilist,
                                 data=fish1, 
                                 family = Gamma(link=log),
                                 prior = priors,
                                 chains = 1,
                                 cores  = 1,
                                 thin = 4,
                                 iter = 1000,
                                 warmup = 100)#, 
#control = list(adapt_delta=0.95))

TotFish.Gamma.brms.full.1e = brm(TotFish ~ DEPTH_c +
                                   POP_STATUS +
                                   poly(SITE_SLOPE_100m_c,2) +
                                   poly(SITE_SLOPE_400m_c,2) +
                                   poly(ISL_SLOPE_c,2) +
                                   DEPTH_c:POP_STATUS +
                                   DEPTH_c:poly(SITE_SLOPE_100m_c,2) +
                                   DEPTH_c:poly(SITE_SLOPE_400m_c,2) +
                                   DEPTH_c:poly(ISL_SLOPE_c,2) +
                                   OBS_YEAR +
                                   (1|REGION)  + 
                                   (1|SITE) +
                                   (1|ISLAND),
                                 inits = function() inilist,
                                 data=fish1, 
                                 family = Gamma(link=log),
                                 prior = c(set_prior('normal(0,2)', class='b'),
                                           set_prior('cauchy(0,5)', class='sd'),
                                           set_prior('normal(4.6,1)', class='Intercept')),
                                 chains = 1,
                                 cores  = 1,
                                 thin = 4,
                                 iter = 1000,
                                 warmup = 100)#, 
#control = list(adapt_delta=0.95))

looB<-loo(TotFish.Gamma.brms.full.1b)
looC<-loo(TotFish.Gamma.brms.full.1c)
looD<-loo(TotFish.Gamma.brms.full.1d)
looE<-loo(TotFish.Gamma.brms.full.1e)

loo_compare(looB,looC,looD,looE)

TotFish.Gamma.brms.full.1d = brm(TotFish ~ DEPTH_c +
                                   POP_STATUS +
                                   poly(SITE_SLOPE_100m_c,2) +
                                   poly(SITE_SLOPE_400m_c,2) +
                                   poly(ISL_SLOPE_c,2) +
                                   DEPTH_c:POP_STATUS +
                                   DEPTH_c:poly(SITE_SLOPE_100m_c,2) +
                                   DEPTH_c:poly(SITE_SLOPE_400m_c,2) +
                                   DEPTH_c:poly(ISL_SLOPE_c,2) +
                                   OBS_YEAR +
                                   (1|REGION)  + 
                                   (1|SITE) +
                                   (1+DEPTH_c||ISLAND),
                                 inits = function() inilist,
                                 data=fish1, 
                                 family = Gamma(link=log),
                                 prior = priors,
                                 chains = 4,
                                 cores  = 4,
                                 thin = 4,
                                 iter = 2500,
                                 warmup = 500,
                                 control = list(adapt_delta=0.95))

#Taking a look at the model:
mcmc_plot(TotFish.Gamma.brms.full.1d, type='trace')
pp = brms::pp_check(TotFish.Gamma.brms.full.1d, type = 'ecdf_overlay', nsamples=200) 
pp + theme_bw()+ xlim(c(0,500)) + xlab('Density (g.m2)') + ylab("Cumulative probability")
brms::pp_check(TotFish.Gamma.brms.full.1d, type = 'loo_pit_qq', nsamples=200) 
brms::pp_check(TotFish.Gamma.brms.full.1d, type = 'error_scatter_avg_vs_x', nsamples=200, x='DEPTH_c') 

summary(TotFish.Gamma.brms.full.1d)
conditional_effects(TotFish.Gamma.brms.full.1d)

#Look at variance of random spatial terms:

mcmc_plot(TotFish.Gamma.brms.full.1d, pars=c('sd_SITE__Intercept','sd_ISLAND__Intercept','sd_REGION__Intercept'), fixed = T) + 
  #xlim(0,2.5) +
  xlab("TotFish model Standard deviation") +
  theme_bw()

mcmc_plot(TotFish.Gamma.brms.full.1d, pars='sd') + 
  #xlim(0,2.5) +
  xlab("TotFish model Standard deviation") +
  theme_bw()

sds <- posterior_samples(TotFish.Gamma.brms.full, pars = 'sd')
hist(sds$sd_OBS_YEAR:REGION__Intercept/sds$sd_OBS_YEAR:ISLAND__Intercept)

pairs(TotFish.Gamma.brms.full.1, pars=c('sd_SITE__Intercept','sd_ISLAND__Intercept','sd_REGION__Intercept'), fixed = T)

sds <- posterior_samples(TotFish.Gamma.brms.full, pars = 'sd')

hist(sds$sd_REGION__Intercept/sds$sd_ISLAND__Intercept)
hist(sds$sd_REGION__Intercept/sds$sd_SITE__Intercept)
hist(sds$sd_ISLAND__Intercept/sds$sd_SITE__Intercept)
#What is the probability that the scale of variation in regions is greater
#than the scale of variation among islands?:
mean(sds$sd_REGION__Intercept > sds$sd_ISLAND__Intercept)
#Is approx 50:50 so might say that "There is no support to say that either of these scales dominates" or if it was more like 0.05
#might make a statement like "there is strong support that the scale of variation for region is not important in driving biomass"
mean(sds$sd_REGION__Intercept > sds$sd_SITE__Intercept)
mean(sds$sd_ISLAND__Intercept > sds$sd_SITE__Intercept)

- sds$sd_REGION__Intercept - sds$sd_ISLAND__Intercept

rel_imp <- apply(sds,1,function(x) x/sum(x))

bayes_R2(TotFish.Gamma.brms.full.1)
performance::r2_bayes(TotFish.Gamma.brms.full.1) #marginal = fixed effs; conditional is random and fixed effects together 
#(so ranef r2 = 0.486%)

#3) Query the model: **WILL DO THIS LATER**
#Probability of increase in X gm2 per 5 (?) m (across 0-30m depth; at unpop islands; with slope held constant?)? And at pop islands?

#Notes from Murray course: 
#Predictive capactiy of these analyses: 
#Can ask, "what proportion of slopes are > 0.808?" must be 50% because we used a normal (gaussian) distribution (centred equally around mean)
#hist(as.matrix(fert.rstanarm)[,2])

#Q. What is probability that rate of change is more than 0.6?
#sum(as.matrix(fert.rstanarm)[,2]>0.6)/4000
#mean(as.matrix(model)[,2]>0.6)


##################################################################################################################
#Trophic group proportions

#PISC:
#Prop: 0.2 mean, 0.5 sd?
hist(rbeta(10000,10,40))
hist(rbeta(10000,2,8))
mean(log(rlnorm(10000, 4.6, 0.21)*rbeta(10000,10,40)))
sd(log(rlnorm(10000, 4.6, 0.21)*rbeta(10000,10,40)))
#2.95,0.35

#PLANK:
#Prop: 0.3 mean, 0.4 sd?
hist(rbeta(10000,15,35))
hist(rbeta(10000,3,7))
mean(log(rlnorm(10000, 4.6, 0.21)*rbeta(10000,15,35)))
sd(log(rlnorm(10000, 4.6, 0.21)*rbeta(10000,15,35)))
#3.36,0.30

#PRIM: 
#(a+b+c+d in MacNeil et al. Extended Data Figure 6: Average reef fish functional group across a biomass gradient)
#Prop: 0.8 mean?
hist(rbeta(10000,40,10))
hist(rbeta(10000,8,2))
mean(log(rlnorm(10000, 4.6, 0.21)*rbeta(10000,40,10)))
sd(log(rlnorm(10000, 4.6, 0.21)*rbeta(10000,40,10)))
#4.37, 0.22

#SEC:
#(f+g+h in MacNeil et al. Extended Data Figure 6: Average reef fish functional group across a biomass gradient)
#Prop: 0.65 mean?
hist(rbeta(10000,32.5,17.5))
hist(rbeta(10000,6.5,3.5))
mean(log(rlnorm(10000, 4.6, 0.21)*rbeta(10000,32.5,17.5)))
sd(log(rlnorm(10000, 4.6, 0.21)*rbeta(10000,32.5,17.5)))
#4.16,0.23

##################################################################################################################

#PISCIVORE

#Remove NAs in slope predictors for polys to work:
PISC1<-PISC[!is.na(PISC$SITE_SLOPE_100m_c),]
PISC1<-PISC1[!is.na(PISC1$SITE_SLOPE_400m_c),]
PISC1<-PISC1[!is.na(PISC1$ISL_SLOPE_c),]
summary(PISC1)

#1) Build the model:

#Trophic group proportions from MacNeil et al. (2015)
#PISC:
#Prop: 0.2 mean, 0.5 sd?
hist(rbeta(10000,10,40))
hist(rbeta(10000,2,8))
mean(log(rlnorm(10000, 4.6, 0.21)*rbeta(10000,10,40)))
sd(log(rlnorm(10000, 4.6, 0.21)*rbeta(10000,10,40)))
#2.95,0.35

#Set priors:
priors_PISC = c(set_prior('normal(0,2)', class='b'),
                set_prior('cauchy(0,5)', class='sd'),
                set_prior('cauchy(0,1)', class='sd', coef = 'DEPTH_c', group = 'ISLAND'),
                set_prior('normal(2.95,0.35)', class='Intercept'))

PISC1$POP_STATUS <- factor(PISC1$POP_STATUS, levels = c('U', 'P'))

PISC.Gamma.brms.prior = brm(PISCIVORE ~ DEPTH_c +
                              POP_STATUS +
                              poly(SITE_SLOPE_100m_c,2) +
                              poly(SITE_SLOPE_400m_c,2) +
                              poly(ISL_SLOPE_c,2) +
                              DEPTH_c:POP_STATUS +
                              DEPTH_c:poly(SITE_SLOPE_100m_c,2) +
                              DEPTH_c:poly(SITE_SLOPE_400m_c,2) +
                              DEPTH_c:poly(ISL_SLOPE_c,2) +
                              OBS_YEAR +
                              (1|OBS_YEAR:REGION) + 
                              (1|OBS_YEAR:ISLAND) + 
                              (1|REGION)  + 
                              (1|SITE) +
                              (1+DEPTH_c||ISLAND),
                            data=PISC1, 
                            family = Gamma(link=log),
                            prior = priors_PISC,
                            chains = 4,
                            cores  = 4,
                            iter = 1000,
                            warmup = 100,
                            sample_prior = "only")

## NEW MAKE INIT LIST 
start_names <- unique(paste0('sd_',PISC.Gamma.brms.prior$prior$group,'__',PISC.Gamma.brms.prior$prior$coef)[PISC.Gamma.brms.prior$prior$class=='sd'])
start_names <- start_names[!grepl('.__$', start_names)]

SLIST <- as.list(start_names)
names(SLIST) <- start_names
SLIST <- lapply(SLIST, function(x) 0.1)

inilist <- c(SLIST,
             list(b_Intercept=4, 
                  b_DEPTH_c=0,
                  `b_DEPTH_c:POP_STATUSP`=0,
                  `b_DEPTH_c:ISL_SLOPE_c`=0,
                  `b_DEPTH_c:SITE_SLOPE_100m_c`=0,
                  `b_DEPTH_c:SITE_SLOPE_400m_c`=0,
                  b_OBS_YEAR_2011=0,
                  b_OBS_YEAR_2012=0,
                  b_OBS_YEAR_2013=0,
                  b_OBS_YEAR_2014=0,
                  b_POP_STATUSP=0,
                  b_ISL_SLOPE_c=0,
                  b_SITE_SLOPE_100m_c=0,
                  b_SITE_SLOPE_400m_c=0,
                  shape=2.5))

PISC.Gamma.brms.full = brm(PISCIVORE ~ DEPTH_c +
                             POP_STATUS +
                             poly(SITE_SLOPE_100m_c,2) +
                             poly(SITE_SLOPE_400m_c,2) +
                             poly(ISL_SLOPE_c,2) +
                             DEPTH_c:POP_STATUS +
                             DEPTH_c:poly(SITE_SLOPE_100m_c,2) +
                             DEPTH_c:poly(SITE_SLOPE_400m_c,2) +
                             DEPTH_c:poly(ISL_SLOPE_c,2) +
                             OBS_YEAR +
                             (1|OBS_YEAR:REGION) + 
                             (1|OBS_YEAR:ISLAND) + 
                             (1|REGION)  + 
                             (1|SITE) +
                             (1+DEPTH_c||ISLAND),
                           inits = function() inilist,
                           data=PISC1, 
                           family = Gamma(link=log),
                           prior = priors_PISC,
                           chains = 4,
                           cores  = 4,
                           thin = 4,
                           iter = 2500,
                           warmup = 500, 
                           control = list(adapt_delta=0.95))

#2) Taking a look at the model:
pairs(PISC.Gamma.brms.full)
mcmc_plot(PISC.Gamma.brms.full, type='trace')
pp = brms::pp_check(PISC.Gamma.brms.full, type = 'ecdf_overlay', nsamples=200) 
pp + theme_bw()+ xlim(c(0,500)) + xlab('Density (g.m2)') + ylab("Cumulative probability")
brms::pp_check(PISC.Gamma.brms.full, type = 'loo_pit_qq', nsamples=200) 
brms::pp_check(PISC.Gamma.brms.full, type = 'error_scatter_avg_vs_x', nsamples=200, x='DEPTH_c') 

summary(PISC.Gamma.brms.full)
conditional_effects(PISC.Gamma.brms.full)

#CHECK TO SEE IF YOU NEED RANDOM DEPTH; THEN SEE IF NEED RANDOM OBS_YEAR terms:
#1. Need random slope for depth?
PISC.Gamma.brms.full.1 = brm(PISCIVORE ~ DEPTH_c +
                                  POP_STATUS +
                                  poly(SITE_SLOPE_100m_c,2) +
                                  poly(SITE_SLOPE_400m_c,2) +
                                  poly(ISL_SLOPE_c,2) +
                                  DEPTH_c:POP_STATUS +
                                  DEPTH_c:poly(SITE_SLOPE_100m_c,2) +
                                  DEPTH_c:poly(SITE_SLOPE_400m_c,2) +
                                  DEPTH_c:poly(ISL_SLOPE_c,2) +
                                  OBS_YEAR +
                                  (1|OBS_YEAR:REGION) + 
                                  (1|OBS_YEAR:ISLAND) + 
                                  (1|REGION)  + 
                                  (1|SITE) +
                                  (1+DEPTH_c||ISLAND),
                                inits = function() inilist,
                                data=PISC1, 
                                family = Gamma(link=log),
                                prior = priors_PISC,
                                chains = 4,
                                cores  = 4,
                                thin = 4,
                                iter = 1000,
                                warmup = 100, 
                                control = list(adapt_delta=0.95))

PISC.Gamma.brms.full.1a = brm(PISCIVORE ~ DEPTH_c +
                                   POP_STATUS +
                                   poly(SITE_SLOPE_100m_c,2) +
                                   poly(SITE_SLOPE_400m_c,2) +
                                   poly(ISL_SLOPE_c,2) +
                                   DEPTH_c:POP_STATUS +
                                   DEPTH_c:poly(SITE_SLOPE_100m_c,2) +
                                   DEPTH_c:poly(SITE_SLOPE_400m_c,2) +
                                   DEPTH_c:poly(ISL_SLOPE_c,2) +
                                   OBS_YEAR +
                                   (1|OBS_YEAR:REGION) + 
                                   (1|OBS_YEAR:ISLAND) + 
                                   (1|REGION)  + 
                                   (1|SITE) +
                                   (1|ISLAND),
                                 inits = function() inilist,
                                 data=PISC1, 
                                 family = Gamma(link=log),
                                 prior = c(set_prior('normal(0,2)', class='b'),
                                           set_prior('cauchy(0,5)', class='sd'),
                                           set_prior('normal(2.95,0.35)', class='Intercept')),
                                 chains = 1,
                                 cores  = 1,
                                 thin = 4,
                                 iter = 1000,
                                 warmup = 100)#, 
#control = list(adapt_delta=0.95))

loofull<-loo(PISC.Gamma.brms.full.1)
looA<-loo(PISC.Gamma.brms.full.1a)

loo_compare(loofull,looA)

mcmc_plot(PISC.Gamma.brms.full.1a, type='trace')
pp = brms::pp_check(PISC.Gamma.brms.full.1a, type = 'ecdf_overlay', nsamples=200) 
pp + theme_bw()+ xlim(c(0,500)) + xlab('Density (g.m2)') + ylab("Cumulative probability")
brms::pp_check(PISC.Gamma.brms.full.1a, type = 'loo_pit_qq', nsamples=200) 

#2. Need random OBS_YEAR terms?
PISC.Gamma.brms.full.1b = brm(PISCIVORE ~ DEPTH_c +
                                   POP_STATUS +
                                   poly(SITE_SLOPE_100m_c,2) +
                                   poly(SITE_SLOPE_400m_c,2) +
                                   poly(ISL_SLOPE_c,2) +
                                   DEPTH_c:POP_STATUS +
                                   DEPTH_c:poly(SITE_SLOPE_100m_c,2) +
                                   DEPTH_c:poly(SITE_SLOPE_400m_c,2) +
                                   DEPTH_c:poly(ISL_SLOPE_c,2) +
                                   OBS_YEAR +
                                   (1|OBS_YEAR:ISLAND) + 
                                   (1|REGION)  + 
                                   (1|SITE) +
                                   (1+DEPTH_c||ISLAND),
                                 inits = function() inilist,
                                 data=PISC1, 
                                 family = Gamma(link=log),
                                 prior = priors_PISC,
                                 chains = 1,
                                 cores  = 1,
                                 thin = 4,
                                 iter = 1000,
                                 warmup = 100)#, 
#control = list(adapt_delta=0.95))

PISC.Gamma.brms.full.1c = brm(PISCIVORE ~ DEPTH_c +
                                   POP_STATUS +
                                   poly(SITE_SLOPE_100m_c,2) +
                                   poly(SITE_SLOPE_400m_c,2) +
                                   poly(ISL_SLOPE_c,2) +
                                   DEPTH_c:POP_STATUS +
                                   DEPTH_c:poly(SITE_SLOPE_100m_c,2) +
                                   DEPTH_c:poly(SITE_SLOPE_400m_c,2) +
                                   DEPTH_c:poly(ISL_SLOPE_c,2) +
                                   OBS_YEAR +
                                   (1|OBS_YEAR:REGION) + 
                                   (1|REGION)  + 
                                   (1|SITE) +
                                   (1+DEPTH_c||ISLAND),
                                 inits = function() inilist,
                                 data=PISC1, 
                                 family = Gamma(link=log),
                                 prior = priors_PISC,
                                 chains = 1,
                                 cores  = 1,
                                 thin = 4,
                                 iter = 1000,
                                 warmup = 100)#, 
#control = list(adapt_delta=0.95))

PISC.Gamma.brms.full.1d = brm(PISCIVORE ~ DEPTH_c +
                                   POP_STATUS +
                                   poly(SITE_SLOPE_100m_c,2) +
                                   poly(SITE_SLOPE_400m_c,2) +
                                   poly(ISL_SLOPE_c,2) +
                                   DEPTH_c:POP_STATUS +
                                   DEPTH_c:poly(SITE_SLOPE_100m_c,2) +
                                   DEPTH_c:poly(SITE_SLOPE_400m_c,2) +
                                   DEPTH_c:poly(ISL_SLOPE_c,2) +
                                   OBS_YEAR +
                                   (1|REGION)  + 
                                   (1|SITE) +
                                   (1+DEPTH_c||ISLAND),
                                 inits = function() inilist,
                                 data=PISC1, 
                                 family = Gamma(link=log),
                                 prior = priors_PISC,
                                 chains = 1,
                                 cores  = 1,
                                 thin = 4,
                                 iter = 1000,
                                 warmup = 100)#, 
#control = list(adapt_delta=0.95))

PISC.Gamma.brms.full.1e = brm(PISC ~ DEPTH_c +
                                   POP_STATUS +
                                   poly(SITE_SLOPE_100m_c,2) +
                                   poly(SITE_SLOPE_400m_c,2) +
                                   poly(ISL_SLOPE_c,2) +
                                   DEPTH_c:POP_STATUS +
                                   DEPTH_c:poly(SITE_SLOPE_100m_c,2) +
                                   DEPTH_c:poly(SITE_SLOPE_400m_c,2) +
                                   DEPTH_c:poly(ISL_SLOPE_c,2) +
                                   OBS_YEAR +
                                   (1|REGION)  + 
                                   (1|SITE) +
                                   (1|ISLAND),
                                 inits = function() inilist,
                                 data=PISC1, 
                                 family = Gamma(link=log),
                                 prior = c(set_prior('normal(0,2)', class='b'),
                                           set_prior('cauchy(0,5)', class='sd'),
                                           set_prior('normal(2.95,0.35)', class='Intercept')),
                                 chains = 1,
                                 cores  = 1,
                                 thin = 4,
                                 iter = 1000,
                                 warmup = 100)#, 
#control = list(adapt_delta=0.95))

looB<-loo(PISC.Gamma.brms.full.1b)
looC<-loo(PISC.Gamma.brms.full.1c)
looD<-loo(PISC.Gamma.brms.full.1d)
looE<-loo(PISC.Gamma.brms.full.1e)

loo_compare(looB,looC,looD,looE)

PISC.Gamma.brms.full.1d = brm(PISCIVORE ~ DEPTH_c +
                                   POP_STATUS +
                                   poly(SITE_SLOPE_100m_c,2) +
                                   poly(SITE_SLOPE_400m_c,2) +
                                   poly(ISL_SLOPE_c,2) +
                                   DEPTH_c:POP_STATUS +
                                   DEPTH_c:poly(SITE_SLOPE_100m_c,2) +
                                   DEPTH_c:poly(SITE_SLOPE_400m_c,2) +
                                   DEPTH_c:poly(ISL_SLOPE_c,2) +
                                   OBS_YEAR +
                                   (1|REGION)  + 
                                   (1|SITE) +
                                   (1+DEPTH_c||ISLAND),
                                 inits = function() inilist,
                                 data=PISC1, 
                                 family = Gamma(link=log),
                                 prior = priors_PISC,
                                 chains = 4,
                                 cores  = 4,
                                 thin = 4,
                                 iter = 2500,
                                 warmup = 500,
                                 control = list(adapt_delta=0.95))

#Taking a look at the model:
mcmc_plot(PISC.Gamma.brms.full.1d, type='trace')
pp = brms::pp_check(PISC.Gamma.brms.full.1d, type = 'ecdf_overlay', nsamples=200) 
pp + theme_bw()+ xlim(c(0,500)) + xlab('Density (g.m2)') + ylab("Cumulative probability")
brms::pp_check(PISC.Gamma.brms.full.1d, type = 'loo_pit_qq', nsamples=200) 
brms::pp_check(PISC.Gamma.brms.full.1d, type = 'error_scatter_avg_vs_x', nsamples=200, x='DEPTH_c') 

summary(PISC.Gamma.brms.full.1d)
conditional_effects(PISC.Gamma.brms.full.1d)




#Plots of variance components:
mcmc_plot(PISC.Gamma.brms.full, pars=c('sd_SITE__Intercept','sd_ISLAND__Intercept','sd_REGION__Intercept'), fixed = T) + 
  #xlim(0,2.5) +
  xlab("PISC model Standard deviation") +
  theme_bw()

mcmc_plot(PISC.Gamma.brms.full, pars='sd') + 
  #xlim(0,2.5) +
  xlab("PISC model Standard deviation") +
  theme_bw()

pairs(PISC.Gamma.brms.full, pars=c('sd_SITE__Intercept','sd_ISLAND__Intercept','sd_REGION__Intercept'), fixed = T)

sds <- posterior_samples(PISC.Gamma.brms.full, pars = 'sd')

hist(sds$sd_REGION__Intercept/sds$sd_ISLAND__Intercept)
mean(sds$sd_REGION__Intercept > sds$sd_ISLAND__Intercept)

- sds$sd_REGION__Intercept - sds$sd_ISLAND__Intercept

rel_imp <- apply(sds,1,function(x) x/sum(x))

bayes_R2(PISC.Gamma.brms.full)
performance::r2_bayes(PISC.Gamma.brms.full) #marginal = fixed effs


##################################################################################################################

#PLANKTIVORE

#Remove NAs in slope predictors for polys to work:
PLANK1<-PLANK[!is.na(PLANK$SITE_SLOPE_100m_c),]
PLANK1<-PLANK1[!is.na(PLANK1$SITE_SLOPE_400m_c),]
PLANK1<-PLANK1[!is.na(PLANK1$ISL_SLOPE_c),]
summary(PLANK1)

#1) Build the model:

#Trophic group proportions from MacNeil et al. (2015)
#PLANK:
#Prop: 0.3 mean, 0.4 sd?
hist(rbeta(10000,15,35))
hist(rbeta(10000,3,7))
mean(log(rlnorm(10000, 4.6, 0.21)*rbeta(10000,15,35)))
sd(log(rlnorm(10000, 4.6, 0.21)*rbeta(10000,15,35)))
#3.36,0.30

#Set priors:
priors_PLANK = c(set_prior('normal(0,2)', class='b'),
                 set_prior('cauchy(0,5)', class='sd'),
                 set_prior('cauchy(0,1)', class='sd', coef = 'DEPTH_c', group = 'ISLAND'),
                 set_prior('normal(3.36,1)', class='Intercept'))

PLANK1$POP_STATUS <- factor(PLANK1$POP_STATUS, levels = c('U', 'P'))

PLANK.Gamma.brms.prior = brm(PLANKTIVORE ~ DEPTH_c +
                               POP_STATUS +
                               poly(SITE_SLOPE_100m_c,2) +
                               poly(SITE_SLOPE_400m_c,2) +
                               poly(ISL_SLOPE_c,2) +
                               DEPTH_c:POP_STATUS +
                               DEPTH_c:poly(SITE_SLOPE_100m_c,2) +
                               DEPTH_c:poly(SITE_SLOPE_400m_c,2) +
                               DEPTH_c:poly(ISL_SLOPE_c,2) +
                               OBS_YEAR +
                               (1|OBS_YEAR:REGION) + 
                               (1|OBS_YEAR:ISLAND) + 
                               (1|REGION)  + 
                               (1|SITE) +
                               (1+DEPTH_c||ISLAND),
                             data=PLANK1, 
                             family = Gamma(link=log),
                             prior = priors_PLANK,
                             chains = 4,
                             cores  = 4,
                             iter = 1000,
                             warmup = 100,
                             sample_prior = "only")

## NEW MAKE INIT LIST 
start_names <- unique(paste0('sd_',PLANK.Gamma.brms.prior$prior$group,'__',PLANK.Gamma.brms.prior$prior$coef)[PLANK.Gamma.brms.prior$prior$class=='sd'])
start_names <- start_names[!grepl('.__$', start_names)]

SLIST <- as.list(start_names)
names(SLIST) <- start_names
SLIST <- lapply(SLIST, function(x) 0.1)

inilist <- c(SLIST,
             list(b_Intercept=4, 
                  b_DEPTH_c=0,
                  `b_DEPTH_c:POP_STATUSP`=0,
                  `b_DEPTH_c:ISL_SLOPE_c`=0,
                  `b_DEPTH_c:SITE_SLOPE_100m_c`=0,
                  `b_DEPTH_c:SITE_SLOPE_400m_c`=0,
                  b_OBS_YEAR_2011=0,
                  b_OBS_YEAR_2012=0,
                  b_OBS_YEAR_2013=0,
                  b_OBS_YEAR_2014=0,
                  b_POP_STATUSP=0,
                  b_ISL_SLOPE_c=0,
                  b_SITE_SLOPE_100m_c=0,
                  b_SITE_SLOPE_400m_c=0,
                  shape=2.5))

PLANK.Gamma.brms.full = brm(PLANKTIVORE ~ DEPTH_c +
                              POP_STATUS +
                              poly(SITE_SLOPE_100m_c,2) +
                              poly(SITE_SLOPE_400m_c,2) +
                              poly(ISL_SLOPE_c,2) +
                              DEPTH_c:POP_STATUS +
                              DEPTH_c:poly(SITE_SLOPE_100m_c,2) +
                              DEPTH_c:poly(SITE_SLOPE_400m_c,2) +
                              DEPTH_c:poly(ISL_SLOPE_c,2) +
                              OBS_YEAR +
                              (1|OBS_YEAR:REGION) + 
                              (1|OBS_YEAR:ISLAND) + 
                              (1|REGION)  + 
                              (1|SITE) +
                              (1+DEPTH_c||ISLAND),
                            inits = function() inilist,
                            data=PLANK1, 
                            family = Gamma(link=log),
                            prior = priors_PLANK,
                            chains = 4,
                            cores  = 4,
                            thin = 4,
                            iter = 2500,
                            warmup = 500, 
                            control = list(adapt_delta=0.95))

#2) Taking a look at the model:
mcmc_plot(PLANK.Gamma.brms.full, type='trace')
pp = brms::pp_check(PLANK.Gamma.brms.full, type = 'ecdf_overlay', nsamples=200) 
pp + theme_bw()+ xlim(c(0,500)) + xlab('Density (g.m2)') + ylab("Cumulative probability")
brms::pp_check(PLANK.Gamma.brms.full, type = 'loo_pit_qq', nsamples=200) 
brms::pp_check(PLANK.Gamma.brms.full, type = 'error_scatter_avg_vs_x', nsamples=200, x='DEPTH_c') 

summary(PLANK.Gamma.brms.full)
conditional_effects(PLANK.Gamma.brms.full)
mcmc_plot(PLANK.Gamma.brms.full, pars=c('sd_SITE__Intercept','sd_ISLAND__Intercept','sd_REGION__Intercept'), fixed = T) + 
  #xlim(0,1.75) +
  xlab("PLANK model Standard deviation") +
  theme_bw()

mcmc_plot(PLANK.Gamma.brms.full, pars='sd') + 
  #xlim(0,2.5) +
  xlab("PLANK model Standard deviation") +
  theme_bw()

pairs(PLANK.Gamma.brms.full, pars=c('sd_SITE__Intercept','sd_ISLAND__Intercept','sd_REGION__Intercept'), fixed = T)

sds <- posterior_samples(PLANK.Gamma.brms.full, pars = 'sd')

hist(sds$sd_REGION__Intercept/sds$sd_ISLAND__Intercept)
mean(sds$sd_REGION__Intercept > sds$sd_ISLAND__Intercept)

- sds$sd_REGION__Intercept - sds$sd_ISLAND__Intercept

rel_imp <- apply(sds,1,function(x) x/sum(x))

bayes_R2(PLANK.Gamma.brms.full)
performance::r2_bayes(PLANK.Gamma.brms.full) #marginal = fixed effs


##################################################################################################################

#PRIMARY CONSUMER

#Remove NAs in slope predictors for polys to work:
PRIM1<-PRIM[!is.na(PRIM$SITE_SLOPE_100m_c),]
PRIM1<-PRIM1[!is.na(PRIM1$SITE_SLOPE_400m_c),]
PRIM1<-PRIM1[!is.na(PRIM1$ISL_SLOPE_c),]
summary(PRIM1)

#1) Build the model:

#Trophic group proportions from MacNeil et al. (2015)
#PRIM: 
#(a+b+c+d in MacNeil et al. Extended Data Figure 6: Average reef fish functional group across a biomass gradient)
#Prop: 0.8 mean?
hist(rbeta(10000,40,10))
hist(rbeta(10000,8,2))
mean(log(rlnorm(10000, 4.6, 0.21)*rbeta(10000,40,10)))
sd(log(rlnorm(10000, 4.6, 0.21)*rbeta(10000,40,10)))
#4.37, 0.22

#Set priors:
priors_PRIM = c(set_prior('normal(0,2)', class='b'),
                set_prior('cauchy(0,5)', class='sd'),
                set_prior('cauchy(0,1)', class='sd', coef = 'DEPTH_c', group = 'ISLAND'),
                set_prior('normal(4.37, 1)', class='Intercept'))

PRIM1$POP_STATUS <- factor(PRIM1$POP_STATUS, levels = c('U', 'P'))

PRIM.Gamma.brms.prior = brm(PRIMARY ~ DEPTH_c +
                              POP_STATUS +
                              poly(SITE_SLOPE_100m_c,2) +
                              poly(SITE_SLOPE_400m_c,2) +
                              poly(ISL_SLOPE_c,2) +
                              DEPTH_c:POP_STATUS +
                              DEPTH_c:poly(SITE_SLOPE_100m_c,2) +
                              DEPTH_c:poly(SITE_SLOPE_400m_c,2) +
                              DEPTH_c:poly(ISL_SLOPE_c,2) +
                              OBS_YEAR +
                              (1|OBS_YEAR:REGION) + 
                              (1|OBS_YEAR:ISLAND) + 
                              (1|REGION)  + 
                              (1|SITE) +
                              (1+DEPTH_c||ISLAND),
                            data=PRIM1, 
                            family = Gamma(link=log),
                            prior = priors_PRIM,
                            chains = 4,
                            cores  = 4,
                            iter = 1000,
                            warmup = 100,
                            sample_prior = "only")

## NEW MAKE INIT LIST 
start_names <- unique(paste0('sd_',PRIM.Gamma.brms.prior$prior$group,'__',PRIM.Gamma.brms.prior$prior$coef)[PRIM.Gamma.brms.prior$prior$class=='sd'])
start_names <- start_names[!grepl('.__$', start_names)]

SLIST <- as.list(start_names)
names(SLIST) <- start_names
SLIST <- lapply(SLIST, function(x) 0.1)

inilist <- c(SLIST,
             list(b_Intercept=4, 
                  b_DEPTH_c=0,
                  `b_DEPTH_c:POP_STATUSP`=0,
                  `b_DEPTH_c:ISL_SLOPE_c`=0,
                  `b_DEPTH_c:SITE_SLOPE_100m_c`=0,
                  `b_DEPTH_c:SITE_SLOPE_400m_c`=0,
                  b_OBS_YEAR_2011=0,
                  b_OBS_YEAR_2012=0,
                  b_OBS_YEAR_2013=0,
                  b_OBS_YEAR_2014=0,
                  b_POP_STATUSP=0,
                  b_ISL_SLOPE_c=0,
                  b_SITE_SLOPE_100m_c=0,
                  b_SITE_SLOPE_400m_c=0,
                  shape=2.5))

PRIM.Gamma.brms.full = brm(PRIMARY ~ DEPTH_c +
                             POP_STATUS +
                             poly(SITE_SLOPE_100m_c,2) +
                             poly(SITE_SLOPE_400m_c,2) +
                             poly(ISL_SLOPE_c,2) +
                             DEPTH_c:POP_STATUS +
                             DEPTH_c:poly(SITE_SLOPE_100m_c,2) +
                             DEPTH_c:poly(SITE_SLOPE_400m_c,2) +
                             DEPTH_c:poly(ISL_SLOPE_c,2) +
                             OBS_YEAR +
                             (1|OBS_YEAR:REGION) + 
                             (1|OBS_YEAR:ISLAND) + 
                             (1|REGION)  + 
                             (1|SITE) +
                             (1+DEPTH_c||ISLAND),
                           inits = function() inilist,
                           data=PRIM1, 
                           family = Gamma(link=log),
                           prior = priors_PRIM,
                           chains = 4,
                           cores  = 4,
                           thin = 4,
                           iter = 2500,
                           warmup = 500, 
                           control = list(adapt_delta=0.95))

#2) Taking a look at the model:
mcmc_plot(PRIM.Gamma.brms.full, type='trace')
pp = brms::pp_check(PRIM.Gamma.brms.full, type = 'ecdf_overlay', nsamples=200) 
pp + theme_bw()+ xlim(c(0,500)) + xlab('Density (g.m2)') + ylab("Cumulative probability")
brms::pp_check(PRIM.Gamma.brms.full, type = 'loo_pit_qq', nsamples=200) 
brms::pp_check(PRIM.Gamma.brms.full, type = 'error_scatter_avg_vs_x', nsamples=200, x='DEPTH_c') 

summary(PRIM.Gamma.brms.full)
conditional_effects(PRIM.Gamma.brms.full)
mcmc_plot(PRIM.Gamma.brms.full, pars=c('sd_SITE__Intercept','sd_ISLAND__Intercept','sd_REGION__Intercept'), fixed = T) + 
  #xlim(0,1.75) +
  xlab("PRIM model Standard deviation") +
  theme_bw()

pairs(PRIM.Gamma.brms.full, pars=c('sd_SITE__Intercept','sd_ISLAND__Intercept','sd_REGION__Intercept'), fixed = T)

sds <- posterior_samples(PRIM.Gamma.brms.full, pars = 'sd')

hist(sds$sd_REGION__Intercept/sds$sd_ISLAND__Intercept)
mean(sds$sd_REGION__Intercept > sds$sd_ISLAND__Intercept)

- sds$sd_REGION__Intercept - sds$sd_ISLAND__Intercept

rel_imp <- apply(sds,1,function(x) x/sum(x))

bayes_R2(PRIM.Gamma.brms.full)
performance::r2_bayes(PRIM.Gamma.brms.full) #marginal = fixed effs


###################################################################

####SECONDARY CONSUMER

#Remove NAs in slope predictors for polys to work:
fish1<-fish[!is.na(fish$SITE_SLOPE_100m_c),]
fish1<-fish1[!is.na(fish1$SITE_SLOPE_400m_c),]
fish1<-fish1[!is.na(fish1$ISL_SLOPE_c),]
summary(fish1)

#SEC:
#(f+g+h in MacNeil et al. Extended Data Figure 6: Average reef fish functional group across a biomass gradient)
#Prop: 0.65 mean?
hist(rbeta(10000,32.5,17.5))
hist(rbeta(10000,6.5,3.5))
mean(log(rlnorm(10000, 4.6, 0.21)*rbeta(10000,32.5,17.5)))
sd(log(rlnorm(10000, 4.6, 0.21)*rbeta(10000,32.5,17.5)))
#4.16,0.23


#1) Build the model:

#Set priors:
priors_SEC = c(set_prior('normal(0,2)', class='b'),
               set_prior('cauchy(0,5)', class='sd'),
               set_prior('cauchy(0,1)', class='sd', coef = 'DEPTH_c', group = 'ISLAND'),
               set_prior('normal(4.16,0.5)', class='Intercept')) 

#Look to see if priors are way out by drawing on the prior distribution with 'sample_prior = "only"':

fish1$POP_STATUS <- factor(fish1$POP_STATUS, levels = c('U', 'P'))
SEC.Gamma.brms.prior = brm(SECONDARY ~ DEPTH_c +
                             POP_STATUS +
                             poly(SITE_SLOPE_100m_c,2) +
                             poly(SITE_SLOPE_400m_c,2) +
                             poly(ISL_SLOPE_c,2) +
                             DEPTH_c:POP_STATUS +
                             DEPTH_c:poly(SITE_SLOPE_100m_c,2) +
                             DEPTH_c:poly(SITE_SLOPE_400m_c,2) +
                             DEPTH_c:poly(ISL_SLOPE_c,2) +
                             OBS_YEAR +
                             (1|OBS_YEAR:REGION) + 
                             (1|OBS_YEAR:ISLAND) + 
                             (1|REGION)  + 
                             (1|SITE) +
                             (1+DEPTH_c||ISLAND),
                           data=fish1, 
                           family = Gamma(link=log),
                           prior = priors_SEC,
                           chains = 4,
                           cores  = 4,
                           iter = 1000,
                           warmup = 100,
                           sample_prior = "only")

prior_summary(SEC.Gamma.brms.prior)
hist(log(fitted(SEC.Gamma.brms.prior)[, 1]))
summary(SEC.Gamma.brms.prior)

#Fit model with chosen priors:

## NEW MAKE INIT LIST 
start_names <- unique(paste0('sd_',SEC.Gamma.brms.prior$prior$group,'__',SEC.Gamma.brms.prior$prior$coef)[SEC.Gamma.brms.prior$prior$class=='sd'])
start_names <- start_names[!grepl('.__$', start_names)]

SLIST <- as.list(start_names)
names(SLIST) <- start_names
SLIST <- lapply(SLIST, function(x) 0.1)

inilist <- c(SLIST,
             list(b_Intercept=4, 
                  b_DEPTH_c=0,
                  `b_DEPTH_c:POP_STATUSP`=0,
                  `b_DEPTH_c:ISL_SLOPE_c`=0,
                  `b_DEPTH_c:SITE_SLOPE_100m_c`=0,
                  `b_DEPTH_c:SITE_SLOPE_400m_c`=0,
                  b_OBS_YEAR_2011=0,
                  b_OBS_YEAR_2012=0,
                  b_OBS_YEAR_2013=0,
                  b_OBS_YEAR_2014=0,
                  b_POP_STATUSP=0,
                  b_ISL_SLOPE_c=0,
                  b_SITE_SLOPE_100m_c=0,
                  b_SITE_SLOPE_400m_c=0,
                  shape=2.5))

SEC.Gamma.brms.full = brm(SECONDARY ~ DEPTH_c +
                            POP_STATUS +
                            poly(SITE_SLOPE_100m_c,2) +
                            poly(SITE_SLOPE_400m_c,2) +
                            poly(ISL_SLOPE_c,2) +
                            DEPTH_c:POP_STATUS +
                            DEPTH_c:poly(SITE_SLOPE_100m_c,2) +
                            DEPTH_c:poly(SITE_SLOPE_400m_c,2) +
                            DEPTH_c:poly(ISL_SLOPE_c,2) +
                            OBS_YEAR +
                            (1|OBS_YEAR:REGION) + 
                            (1|OBS_YEAR:ISLAND) + 
                            (1|REGION)  + 
                            (1|SITE) +
                            (1+DEPTH_c||ISLAND),
                          inits = function() inilist,
                          data=fish1, 
                          family = Gamma(link=log),
                          prior = priors_SEC,
                          chains = 4,
                          cores  = 4,
                          thin = 4,
                          iter = 2500,
                          warmup = 500, 
                          control = list(adapt_delta=0.95))


#2) Taking a look at the model:
mcmc_plot(SEC.Gamma.brms.full, type='trace')
pp = brms::pp_check(SEC.Gamma.brms.full, type = 'ecdf_overlay', nsamples=200) 
pp + theme_bw()+ xlim(c(0,500)) + xlab('Density (g.m2)') + ylab("Cumulative probability")
brms::pp_check(SEC.Gamma.brms.full, type = 'loo_pit_qq', nsamples=200) 
brms::pp_check(SEC.Gamma.brms.full, type = 'error_scatter_avg_vs_x', nsamples=200, x='DEPTH_c') 

summary(SEC.Gamma.brms.full)
conditional_effects(SEC.Gamma.brms.full)
mcmc_plot(SEC.Gamma.brms.full, pars=c('sd_SITE__Intercept','sd_ISLAND__Intercept','sd_REGION__Intercept'), fixed = T) + 
  #xlim(0,2.5) +
  xlab("SEC model Standard deviation") +
  theme_bw()

pairs(SEC.Gamma.brms.full, pars=c('sd_SITE__Intercept','sd_ISLAND__Intercept','sd_REGION__Intercept'), fixed = T)

sds <- posterior_samples(SEC.Gamma.brms.full, pars = 'sd')

hist(sds$sd_REGION__Intercept/sds$sd_ISLAND__Intercept)
mean(sds$sd_REGION__Intercept > sds$sd_ISLAND__Intercept)

- sds$sd_REGION__Intercept - sds$sd_ISLAND__Intercept

rel_imp <- apply(sds,1,function(x) x/sum(x))

bayes_R2(SEC.Gamma.brms.full)
performance::r2_bayes(SEC.Gamma.brms.full) #marginal = fixed effs


##################################################################################################################
#end