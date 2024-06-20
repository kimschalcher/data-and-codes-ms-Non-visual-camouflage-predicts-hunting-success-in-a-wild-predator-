# load packages ####
#install.packages('see')
library(report)
library(ggpubr)
library(lubridate)
library(tidyverse)
library(data.table)
library(lme4)
library(lmerTest)
library(MuMIn)
library(DHARMa)
library(mgcv)
library(viridis)
library(hrbrthemes)
library(suncalc)
library(NISTunits)
library(sp)
library(sf)
library(raster)
library(spatial)
library(rgdal)
library(rgeos)
library(fasterize)
library(sjPlot)
library(DHARMa)
library(geosphere)
library(gratia)
library(performance)
library(zoo)
library(slider)
# Load data #### 

## land data ####

all_indiv_land_complete_final_data<-fread("/Volumes/aroulin1/chouette/D2c/papers/data/all_indiv_land_complete_final_data.csv") 


dim(all_indiv_land_complete_final_data) 

#all_indiv_land_complete_final_data<- distinct(all_indiv_land_complete_final_data)
unique(all_indiv_land_complete_final_data$TagID)
indiv_per_year<-all_indiv_land_complete_final_data %>% 
  group_by(TagID) %>% 
  slice(1)
table(indiv_per_year$year)



#model Landin force

all_indiv_land_complete_final_data2<-all_indiv_land_complete_final_data %>% 
  drop_na(WindSpeed)

model_strike_force_divecat<- lmer(log(impact_force) ~  
                                    ###landing_cat*
                                    
                                    ###sex+
                                    #WindSpeed+
                                    #landing_cat*sex+
                                    
                                    landing_cat/(1+sex)-1+
                                    #landing_cat:WindSpeed+
                                    ##landing_cat:WindSpeed:sex+
                                    (1|TagID/night_nb_real),
                                  data = all_indiv_land_complete_final_data, REML = F)

summary(model_strike_force_divecat)
emmeanns<-emmeans(model_strike_force_divecat, pairwise~ landing_cat*sex, type="resp")

##
dredge_model_strike_force_divecat<-dredge(model_strike_force_divecat, rank="AIC", extra= "R^2")
dredge_model_strike_force_divecat
summary(model.avg(dredge_model_strike_force_divecat, subset = delta <= 2))

plot_model(model_strike_force_divecat,type="emm", terms = c("landing_cat","sex"))


#### summary observation model data ####

night_nb_landing_cat<-all_indiv_land_complete_final_data %>% 
  as_tibble() %>% 
  dplyr::group_by(TagID, landing_cat) %>% 
  dplyr::summarize(n=n(),
                   obs_per_night= n()/max(night_nb_real, na.rm = T),
                   sex=first(sex)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(sex,landing_cat) %>% 
  dplyr::summarize(night_nb_mean=mean(obs_per_night),
                   nighht_nb_sd=sd(obs_per_night))



## summary table results ####
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(webshot)

tab_model(model_strike_force_divecat,string.est = "Estimate", 
          auto.label = T,
          string.ci = "Conf. Int.", string.p = "p-value",transform = "exp" )

table(all_indiv_land_complete_final_data$landing_cat)
###export tab mode

tab_model(model_strike_force_divecat,string.est = "Estimate", 
          auto.label = T,
          string.ci = "Conf. Int.", string.p = "p-value",transform = "exp" , file = "/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/table_result1.html" )


webshot("/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/table_result1.html", "/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/table_result1.pdf")


## plots ####

simulationOutput <- simulateResiduals(fittedModel = model_strike_force_divecat, plot = F,refit = F)
plot(simulationOutput)
summary( model_strike_force_divecat)



# plot final ####
#model refit

model_strike_force_divecat_refit<- lmer(log(impact_force) ~  
                                          landing_cat*
                                          #WindSpeed+
                                          sex+
                                          #landing_cat:sex+
                                          #landing_cat/(1+sex)-1+
                                          #landing_cat:WindSpeed+
                                          ##landing_cat:WindSpeed:sex+
                                          (1|TagID),
                                        data = all_indiv_land_complete_final_data, REML = F)


summary(model_strike_force_divecat_refit)
#predict_data
library(arm)
nsim <- 2000
bsim <- sim(model_strike_force_divecat_refit, n.sim=nsim)
str(bsim)
colnames (bsim @ fixef) <- names(fixef(model_strike_force_divecat_refit)) # for arm 1.6e10
fixef(model_strike_force_divecat_refit)
apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))

## get model prediction

newdata_landingcat<- expand.grid (landing_cat = c("dive","landing"),
                                  sex=c("Female", "Male"))

Xmat_landingcat <- model.matrix(~ landing_cat*sex, data=newdata_landingcat)   # the model matrix
newdata_landingcat$pred <- Xmat_landingcat%*%fixef(model_strike_force_divecat_refit)  

predmat_landingcat <- matrix(nrow=nrow(newdata_landingcat),ncol=nsim)
for(i in 1:nsim) predmat_landingcat[,i] <-Xmat_landingcat  %*% bsim@fixef[i,]
newdata_landingcat$lwr <- apply(predmat_landingcat,1,quantile,probs=0.025)
newdata_landingcat$upr <- apply(predmat_landingcat,1,quantile,probs=0.975)

#predmat_impact_cat$variable scale.bt = (newdatdistnb$distnb.sc* sd(tablesansshrew$distnb)) + mean(tablesansshrew$distnb) # back transformer variable réponse

# summary(tablefinale$distnb)
# summary(newdatdistnb$distnb.sc.bt)
###

# plot ####
newdata_landingcat2<-newdata_landingcat %>% 
  arrange(landing_cat)

summary(log(all_indiv_land_complete_final_data$impact_force))



pdf(file="/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/Result_1.pdf")


plot(as.numeric(newdata_landingcat2$landing_cat),newdata_landingcat2$pred,type="n",las=1,xlab="",
     ylab="",xaxt="n",xlim=c(0.2,0.8),
     cex.lab=1, cex.axis=1, ylim=c(-10,142))

mtext(side = 2, "impact force (N)", line = 2, cex = 1)

axis(1, c(0.3,0.7), labels=c("dive","landing"),cex.axis=1)

dive.fem<-exp(bsim@fixef[,1])
#points(jitter(rep(1,2000),2), dive.fem, pch = 16,cex=0.5,col=alpha("#3080B5",0.1))

#dives
points(jitter(rep(0.28,length(all_indiv_land_complete_final_data$impact_force[which(all_indiv_land_complete_final_data$landing_cat=="dive" &
                                                                                      all_indiv_land_complete_final_data$sex =="Female")])),2.2),
       all_indiv_land_complete_final_data$impact_force[which(all_indiv_land_complete_final_data$landing_cat=="dive" &
                                                               all_indiv_land_complete_final_data$sex =="Female")], 
       pch = 16,
       cex=0.5,col=alpha("#3080B5",0.08))

points(jitter(rep(0.32,length(all_indiv_land_complete_final_data$impact_force[which(all_indiv_land_complete_final_data$landing_cat=="dive" &
                                                                                      all_indiv_land_complete_final_data$sex =="Male")])),2),
       all_indiv_land_complete_final_data$impact_force[which(all_indiv_land_complete_final_data$landing_cat=="dive" &
                                                               all_indiv_land_complete_final_data$sex =="Male")], 
       pch = 16,
       cex=0.5,col=alpha("#3080B5",0.08))

# landings

points(jitter(rep(0.68,length(all_indiv_land_complete_final_data$impact_force[which(all_indiv_land_complete_final_data$landing_cat=="landing" &
                                                                                      all_indiv_land_complete_final_data$sex =="Female")])),0.8),
       all_indiv_land_complete_final_data$impact_force[which(all_indiv_land_complete_final_data$landing_cat=="landing" &
                                                               all_indiv_land_complete_final_data$sex =="Female")], 
       pch = 16,
       cex=0.5,col=alpha("#F7C02D",0.08))

points(jitter(rep(0.72,length(all_indiv_land_complete_final_data$impact_force[which(all_indiv_land_complete_final_data$landing_cat=="landing" &
                                                                                      all_indiv_land_complete_final_data$sex =="Male")])),0.8),
       all_indiv_land_complete_final_data$impact_force[which(all_indiv_land_complete_final_data$landing_cat=="landing" &
                                                               all_indiv_land_complete_final_data$sex =="Male")], 
       pch = 16,
       cex=0.5,col=alpha("#F7C02D",0.08))

segments(x0 = c(0.28,0.32,0.68,0.72), y0 = exp(newdata_landingcat2$lwr), 
         x1 = c(0.28,0.32,0.68,0.72), y1 = exp(newdata_landingcat2$upr), 
         col = "white",
         lwd = 2)

points(c(0.28,0.32,0.68,0.72),t(exp(newdata_landingcat2$pred)),pch = 19,col="white", cex=0.5)

#points(sex.symbols(c(0.25,0.35,0.65,0.75), c(-.7,-.7,-.7,-.7), sex = c(2,1,2,1), col = par("fg"), lwd = par("lwd"), cex = 1))

#end
# install.packages("swfscMisc")
# library(swfscMisc)


#sex.symbols(c(0.28,0.32,0.68,0.72),c(140,140,140,140),sex=c(2,1,2,1), cex = 1.5)
#sex.symbols(c(0.28,0.32,0.68,0.72),c(0,0,0,0),sex=c(2,1,2,1), cex = 1.5)
sex.symbols(c(0.28,0.32,0.68,0.72),c(-5,-5,-5,-5),sex=c(2,1,2,1), cex = 1.5)

dev.off()




# GENERATION LAND TABLE FINAL ####
## extract time to hunt variable ####
names(all_indiv_land_complete_final_data)

lands<- all_indiv_land_complete_final_data %>% 
  ungroup() %>% 
  mutate(hunting_strategy_20s= case_when(hunt_strat_flyprop_20sec>0.1 ~ "flying",
                                         T ~ "perching"),
         hunting_strategy_0s= case_when(hunt_strat_flyprop_0sec==1 ~ "flying",
                                        T ~ "perching")) %>% 
  mutate(hunting_strategy_20s= as.factor(hunting_strategy_20s),
         hunting_strategy_0s= as.factor(hunting_strategy_0s)) %>% 
  filter(landing_cat=="landing") %>% 
  arrange(TagID,timestamp) %>% 
  group_by(TagID) %>% 
  mutate(ID = cur_group_id()) 


dives<-all_indiv_land_complete_final_data  %>% 
  ungroup() %>% 
  mutate(hunting_strategy_20s= case_when(hunt_strat_flyprop_20sec>0.1 ~ "flying",
                                         T ~ "perching"),
         hunting_strategy_0s= case_when(hunt_strat_flyprop_0sec==1 ~ "flying",
                                        T ~ "perching")) %>% 
  mutate(hunting_strategy_20s= as.factor(hunting_strategy_20s),
         hunting_strategy_0s= as.factor(hunting_strategy_0s)) %>% 
  filter(landing_cat=="dive") %>% 
  arrange(TagID,timestamp) %>% 
  group_by(TagID) %>% 
  mutate(ID = cur_group_id()) 




time_to_hunt<-rep(NA, nrow(lands))
time_to_hunt_units<- rep(NA, nrow(lands))
dist_to_hunt<-rep(NA, nrow(lands))



j=1

for(i in 1: nrow(lands)){
  print(paste(c("i=",i), sep=""))
  print(paste(c("j=",j), sep=""))
  # ?paste
  if(lands$TagID[i] == dives$TagID[j]){
    
    if(lands$timestamp[i] > dives$timestamp[j]){
      
      repeat {
        
        j=j+1 
        
        if(lands$TagID[i] != dives$TagID[j]){
          break
        }
        
        if(lands$timestamp[i] < dives$timestamp[j]) {
          
          time_to_hunt[i] =  dives$timestamp[j] - lands$timestamp[i]
          time_to_hunt_units[i] = units(dives$timestamp[j] - lands$timestamp[i])
          
          dist_to_hunt[i] = distm(c(lands$location.lon[i],lands$location.lat[i]),
                                  c(dives$location.lon[j],dives$location.lat[j]))
          
          break
        }
      }
      
    } else{
      
      time_to_hunt[i]=  dives$timestamp[j] - lands$timestamp[i]
      time_to_hunt_units[i] = units(dives$timestamp[j] - lands$timestamp[i])
      
      dist_to_hunt[i] = distm(c(lands$location.lon[i],lands$location.lat[i]),
                              c(dives$location.lon[j],dives$location.lat[j]))
    }
    
  } else if( lands$ID[i] < dives$ID[j]) {
    
    next
    
  } else {
    
    repeat{
      j=j+1 
      
      if(lands$TagID[i] == dives$TagID[j]){
        time_to_hunt[i]=  dives$timestamp[j] - lands$timestamp[i]
        time_to_hunt_units[i] = units(dives$timestamp[j] - lands$timestamp[i])
        
        dist_to_hunt[i] = distm(c(lands$location.lon[i],lands$location.lat[i]),
                                c(dives$location.lon[j],dives$location.lat[j]))
        break
      }
      
    }
  }
}

lands$time_to_hunt<- time_to_hunt
lands$time_to_hunt_units<- time_to_hunt_units
lands$dist_to_hunt<-dist_to_hunt
#put all time at the same scale (seconds)
lands$time_to_hunt_corrected<- NA
lands$time_to_hunt_corrected[which(lands$time_to_hunt_units=="secs")]= lands$time_to_hunt[which(lands$time_to_hunt_units=="secs")]
lands$time_to_hunt_corrected[which(lands$time_to_hunt_units=="mins")]= lands$time_to_hunt[which(lands$time_to_hunt_units=="mins")]*60
lands$time_to_hunt_corrected[which(lands$time_to_hunt_units=="hours")]= lands$time_to_hunt[which(lands$time_to_hunt_units=="hours")]*3600
lands$time_to_hunt_corrected[which(lands$time_to_hunt_units=="days")]= lands$time_to_hunt[which(lands$time_to_hunt_units=="days")]*86400


lands<-distinct(lands)

# MODEL Fit: GAMM landing force ####
names(lands2)

lands2<-lands  %>% 
  #filter(night_nb_roost!=0) %>% 
  #filter(time_to_hunt_corrected<30000) %>% 
  filter(time_to_hunt_corrected<18000) %>% # to exclude when they don't hunt in one night
  #filter(time_to_hunt_corrected<1800) %>% 
  filter(time_to_hunt_corrected>0) %>% 
  mutate(#time_to_hunt_log=log(time_to_hunt_corrected),
    has_prey=as.factor(has_prey)) %>% 
  filter(has_prey==0) %>% 
  mutate(TagID=as.factor(TagID),
         sex=as.factor(sex),
         RingId=as.factor(RingId)) %>% 
  filter(!is.na(time_to_hunt_corrected)) %>% 
  as_tibble() %>% 
  mutate(noise_pollution2=ifelse(is.na(noise_pollution),0,noise_pollution)) %>% 
  #filter(noise_pollution2<50) %>% 
  filter(max_vecsum<5) %>% ### warning check here
  filter(max_vecsum>1.5) %>% ### warning chheck here
  mutate(Year=as.factor(as.character(Year))) %>% 
  mutate(perch_type=case_when(distance_to_forest<2 | distance_to_hedges < 2 | distance_to_singletree < 2 ~ "tree",
                              distance_to_urban<2 ~ "urban" ,
                              # distance_to_highway <2 ~ "highway",
                              T ~"road_past")) %>% 
  mutate(perch_type= as.factor(perch_type)) 




# MODEL GAMM ####
## without sex interaction #### with windspeed k=9 ok --> without windspeed k9 ok # withhou perch type as linear predictor, k=8

# test representation raw data
#plot_dens<-
ggplot(sample_n(lands2, 3000), aes(log(time_to_hunt_corrected),impact_force))+
  geom_point()+
  geom_density_2d_filled(show.legend = F, adjust=0.9, alpha=0.7)+
  coord_cartesian(expand = F)+
  labs(x= "time to nenxt hunt (s)", y="vGRF (N)")


library(mgvc)
unique(lands2$sex)
lands2_nona<-lands2 %>% 
  drop_na(perch_type)

lands2$sex_f=as.factor(lands2$sex)
gamtest_allsex0<- gam(log(impact_force) ~ 
                        perch_type+
                        WindSpeed+
                        sex+
                        #perch_type/(1+WindSpeed)-1+
                        #sex+
                        s(time_to_hunt_corrected,by =perch_type,k=9)+
                        #s(Tag_night,bs="re") ,
                        s(TagID,bs="re"),## ajouter time to hunt dans "re" pour avoir un random slope
                      data = lands2,
                      #method = "REML",
                      family=gaussian(link="identity"))

summary(gamtest_allsex0)
plot(gamtest_allsex0)
table(lands2$perch_type)

tab_model(gamtest_allsex0, transform = "exp")

gam.check(gamtest_allsex0)

plot_model(gamtest_allsex0, type = "pred", terms = c("time_to_hunt_corrected[all]","perch_type"), show.data = F)+coord_cartesian(xlim = c(0, 10000), ylim = c(7,11))
plot_model(gamtest_allsex0, type = "pred", terms = c("time_to_hunt_corrected[all]","sex_f"), show.data = F)+
  coord_cartesian(xlim = c(0, 10000))
?plot_model
plot_model(gamtest_allsex0, type = "pred", terms = c("sex"), show.data = F)
##

options(na.action = "na.fail")
options(na.action = "na.omit")
dredge_gam=dredge(gamtest_allsex0, extra="R^2")
# so we take the model with all predictors !!!


library(modelsummary)
library(itsadug)
library(emmeans)
gamtabs(gamtest_allsex0)

modelsummary(gamtest_allsex0,exponentiate = TRUE, estimate = "{estimate} [{conf.low}, {conf.high}]",statistic = "{p.value}")


tab_model(gamtest_allsex0, string.est = "Estimate", 
          string.ci = "Conf.Int.", string.p = "p-value",transform = "exp")

emmeanns<-emmeans(gamtest_allsex0, pairwise~ sex , type="resp",nesting = NULL)

emmeanns<-emmeans(gamtest_allsex0, pairwise~   perch_type,
                  type="resp",nesting = NULL, at=c(0,0.9))


library(suncalc)

####
install.packages("modelbased")
install.packages("datawizard")
unloadNamespace("bayestestR")
library(datawizard)
library(modelbased)
unload


summary(gamtest_allsex0)
summary(gamtest_allsex_norand)
summary(gamtest_allsex2)
summary(gamtest_allsex3)
AIC(gamtest_allsex)
AIC(gamtest_allsex2)
AIC(gamtest_allsex3)
anova.gam(gamtest_allsex,gamtest_allsex2,gamtest_allsex3)
gam.check(gamtest_allsex0)

plot.gam(gamtest_allsex0)
plot(gamtest_allsex0)
simulationOutput <- simulateResiduals(fittedModel = gamtest_allsex0, plot = T)
plot_model(gamtest_allsex0, type = "diag")

y_lim<-c(9,11.5)
x_lim<-c(0,7200)
plot_model(gamtest_allsex0, type = "pred", terms = c("time_to_hunt_corrected[all]","perch_type"), show.data = F, axis.lim = list(x_lim,y_lim))
plot_model(gamtest_allsex0, type = "pred", terms = c("time_to_hunt_corrected[all]"), show.data = F)
plot_model(gamtest_allsex0, type = "emm", terms = c("sex"), show.data = F,nesting = NULL)
plot_model(gamtest_allsex0, type = "pred", terms = c("WindSpeed","perch_type"), show.data = F)

plot_model(gamtest_allsex, type = "pred", terms = c("perch_type"))

lands2$dist_to_hunt

### derivative of fthe gam ####

library(gratia)
## perches type ## urban
der_urban<-fderiv(gamtest_allsex0)
fder_urban<-der_urban$derivatives$`s(time_to_hunt_corrected):perch_typeurban`$deriv
time_to_hunt_fit_urban<-der_urban$eval$time_to_hunt_corrected

se.deriv_urban<- der_urban$derivatives$`s(time_to_hunt_corrected):perch_typeurban`$se.deriv
upper_urban<- fder_urban + (1.96* se.deriv_urban)
lower_urban<- fder_urban - (1.96* se.deriv_urban)


inflex_pt_urban<- time_to_hunt_fit_urban[which.max(lower_urban<=0)]
inflex_pos_urban<- which.max(lower_urban<0)
par(mfrow=c(1,1))

plot(fder_urban~time_to_hunt_fit_urban, type="l", col="red" ,ylim=c(-0.00003,0.00004),xlim=rev(range(1,5400)), # focus on the last 90 min
     
     ylab="first derivative",
     xlab="Time to the next hunt (secs)")
lines(time_to_hunt_fit_urban, upper_urban, lty = 'dashed', col = 'grey')
lines(time_to_hunt_fit_urban, lower_urban, lty = 'dashed', col = 'grey')
abline(0,0, col="black")

time.signif.urban_increase<-time_to_hunt_fit_urban[which(upper_urban>0 & lower_urban>0)]
time.signif.urban_decrease<-time_to_hunt_fit_urban[which(upper_urban<0 & lower_urban<0)]

#dev.off()

####################
####################

### first derivative - tsgam ####
## parameters for testing


install.packages("gratia")
library(gratia)

UNCONDITIONAL <- FALSE # unconditional or conditional on estimating smooth params?
N <- 10000             # number of posterior draws
n <- 1000               # number of newdata values
EPS <- 1e-07           # finite difference

# first derivative calculation

# pasture pole

gamtest_all_pasture<- gam(log(impact_force) ~ 
                            perch_type+
                            WindSpeed+
                            #perch_type/(1+WindSpeed)-1+
                            sex+
                            s(time_to_hunt_corrected,by=perch_type,k=9)+
                            s(TagID,bs="re") , ## ajouter time to hunt dans "re" pour avoir un random slope
                          data = lands2[-1,],
                          # method = "REML",
                          family=gaussian(link="identity"))

summary(gamtest_all_pasture)
summary(gamtest_allsex0)
der_pasture<-fderiv(gamtest_all_pasture, n=200)
#?fderiv
fder_pasture<-der_pasture$derivatives$`s(time_to_hunt_corrected):perch_typeroad_past`$deriv
time_to_hunt_fit_pasture<-der_pasture$eval$time_to_hunt_corrected

se.deriv_pasture<- der_pasture$derivatives$`s(time_to_hunt_corrected):perch_typeroad_past`$se.deriv
upper_pasture<- fder_pasture + (1.96* se.deriv_pasture)
lower_pasture<- fder_pasture - (1.96* se.deriv_pasture)

inflex_pt_pasture<- time_to_hunt_fit_pasture[which.max(lower_pasture<=0)]
inflex_pos_pasture<- which.max(lower_pasture<0)
par(mfrow=c(1,1))
abline(0,0, col="black")


pos.signif.pasture_increase<-which(upper_pasture>0 & lower_pasture>0)
der.signif.pasture_increase<-fder_pasture[,1][which(upper_pasture>0 & lower_pasture>0)]

time.signif.pasture_increase<-time_to_hunt_fit_pasture[which(upper_pasture>0 & lower_pasture>0)]
time.signif.pasture_decrease<-time_to_hunt_fit_pasture[which(upper_pasture<0 & lower_pasture<0)]

significant_pasture<-rep(NA, 200)
significant_pasture[pos.signif.pasture_increase]<-fder_pasture[which(upper_pasture>0 & lower_pasture>0)]


plot(fder_pasture~time_to_hunt_fit_pasture, type="l", col="grey" ,ylim=c(-0.00003,0.00004),xlim=rev(range(1,5400)), # for te focus on the last 90 min
     
     ylab="first derivative",
     xlab="Time to the next hunt (secs)")
abline(0,0, col="black")
lines(time_to_hunt_fit_pasture, upper_pasture, lty = 'dashed', col = 'grey')
lines(time_to_hunt_fit_pasture, lower_pasture, lty = 'dashed', col = 'grey')
lines(time_to_hunt_fit_pasture,significant_pasture, col = 'red')   




# trees

gamtest_all_tree<- gam(log(impact_force) ~ 
                         perch_type+
                         WindSpeed+
                         #perch_type/(1+WindSpeed)-1+
                         sex+
                         s(time_to_hunt_corrected,by=perch_type,k=9)+
                         s(TagID,bs="re") , ## ajouter time to hunt dans "re" pour avoir un random slope
                       data = lands2[-c(1:2),],
                       # method = "REML",
                       family=gaussian(link="identity"))
summary(gamtest_all_tree)

der_tree<-fderiv(gamtest_all_tree, n = 100)
fder_tree<-der_tree$derivatives$`s(time_to_hunt_corrected):perch_typetree`$deriv
time_to_hunt_fit_tree<-der_tree$eval$time_to_hunt_corrected

se.deriv_tree<- der_tree$derivatives$`s(time_to_hunt_corrected):perch_typetree`$se.deriv
upper_tree<- fder_tree + (1.96* se.deriv_tree)
lower_tree<- fder_tree - (1.96* se.deriv_tree)

inflex_pt_tree<- time_to_hunt_fit_tree[which.max(lower_tree<=0)]
inflex_pos_tree<- which(lower_tree>0)
par(mfrow=c(1,1))

pos.signif.tree_increase<-which(upper_tree>0 & lower_tree>0)
der.signif.tree_increase<-fder_tree[,1][which(upper_tree>0 & lower_tree>0)]

pos.signif.tree_decrease<-which(upper_tree<0 & lower_tree<0)
der.signif.tree_decrease<-fder_tree[,1][which(upper_tree<0 & lower_tree<0)]


time.signif.tree_increase<-time_to_hunt_fit_tree[which(upper_tree>0 & lower_tree>0)]
time.signif.tree_decrease<-time_to_hunt_fit_tree[which(upper_tree<0 & lower_tree<0)]

significant_tree_increase<-rep(NA, 100)
significant_tree_increase[pos.signif.tree_increase]<-fder_tree[which(upper_tree>0 & lower_tree>0)]

significant_tree_decrease<-rep(NA, 100)
significant_tree_decrease[pos.signif.tree_decrease]<-fder_tree[which(upper_tree<0 & lower_tree<0)]

#pdf(file="/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/derivative_tree.pdf")
pdf(file="/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/derivative_tree_zoom.pdf")

#plot(fder_tree~time_to_hunt_fit_tree, type="l", col="grey" , ylim=c(-0.00003,0.00004),xlim=rev(range(time_to_hunt_fit_tree)),#full range
plot(fder_tree~time_to_hunt_fit_tree, type="l", col="grey" , ylim=c(-0.00003,0.00004),xlim=rev(range(1,5400)),     # focus last 90 min
     ylab="first derivative",
     xlab="Time to the next hunt (secs)")
abline(0,0, col="black")
lines(time_to_hunt_fit_tree, upper_tree, lty = 'dashed', col = 'grey')
lines(time_to_hunt_fit_tree, lower_tree, lty = 'dashed', col = 'grey')
lines(time_to_hunt_fit_tree,significant_tree_increase, col = 'red')   
lines(time_to_hunt_fit_tree,significant_tree_decrease, col = 'blue')   

## plot final ####


plot_mod_gam2_response_all = ggplot(aes(x=(time_to_hunt_corrected*-1),y= exp_fit), group=perch_type, data=predicts[predicts$sex=="Male",]) + # geom_point(aes(x=time_to_hunt_corrected_min, y= impact_force) ,data=lands2_min_male)+
  geom_ribbon(aes(group=perch_type,fill=perch_type,ymin = exp(lower), ymax=exp(upper)), alpha=0.3) +
  #scale_fill_manual(values = c( "#faaa20" ,"#2F7D3B","#FFD027"))+
  scale_fill_manual(values = c( "#959083" ,"#2F7D3B","#FCBC03"))+ # best pour l'instant
  #scale_fill_manual(values = c( "#E3D29F" ,"#2F7D3B","#FFD027"))+
  # geom_line(aes(group=perch_type,colour= significant)) +
  # scale_colour_manual(values=c("grey", "white"))+
  geom_line(aes(group=perch_type)) +
  scale_colour_manual(values= "red")+
  
  scale_x_continuous(limits = c(-5400,0), breaks = c(0,-1800,-3600,-5400) ,
                     labels =c("0","-30","-60", "-90"))+

  scale_y_continuous(breaks = seq(6, 13, by = 1))+

  ylim(7.2,9)+ 
 
  ggtitle("Predicted values of impact_force")+
    labs(y= "Landing force (N)", x = "Pre-hunt time (min)")+
  
 
  theme_bw() + theme( # dans les autres cas
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18))
plot_mod_gam2_response_all
figpath<-"/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/"
ggsave(paste0(figpath, "Figure_GAM_68.pdf"), plot = plot_mod_gam2_response_all, width = 7, height = 6, units = "in", dpi = 300) 




xdensity <- ggplot(lands2[lands2$sex=="Female",], aes(x=((time_to_hunt_corrected)*-1), fill=perch_type)) + 
  geom_histogram(aes(y=..density..), alpha=.85)+
  #geom_histogram(aes(y = after_stat(count)), alpha=.85)+
  # geom_density(alpha=.3,) + 
  #facet_wrap(~sex, ncol = 1)+
  scale_fill_manual(values = c("#959083" ,"#2F7D3B","#FCBC03")) + 
  theme(legend.position = "none")+
  scale_x_continuous(limits = c(-2000,0))+
  theme_bw() + theme( # dans les autres cas
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18))

xdensity

figpath<-"/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/"
ggsave(paste0(figpath, "Figure_perchtype_density.pdf"), plot = xdensity , width = 6, height = 6, units = "in", dpi = 300) 

##
ydensity <- ggplot(lands2[lands2$sex=="Female",], aes(x=log(impact_force), fill=perch_type)) + 
  geom_density(alpha=.3) + 
  scale_fill_manual(values = c("#959083" ,"#2F7D3B","#FCBC03")) + 
  theme(legend.position = "none")+
  scale_x_continuous(limits = c(1,3))+
  theme_bw() + theme( # dans les autres cas
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18))

ydensity




ydensity <- ggplot(lands2, aes(x=impact_force, fill=perch_type)) + 
  geom_boxplot(alpha=.3) + 
  facet_wrap(~sex, ncol=1)+
  scale_fill_manual(values = c("#959083" ,"#2F7D3B","#FCBC03")) + 
  theme(legend.position = "none")+
  #scale_x_continuous(limits = c(1,3))+
  theme_bw() + theme( # dans les autres cas
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18))

ydensity




# all together

library("gridExtra")

blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
  )

grid.arrange(xdensity, 
             #blankPlot, 
             plot_mod_gam2_response_all, 
             #“blankPlot, 
             ncol=2, nrow=1)
group_by(perch_type) 

table(lands_perchetype$perch_type)
str(lands_perchetype$perch_type)

lands_perchetype$noise_pollution
table(lands_perchetype$perch_type)

prop_closetohunt<-lands2 %>% 
  filter(time_to_hunt_corrected<10) %>% 
  group_by(perch_type) %>% 
  summarise(n=n()) %>% 
  mutate(prop=n/sum(n))

library(sjstats)

prop_ffartohunt<-lands2 %>% 
  group_by(perch_type) %>% 
  filter(time_to_hunt_corrected>1800) %>% 
  summarise(n=n()) %>% 
  mutate(prop=n/sum(n))

p3 <- ggplot(data=lands2, aes(x=time_to_hunt_corrected,
                              group=interaction(perch_type, Year) , 
                              fill=perch_type, shape=Year)) +
  geom_density(adjust=3, alpha=.4) +
  xlim(0,1000)+
  theme_ipsum()



# Hunt strike success ####
##  hunting success according to the landing force ####



tbl_timetohunt<-lands2 %>% 
  dplyr::select(TagID,landing_cat,landing_id,time_to_hunt_corrected,dist_to_hunt,perch_type)



test<- tbl_timetohunt %>% 
  as_tibble() %>% 
  right_join(all_indiv_land_complete_final_data) %>% 
  group_by(TagID, landing_cat, landing_id) %>% 
  arrange(TagID,landing_id)



#calculate average impact force over the n previous landings ## add distance to the next hunt?
library(geosphere)

landing_force_previous_landing<-rep(NA, length(test$TagID)) 
time_to_hunt_from_previous_landing<-rep(NA, length(test$TagID))
dist_from_previous_landing<- rep(NA,length(test$TagID))
perch_type_previous_landing<-rep(NA,length(test$TagID))

i=1

for(i in 2:nrow(test)){
  
  if(test$TagID[i]!=test$TagID[i-1]){
    next
  } else if(test$landing_cat[i]=="landing"){
    next
  }else if( test$landing_cat[i]=="dive"){
    
    if (test$landing_cat[i-1]=="dive"){
      next
    }
    
    landing_force_previous_landing[i]= test$impact_force[i-1]
    time_to_hunt_from_previous_landing[i]=test$time_to_hunt_corrected[i-1]
    dist_from_previous_landing[i]=  abs(distm(c(test$location.lon[i],test$location.lat[i]),
                                              c(test$location.lon[i-1],test$location.lat[i-1])))
    perch_type_previous_landing[i] = as.character(test$perch_type[i-1])
  }
}


dive_tbl<-test %>% 
  as_tibble() %>% 
  filter(landing_cat2!="landing") %>% 

  mutate(dive_bmk_sucess= as.factor(ceiling(as.numeric(as.character(dive_bmk_sucess))))) %>% 
  mutate(sex=as.factor(sex)) %>% 
  mutate(TagID=as.factor(TagID)) %>% 
  #filter(time_to_hunt_from_previous_landing>0) %>% 
  mutate(noise_pollution2=ifelse(is.na(noise_pollution),0,noise_pollution)) %>% 
  #mutate(hunting_strategy_20s=as.factor(hunting_strategy_20s)) %>% 
  mutate(hunting_strategy_20s_05=case_when(hunt_strat_flyprop_20sec>=0.5~"flying",
                                           T ~"perching")) %>% 
  mutate(hunting_strategy_20s_0_1=case_when(hunt_strat_flyprop_20sec==1~"flying",
                                            hunt_strat_flyprop_20sec==0~"perching",
                                            T ~NA_character_)) %>% 
  mutate(hunting_strategy_5s_0_1=case_when(hunt_strat_flyprop_5sec==1~"flying",
                                           hunt_strat_flyprop_5sec==0~"perching",
                                           T ~NA_character_)) %>% 
  mutate(hunting_strategy_5s_05=case_when(hunt_strat_flyprop_5sec==1~"flying",
                                          hunt_strat_flyprop_0sec==1.~"perching",
                                          T ~ NA_character_)) %>% 
  mutate(hunting_strategy_final=case_when(hunt_strat_flyprop_5sec==1~"flying",
                                          hunt_strat_flyprop_0sec==0~"perching",
                                          T ~ NA_character_)) %>% 
  mutate(hunting_strategy_final2=case_when(hunt_strat_flyprop_5sec==1~"flying",
                                           hunt_strat_flyprop_5sec==0~"perching",
                                           T ~ "hybrid")) %>% 
  mutate(hunting_strategy_20s_05= as.factor(hunting_strategy_20s_05)) %>% 
  mutate(hunting_strategy_20s_0_1= as.factor(hunting_strategy_20s_0_1)) %>% 
  mutate(hunting_strategy_final= as.character(hunting_strategy_final)) %>% 
  #filter(hunting_strategy_final=="perching") %>% 
  filter(time_to_hunt_from_previous_landing<90) %>%  # %>% when tripple interaction, fixed to 500 other wise select only time to hunt within 60 sec
  mutate(Genus= ifelse(Genus == "", NA, Genus)) %>% 
  drop_na(WindSpeed)
#filter(sex=="Male")


#### dive_tbl_all
dive_tbl_full<-test %>% 
  as_tibble() %>% 
  filter(landing_cat2!="landing") %>% 
  #filter(max_vecsum<9) %>% ### warning check here
  #filter(max_vecsum>1.5) %>% ### warning chheck here
  #filter(hunting_strategy_20s=="perching") %>% 
  #filter(landing_force_previous_landing<30) %>% ## tripple fixed to 40
  mutate(dive_bmk_sucess= as.factor(ceiling(as.numeric(as.character(dive_bmk_sucess))))) %>% 
  mutate(sex=as.factor(sex)) %>% 
  mutate(TagID=as.factor(TagID)) %>% 
  #filter(time_to_hunt_from_previous_landing>0) %>% 
  mutate(noise_pollution2=ifelse(is.na(noise_pollution),0,noise_pollution)) %>% 
  #mutate(hunting_strategy_20s=as.factor(hunting_strategy_20s)) %>% 
  mutate(hunting_strategy_20s_05=case_when(hunt_strat_flyprop_20sec>=0.5~"flying",
                                           T ~"perching")) %>% 
  mutate(hunting_strategy_20s_0_1=case_when(hunt_strat_flyprop_20sec==1~"flying",
                                            hunt_strat_flyprop_20sec==0~"perching",
                                            T ~NA_character_)) %>% 
  mutate(hunting_strategy_5s_0_1=case_when(hunt_strat_flyprop_5sec==1~"flying",
                                           hunt_strat_flyprop_5sec==0~"perching",
                                           T ~NA_character_)) %>% 
  mutate(hunting_strategy_5s_05=case_when(hunt_strat_flyprop_5sec==1~"flying",
                                          hunt_strat_flyprop_0sec==1.~"perching",
                                          T ~ NA_character_)) %>% 
  mutate(hunting_strategy_final=case_when(hunt_strat_flyprop_5sec==1~"flying",
                                          hunt_strat_flyprop_0sec==0~"perching",
                                          T ~ NA_character_)) %>% 
  mutate(hunting_strategy_final2=case_when(hunt_strat_flyprop_5sec==1~"flying",
                                           hunt_strat_flyprop_5sec==0~"perching",
                                           T ~ "hybrid")) %>% 
  mutate(hunting_strategy_20s_05= as.factor(hunting_strategy_20s_05)) %>% 
  mutate(hunting_strategy_20s_0_1= as.factor(hunting_strategy_20s_0_1)) %>% 
  mutate(hunting_strategy_final= as.character(hunting_strategy_final)) %>% 
  #filter(hunting_strategy_final=="perching") %>% 
  #filter(time_to_hunt_from_previous_landing<90) %>%  # %>% when tripple interaction, fixed to 500 other wise select only time to hunt within 60 sec
  mutate(Genus= ifelse(Genus == "", NA, Genus))


######
## model hhunting success ####
names(dive_tbl)
str(dive_tbl)
### test 

dive_tbl$sex_relevel<-relevel(dive_tbl$sex,"Male")


library(lme4)
library(lmerTest)
library(performance)
# check<-dive_tbl %>% 
#   distinct()


dive_tbl<-dive_tbl %>% 
  drop_na((WindSpeed))

glm_test<-glmer(dive_bmk_sucess ~ 
                  #scale(landing_force_previous_landing)+
                  ##scale(I(landing_force_previous_landing^2))+
                  #hunting_strategy_final*
                  #sex+
                  sex*hunting_strategy_final+
                  #scale(time_to_hunt_from_previous_landing)+
                  #scale(dist_from_previous_landing)+ # check savoir si je l'ai lassé a la fin ou pas
                  ##scale(I(dist_from_previous_landing^2))+
                  #scale(WindSpeed)+
                  ##scale(impact_force)+
                  hunting_strategy_final/(1+scale(landing_force_previous_landing))-1+
                  #hunting_strategy_5s_0_1/(1+scale(landing_force_previous_landing))-1+ # this gives same result but better effect
                  #perch_type_previous_landing+
                  ##scale(dist_from_previous_landing):scale(landing_force_previous_landing)+
                  (1|TagID/night_nb_real),
                family =binomial(link = "logit"), data = dive_tbl)

summary(glm_test)


# sex comp hunntin success ####
try<-glmer(dive_bmk_sucess ~ 
             sex*hunting_strategy_final+
             #hunting_strategy_final/(1+sex) -1 +
             (1|TagID/night_nb_real),
           family =binomial(link = "logit"), data = dive_tbl_full)
check_model(try)
summary(try)



emm.test<-emmeans(try, pairwise~ hunting_strategy_final*sex, 
                  type="response")

dat.emm.test<-emm.test$emmeans


emm.test<-emmeans(try, pairwise~ hunting_strategy_final, 
                  type="response")


tab_model(glm_test)

glm_test<-glmer(dive_bmk_sucess ~ 
                  #scale(landing_force_previous_landing)*
                  ##scale(I(landing_force_previous_landing^2))+
                  #hunting_strategy_final*
                  sex*hunting_strategy_final+
                  #scale(time_to_hunt_from_previous_landing)+
                  #scale(dist_from_previous_landing)+ # check savoir si je l'ai lassé a la fin ou pas
                  ##scale(I(dist_from_previous_landing^2))+
                  #scale(WindSpeed)+
                  ##scale(impact_force)+
                  
                  hunting_strategy_final*scale(landing_force_previous_landing)+
                  #sex*hunting_strategy_final+
                  #hunting_strategy_5s_0_1/(1+scale(landing_force_previous_landing))-1+ # this gives same result but better effect
                  #perch_type_previous_landing+
                  ##scale(dist_from_previous_landing):scale(landing_force_previous_landing)+
                  (1|TagID/night_nb_real),
                family =binomial(link = "logit"), data = dive_tbl_no_na)


summary(glm_test)
dive_tbl$perch_type_previous_landing
### dredge
dive_tbl_no_na<-dive_tbl %>% 
  dplyr::select(sex,WindSpeed,hunting_strategy_final,landing_force_previous_landing,TagID,night_nb_real, dive_bmk_sucess, perch_type_previous_landing,time_to_hunt_from_previous_landing) %>% 
  drop_na()

mod_fit_dredge<-glm_test<-glmer(dive_bmk_sucess ~ 
                                  #scale(landing_force_previous_landing)+
                                  ##scale(I(landing_force_previous_landing^2))+
                                  #hunting_strategy_final+
                                  # scale(time_to_hunt_from_previous_landing)+
                                  #scale(dist_from_previous_landing)+ # check savoir si je l'ai lassé a la fin ou pas
                                  ##scale(I(dist_from_previous_landing^2))+
                                  ##scale(impact_force)+
                                  #sex+
                                  sex*hunting_strategy_final*scale(landing_force_previous_landing)+
                                  ###hunting_strategy_final*scale(landing_force_previous_landing)+
                                  ### hunting_strategy_final/(1+scale(landing_force_previous_landing))-1+
                                  #hunting_strategy_final*scale(landing_force_previous_landing)+#*sex+
                                  # perch_type_previous_landing+
                                  scale(WindSpeed)+
                                  #scale(landing_force_previous_landing)*hunting_strategy_final+# remove 
                                  #hunting_strategy_5s_0_1/(1+scale(landing_force_previous_landing))-1+ # this gives same result but better effect
                                  perch_type_previous_landing+
                                  ##scale(dist_from_previous_landing):scale(landing_force_previous_landing)+
                                  (1|TagID/night_nb_real),
                                family =binomial(link = "logit"), data = dive_tbl_no_na)


summary(mod_fit_dredge)


library(MuMIn)
options(na.action = "na.fail")

options(na.action = "na.omit")
dredge<- dredge(mod_fit_dredge, extra="R^2" )
dredge


AIC(glm_test)
summary(glm_test)
check_model(glm)

library(sjPlot)
str(dive_tbl$hunting_strategy_final)

plot_model(mod_fit_dredge, terms=c("landing_force_previous_landing","perch_type_previous_landing"), type="emm")

### check how much data are dicarded with NA hunting strategy
dive_tbl_check<-dive_tbl %>% 
  dplyr::select(dive_bmk_sucess,hunting_strategy_5s_05,landing_force_previous_landing,WindSpeed,sex) %>% 
  drop_na(landing_force_previous_landing) %>% 
  drop_na(WindSpeed) %>% 
  drop_na(sex) %>% 
  drop_na(dive_bmk_sucess)


##prediction

newdat<-expand.grid(sex=c("Male","Female"),
                    WindSpeed=mean(dive_tbl$WindSpeed),
                    hunting_strategy_final(c("flying","perching")),
)


# plot 
library(plyr)
library(ggpointdensity)


dive_tbl_plot<-dive_tbl_no_na %>% 
  mutate(new_resp=round_any(landing_force_previous_landing,0.5)) %>% 
  mutate(new_resp=new_resp) %>% 
  drop_na(hunting_strategy_final) 
#drop_na(hunting_strategy_5s_0_1) 


table_bubble<- as.data.frame(table(dive_tbl_plot$new_resp,dive_tbl_plot$hunting_strategy_final, dive_tbl_plot$dive_bmk_sucess)) 
str(table_bubble)
colnames(table_bubble)<- c("new_resp", "hunting_strategy_final", "success" ,"n")
#colnames(table_bubble)<- c("new_resp", "hunting_strategy_5s_0_1", "success" ,"n")

bubble_fly<-table_bubble %>% 
  filter(hunting_strategy_final=="flying") %>% 
  #filter(hunting_strategy_5s_0_1=="flying") %>% 
  mutate(new_resp=as.numeric(as.character(new_resp)))


bubble_perch<-table_bubble %>% 
  filter(hunting_strategy_final=="perching") %>% 
  #filter(hunting_strategy_5s_0_1=="perching") %>% 
  mutate(new_resp=as.numeric(as.character(new_resp)))

str(bubble_fly)


df<-ggeffect(mod_fit_dredge,terms = c("landing_force_previous_landing[all]","hunting_strategy_final"))
#df<-ggemmeans(glm_test,terms = c("landing_force_previous_landing[all]","hunting_strategy_5s_0_1"))

plot_hunt_succ<-ggplot(df, aes(x, predicted))+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, colour=group), alpha=0.08) +
  geom_line(aes(linetype=group, colour=group)) +
  #scale_color_manual(values=c('#F6C987','#2D6E7F'))+
  scale_color_manual(values=c('#EF6B22','#29686B'))+
  
  geom_point(data = bubble_fly,
             #colour='#F6C987',
             colour='#EF6B22',
             alpha=0.15,
             aes(x= as.numeric(new_resp), 
                 y=as.numeric(success)-0.98,
                 size=n)) +
  geom_point(data = bubble_perch,
             #colour='#2D6E7F',
             colour='#29686B',
             alpha=0.15,
             aes(x= as.numeric(new_resp), 
                 y=as.numeric(success)-1.02,
                 size=n))+
  scale_size_binned(range = c(3, 16),
                    breaks = breaks_extended(n=6),
                    limits = c(1,150))+
  
  xlab("pre-hunt landing force (N)") + 
  ylab("Hunnting success (%)")+
  scale_y_continuous(labels = percent)+
  theme_bw() + theme(panel.border = element_blank(), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     axis.title.x = element_text(size = 22),
                     axis.title.y = element_text(size = 22),
                     axis.text.x = element_text(size = 18),
                     axis.text.y = element_text(size = 18))
plot_hunt_succ

# Difference in distance between last land and hunt strike according to both hunting strategy ####

table_dist<-dive_tbl %>% 
  dplyr::select(hunting_strategy_final,dist_from_previous_landing, WindSpeed, landing_force_previous_landing) %>% 
  na.omit() %>% 
  group_by(hunting_strategy_final) %>%
  summarise(
    count = n(),
    median = median(dist_from_previous_landing, na.rm = TRUE),
    sd = sd(dist_from_previous_landing, na.rm = TRUE)) 


plot_dist<-ggboxplot(dive_tbl %>% filter(!is.na(hunting_strategy_final)), x = "hunting_strategy_final", y = "dist_from_previous_landing", 
                     color = "hunting_strategy_final", palette = c('#EF6B22','#29686B' ),alpha=0.15,
                     fill= "hunting_strategy_final")+
  labs(x="Hunting strategy", y="Distance from pre-hunt perch to strike point (m)")+
  #stat_compare_means(method = "wilcox.test",hide.ns = FALSE)+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")


# MODEL FIT DIVES ####

table_dives<-all_indiv_land_complete_final_data  %>% 
  ungroup() %>% 
  mutate(hunting_strategy_20s= case_when(hunt_strat_flyprop_20sec>0.1 ~ "flying",
                                         T ~ "perching"),
         hunting_strategy_0s= case_when(hunt_strat_flyprop_0sec==1 ~ "flying",
                                        T ~ "perching")) %>% 
  mutate(hunting_strategy_20s= as.factor(hunting_strategy_20s),
         hunting_strategy_0s= as.factor(hunting_strategy_0s)) %>% 
  mutate(hunting_strategy_final=case_when(hunt_strat_flyprop_5sec==1~"flying",
                                          hunt_strat_flyprop_0sec==0~"perching",
                                          T ~ NA_character_)) %>% 
  mutate(hunting_strategy_5s_0_1=case_when(hunt_strat_flyprop_5sec==1~"flying",
                                           hunt_strat_flyprop_5sec==0~"perching",
                                           T ~NA_character_)) %>% 
  filter(landing_cat=="dive") %>% 
  arrange(TagID,timestamp) %>% 
  group_by(TagID) 


##### modelling variation on dives impact
plot(check_distribution(log(table_dives$impact_force)))
names(table_dives)
table_dives$f

model_strike_force_divesucc<- lmer(log(impact_force) ~ landing_cat2+
                                     sex+
                                     hunting_strategy_final+
                                     
                                     #hunting_strategy_5s_0_1+
                                     #WindSpeed+
                                     # landing_cat2:sex + 
                                     #sex:hunting_strategy_final + 
                                     #sex:landing_cat2+
                                     
                                     landing_cat2:hunting_strategy_final+
                                     #landing_cat2:sex+
                                     #hunting_strategy_final:sex+
                                     #landing_cat2:hunting_strategy_5s_0_1+
                                     (1|TagID/night_nb_real),
                                   data = table_dives, REML = F)
# ,
#                                                        control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e+06)))

dredge_strikeforce<- dredge(model_strike_force_divesucc)


tab_model(model_strike_force_divesucc ,string.est = "Estimate",show.fstat = T,
          string.ci = "Conf. Int.", string.p = "p-value",digits=2, transform = "exp")


tab_model(model_strike_force_divesucc  ,string.est = "Estimate", show.df = F,show.fstat = T,
          string.ci = "Conf. Int.", string.p = "p-value",digits=2,transform = "exp",file = "/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/table_dives_succ_lastresult.html" )

m.emm_succes<-emmeans(model_strike_force_divesucc, pairwise~ sex, 
                      type="response", transform="log" )

m.emm_succes<-emmeans(model_strike_force_divesucc, pairwise~ landing_cat2*hunting_strategy_final, 
                      type="response", transform="log" )


check_model(model_strike_force_divesucc)
r.squaredGLMM(model_strike_force_divesucc)
m.emm_dives<-emmeans(model_strike_force_divesucc, pairwise~landing_cat2*sex, 
                     transform="log", type="response" )
summary(m.emm_dives)
plot(m.emm_dives)

m.emm_dives<-emmeans(model_strike_force_divesucc, pairwise~landing_cat2*hunting_strategy_final, 
                     transform="log", type="response" )
summary(m.emm_dives)



# plot final ####

nsim <- 2000
bsim <- sim(model_strike_force_divesucc, n.sim=nsim)
str(bsim)
colnames (bsim @ fixef) <- names(fixef(model_strike_force_divesucc)) # for arm 1.6e10
fixef(model_strike_force_divesucc)
apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))


summary(model_strike_force_divesucc)


## get model prediction
newdata_impactcat<- expand.grid (landing_cat2 = c("dive_successfull","dive_unsuccessfull"),
                                 hunting_strategy_final=c("perching","flying"),
                                 sex=c("Female", "Male"))


Xmat_impactcat <- model.matrix(~landing_cat2+sex+ hunting_strategy_final+landing_cat2:hunting_strategy_final, data=newdata_impactcat)   # the model matrix

newdata_impactcat$pred <- Xmat_impactcat%*%fixef(model_strike_force_divesucc)  

predmat_impactcat <- matrix(nrow=nrow(newdata_impactcat),ncol=nsim)
for(i in 1:nsim) predmat_impactcat[,i] <-Xmat_impactcat  %*% bsim@fixef[i,]
newdata_impactcat$lwr <- apply(predmat_impactcat,1,quantile,probs=0.025)
newdata_impactcat$upr <- apply(predmat_impactcat,1,quantile,probs=0.975)


# plot 
newdata_impactcat2<-newdata_impactcat %>% 
  arrange(landing_cat2)

#summary(log(all_indiv_land_complete_final_data$impact_force))

pdf(file="/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/Figure_dive_succ_lastresult_zoomed.pdf")

plot(as.numeric(newdata_impactcat2$landing_cat2),newdata_impactcat2$pred,type="n",las=1,xlab="Hunting Strategy",
     ylab="",xaxt="n",xlim=c(0.2,0.8),
     cex.lab=2, cex.axis=1, ylim=c(37,43))
#cex.lab=2, cex.axis=1, ylim=c(0,143))
mtext(side = 2, "impact force (N)", line = 2.5, cex = 2)

axis(1, c(0.3,0.7), labels=c("Flying","Perching"),cex.axis=1.6)

points(jitter(rep(0.28,length(table_dives$impact_force[which(table_dives$landing_cat2=="dive_successful" &
                                                               table_dives$hunting_strategy_final =="flying"&
                                                               table_dives$sex=="Female")])),2.2),
       table_dives$impact_force[which(table_dives$landing_cat2=="dive_successful" &
                                        table_dives$hunting_strategy_final =="flying"&
                                        table_dives$sex=="Female")], 
       pch = 16,
       cex=0.5,col=alpha("#3080B5",0.1))

points(jitter(rep(0.32,length(table_dives$impact_force[which(table_dives$landing_cat2=="dive_unsuccessful" &
                                                               table_dives$hunting_strategy_final =="flying"&
                                                               table_dives$sex=="Female")])),2),
       table_dives$impact_force[which(table_dives$landing_cat2=="dive_unsuccessful" &
                                        table_dives$hunting_strategy_final =="flying"&
                                        table_dives$sex=="Female")], 
       pch = 16,
       cex=0.5,col=alpha("#3080B5",0.1))

# perching

points(jitter(rep(0.68,length(table_dives$impact_force[which(table_dives$landing_cat2=="dive_successful" &
                                                               table_dives$hunting_strategy_final =="perching" &
                                                               table_dives$sex=="Female")])),1),
       table_dives$impact_force[which(table_dives$landing_cat2=="dive_successful" &
                                        table_dives$hunting_strategy_final =="perching" &
                                        table_dives$sex=="Female")], 
       pch = 16,
       cex=0.5,col=alpha("#F7C02D",0.5))



points(jitter(rep(0.72,length(table_dives$impact_force[which(table_dives$landing_cat2=="dive_unsuccessful" &
                                                               table_dives$hunting_strategy_final =="perching"&
                                                               table_dives$sex=="Female")])),1),
       table_dives$impact_force[which(table_dives$landing_cat2=="dive_unsuccessful" &
                                        table_dives$hunting_strategy_final =="perching"&
                                        table_dives$sex=="Female")],
       pch = 16,
       cex=0.5,col=alpha("#F7C02D",0.5))


## confidence interval arounnd mean
segments(x0 = c(0.28,0.32,0.68,0.72), y0 = c((40.5- (1.96*0.497)),
                                             ((38.3- (1.96*0.401))),
                                             ((38.6- (1.96*0.630))),
                                             ((38.2- (1.96*0.519)))), # dive Female ci around the mean
         
         x1 = c(0.28,0.32,0.68,0.72), y1 = c((40.5+ (1.96*0.497)),
                                             ((38.3+ (1.96*0.401))),
                                             ((38.6+ (1.96*0.630))),
                                             ((38.2+ (1.96*0.519)))), # dive Female ci around the mean
         col = "black",
         lwd = 2)

points(c(0.28,0.32,0.68,0.72),c(40.5,38.3,38.6,38.2),pch = 19,col="black", cex=0.5)


dev.off()


########################
########################
# Figure Paper ####
packages<-function(x){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x,repos="http://cran.r-project.org")
    require(x,character.only=TRUE)
  }
}
packages(tidyverse) 
packages(tidylog) 
packages(lubridate)
packages(move)
packages(maptools)
packages(maps)
packages(recurse)
packages(remotes)
packages(ggmap)
packages(mapproj)
packages(moveVis)
packages(leaflet)
packages(mapview)
packages(devtools)
packages(viridis)
packages(data.table)
#packages(elevatr)
packages(raster)
packages(rgdal)
packages(rnaturalearth)
packages(rnaturalearthdata)  

## Figure1_map####
## figure 1. map of all tracks and zoom in on one sequence of landing and dive

all_indiv_land_complete_final_data<-fread("/Volumes/aroulin1/chouette/D2c/papers/data/all_indiv_land_complete_final_data.csv")

# import all indiv tracks
tag_id<- unique(all_indiv_land_complete_final_data$TagID)
tag_id_2019<-unique(all_indiv_land_complete_final_data$TagID[all_indiv_land_complete_final_data$year=="2019"])
tag_id_2020<-unique(all_indiv_land_complete_final_data$TagID[all_indiv_land_complete_final_data$year=="2020"])

full_data_2019<- fread("/Volumes/aroulin1/chouette/D2c/Annotated_data_complete/annotated_individuals_2019_complete_withland.csv") %>% 
  filter(TagID %in%tag_id_2019)
full_data_2020<- fread("/Volumes/aroulin1/chouette/D2c/Annotated_data_complete/annotated_individuals_2020_complete_withland.csv") %>% 
  filter(TagID %in%tag_id_2020)

names(full_data_2019)


## remove duplicate timestamps
#2019
duplicates_2019 <- move::getDuplicatedTimestamps(x = full_data_2019$TagID,timestamps = full_data_2019$timestamp)
dupl.rows.2019 <- vector("integer")
for(i in 1:length(duplicates_2019)){
  x <- duplicates_2019[[i]][2]
  dupl.rows.2019[i] <- x
  rm(x)
}

full_data_2019_corr <- full_data_2019[-dupl.rows.2019,] %>% 
  rename_all(funs(make.names(.))) %>% 
  group_by(TagID) %>% 
  data..table::setorder(timestamp) %>% 
  as_tibble() 

#2020
duplicates_2020 <- move::getDuplicatedTimestamps(x = full_data_2020$TagID,timestamps = full_data_2020$timestamp)
dupl.rows.2020 <- vector("integer")
for(i in 1:length(duplicates_2020)){
  x <- duplicates_2020[[i]][2]
  dupl.rows.2020[i] <- x
  rm(x)
}
full_data_2020_corr <- full_data_2020[-dupl.rows.2020,] %>% 
  #as_tibble() %>% 
  arrange(timestamp) %>% 
  arrange(TagID) %>% 
  rename_all(funs(make.names(.))) %>% 
  as_tibble()



## convnersion move format
# 2019
full_data_2019.move <- move(x = full_data_2019_corr$location.lon,
                            y =full_data_2019_corr$location.lat,
                            time = full_data_2019_corr$timestamp,
                            proj = CRS("+proj=longlat +ellps=WGS84"),
                            animal = full_data_2019_corr$TagID,
                            data = full_data_2019_corr)

# 2020
full_data_2020.move <- move(x = full_data_2020_corr$location.lon,
                            y =full_data_2020_corr$location.lat,
                            time = full_data_2020_corr$timestamp,
                            proj = CRS("+proj=longlat +ellps=WGS84"),
                            animal = full_data_2020_corr$TagID,
                            data = full_data_2020_corr)
## GPS error correction

setMethod("distance", 
          signature=c(x=".MoveTrackSingle",y="missing"),
          definition=function(x){ 
            Dists<-raster::pointDistance(x[-sum(n.locs(x)),], x[-1,], longlat=isLonLat(x))
            return(Dists)
          })
setMethod("distance", 
          signature=c(".MoveTrackStack",y="missing"), 
          definition=function(x){
            lst <- lapply(split(x), distance)
            return(lst)
          })

# 2019
full_data_2019.move$distance.m <- unlist(lapply(distance(full_data_2019.move), c, NA))
full_data_2019.back_df<-as.data.frame(full_data_2019.move)
full_data_2019.err_corr<-subset(full_data_2019.back_df, distance.m<100)
full_data_2019.err_corr.move<-move(x = full_data_2019.err_corr$coords.x1,
                                   y =full_data_2019.err_corr$coords.x2,
                                   time = full_data_2019.err_corr$timestamps,
                                   proj = CRS("+proj=longlat +ellps=WGS84"),
                                   animal = full_data_2019.err_corr$trackId,
                                   data = full_data_2019.err_corr)

setwd("/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures")
save(full_data_2019.err_corr.move,file = "full_data_2019.err_corr.move.RData")
save(full_data_2019.err_corr,file = "full_data_2019.err_corr.RData")



# 2020
full_data_2020.move$distance.m <- unlist(lapply(distance(full_data_2020.move), c, NA))
full_data_2020.back_df<-as.data.frame(full_data_2020.move)
full_data_2020.err_corr<-subset(full_data_2020.back_df, distance.m<100)
full_data_2020.err_corr.move<-move(x = full_data_2020.err_corr$coords.x1,
                                   y =full_data_2020.err_corr$coords.x2,
                                   time = full_data_2020.err_corr$timestamps,
                                   proj = CRS("+proj=longlat +ellps=WGS84"),
                                   animal = full_data_2020.err_corr$trackId,
                                   data = full_data_2020.err_corr)


setwd("/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures")
save(full_data_2020.err_corr,file = "full_data_2020.err_corr.RData")
save(full_data_2020.err_corr.move,file = "full_data_2020.err_corr.move.RData")



# pour le front end du cluster

setwd("/users/kschalch/work/Figures")

###
load("full_data_2019.err_corr.RData")
load("full_data_2019.err_corr.move.RData")

full_data_2019.err_corr<- test_124
full_data_2019.err_corr.move<- test_126


load("full_data_2020.err_corr.RData" )
load("full_data_2020.err_corr.move.RData")


## mapping
# 2019
m_2019 <- get_stamenmap(bbox(extent(full_data_2019.err_corr.move)*1.05), source="stamen", maptype = "terrain-background", zoom=12)
map.tracks.2019 <- ggmap(m_2019, darken = 0) +
  theme_bw() +
  # geom_path(data = females.data, aes(x=location.long, y=location.lat, group = individual.local.identifier, colour = avg.rev.night),
  #           alpha = 0.5, size = 0.1) +
  geom_path(data = full_data_2019.err_corr, aes(x=location.lon, y=location.lat, group = TagID),
            alpha = 0.5, size = 0.1) +
    theme(legend.position = "right",
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(colour = "black", size = 15, face = "bold")) +
  labs(x = "Longitude",
       y = "Latitude",
       colour = "Revisits")

figpath<-"/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures"
ggsave(paste0(figpath, "Figure_1.map.pdf"), plot = map.tracks, width = 6, height = 6, units = "in", dpi = 300) 

#2020
m_2020 <- get_stamenmap(bbox(extent(full_data_2020.err_corr.move)*1.05), source="stamen", maptype = "terrain-background", zoom=12)
map.tracks.2020 <- ggmap(m_2020, darken = 0) +
  theme_bw() +
  # geom_path(data = females.data, aes(x=location.long, y=location.lat, group = individual.local.identifier, colour = avg.rev.night),
  #           alpha = 0.5, size = 0.1) +
  geom_path(data = full_data_2020.err_corr, aes(x=location.lon, y=location.lat, group = TagID),
            alpha = 0.5, size = 0.1) +
  theme(legend.position = "right",
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(colour = "black", size = 15, face = "bold")) +
  labs(x = "Longitude",
       y = "Latitude",
       colour = "Revisits")
# 2019 + 2020

m_2019 <- get_stamenmap(bbox(extent(full_data_2019.err_corr.move)*1.05), source="stamen", maptype = "terrain-background", zoom=12)
m_2019 <- get_stamenmap(bbox(extent(full_data_2019.err_corr.move)*1.05), source="stamen", maptype = "watercolor", zoom=12)

map.tracks.2019_2020 <- ggmap(m_2019, darken = 0) +
  theme_bw() +
  # geom_path(data = females.data, aes(x=location.long, y=location.lat, group = individual.local.identifier, colour = avg.rev.night),
  #           alpha = 0.5, size = 0.1) +
  geom_path(data = full_data_2019.err_corr, aes(x=location.lon, y=location.lat, group = TagID),
            alpha = 0.5, size = 0.1) +
  geom_path(data = full_data_2020.err_corr, aes(x=location.lon, y=location.lat, group = TagID),
            alpha = 0.5, size = 0.1) +
   theme(legend.position = "right",
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(colour = "black", size = 15, face = "bold")) +
  labs(x = "Longitude",
       y = "Latitude",
       colour = "Revisits")

#figpath<-"/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures"
figpath<-"/users/kschalch/work/Figures/"
ggsave(paste0(figpath, "Figure_1.2.map.pdf"), plot = map.tracks.2019_2020, width = 6, height = 6, units = "in", dpi = 300) 


##############################
##############################
## Figure1_perch_sequence #### 3rd version Final
##############################
##############################

#owl track
track<-fread("/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/female_M038111_170719_behaviour_GPS_meteo.csv")
track_sub<-track %>% 
  filter(night_nb_real==1) %>% 
  dplyr::filter(timestamp>"2019-07-18 02:57:55" & timestamp<"2019-07-18 03:13:10") 
#filter(hourafterdusk > 5.22444166666667 & hourafterdusk<5.29471944444445)
track_sub$hourafterdusk

track_sf<-st_as_sf(track_sub,coords = c('location-lon', 'location-lat'), crs=st_crs(4326))
track_sf<-sf::st_transform(track_sf,crs = st_crs(2056))
track_sf_df<- as.data.frame(track_sf)
coord<-data.frame(st_coordinates(track_sf))
track_sf_df$coords.x1<-coord$X
track_sf_df$coords.x2<-coord$Y


#download raster image, crop and merge thhem together with stack to have al colour layer

img1<-raster::stack("/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/image aériene figure 1/swissimage-dop10_2020_2528-1153_0.1_2056.tif")

LAT1 = 1153000 ; LAT2 = 1153550
LON1 = 2528460 ; LON2 = 2528860

ext<-extent(c(LON1,LON2,LAT1,LAT2))
#plot(ext)
img1_crop<-crop(img1,ext)
#img2_crop<-crop(img2,ext)
fond_de_carte=img1_crop
#fond_de_carte<-merge(img1_crop,img2_crop)

#install.packages("rasterVis")
#library(rasterVis)

df <- as.data.frame(fond_de_carte, xy= TRUE)
head(df)
df <- df %>% rename(Red =  swissimage.dop10_2020_2528.1153_0.1_2056_1,   #Rename bands
                    Green =  swissimage.dop10_2020_2528.1153_0.1_2056_2,
                    Blue =  swissimage.dop10_2020_2528.1153_0.1_2056_3)

# creat the path and percinng poitns


mov.pts<-track_sf_df%>% 
  filter(fly_bmk_category=="fly") #%>% 
#filter(GPS_id!=90318 & GPS_id !=90319)

perch.pts<-track_sf_df%>% 
  filter(fly_bmk_category=="Perching") %>% 
  group_by(fly_bmk_id) %>% 
  summarize(longitudes=median(coords.x1),
            latitudes=median(coords.x2),
            minutes=n()/60,
            acc_max=first(vectorial_sum_acc_max),
            timestamp=timestamp[1]) #%>% 
#filter(!fly_bmk_id %in% c(212))


#perch.pts.all<-bind_rows(perch.pts,strike.pts)
perch.pts.all<-perch.pts


#movement tracks ffor landing force
owl.track<-fread("/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/all_individuals_complete_land.csv")
owl.track_sub<- owl.track[owl.track$TagID=="Female_M038111_170719",] %>% 
  dplyr::filter(night_nb_real=="1") %>% 
  #dplyr::  filter(hourafterdusk >= 5.2244& hourafterdusk<5.294719) %>% 
  dplyr::filter(timestamp>"2019-07-18 02:57:50" & timestamp<"2019-07-18 03:13:05")  %>% 
  dplyr::distinct() %>% 
  as.data.frame()


### anontating to the table
# perch.pts.all$landing_impact<-owl.track_sub$max_vecsum
# perch.pts.all$landing_force<-(owl.track_sub$max_vecsum*9.81)*(290/1000)
perch.pts.all$landing_impact<-owl.track_sub$max_vecsum
perch.pts.all$landing_force<-(owl.track_sub$max_vecsum*9.81)*(290/1000)



####
# final mapping
###
install.packages("ggspatial")
library(raster)
library(RColorBrewer)
library(ggspatial)
black<-raster(extent(df))
black<-raster::setValues(black,1)
df_black <- as.data.frame(black, xy= TRUE)




map_final<- ggplot(data = df, aes(x = x, y =y))+                   #plot map
  geom_raster(fill = rgb(red = df$Red,
                         green = df$Green,
                         blue = df$Blue,
                         maxColorValue = 255),show.legend = FALSE) +
  geom_raster(data=df_black, aes(x = x, y =y), fill="black", color="black", alpha=0.3)+
  labs(title = "",x = "Longitude", y="Latitude")+
  annotation_scale(location = "bl", width_hint = 0.1, text_cex = 0.75, bar_cols = c("black", "white")) +
  theme_minimal()+
  #geom_label(data=bias.dtf, aes(x=LON ,y=LAT, label=ID, fontface=7), hjust=1, vjust=0, size=4) +
  # geom_path(data = track_sub.move.reproj.df, aes(x=coords.x1, y= coords.x2, group = TagID),
  #           alpha = 0.8, size = 1)+
  geom_path(data = mov.pts, aes(x=coords.x1, y= coords.x2,  group = TagID),
            colour="white",alpha = 0.8 ,size =0.8)+
  # geom_path(data = mov.pts, aes(x=coords.x1, y= coords.x2,  group = TagID),
  #           colour="white",alpha = 0.5 ,size =1)+
  # geom_point(data = mov.pts, aes(x=coords.x1, y= coords.x2,  group = TagID),
  #           alpha = 0.7, size = 1, colour="orange")+
  # geom_point(data = perch.pts.all, aes(x=longitudes, y= latitudes, group = fly_bmk_id, size=minutes,  colour=landing_force),
  #            alpha = 1)+
  geom_point(data = perch.pts.all, aes(x=longitudes, y= latitudes, group = fly_bmk_id, size=landing_force,colour=landing_force),
             alpha = 1)+
  
  geom_point(data = perch.pts.all, aes(x=longitudes, y= latitudes, group = fly_bmk_id, size=landing_force),
             shape = 21,colour = "white",alpha = 0.8)+
  # # geom_point(data = perch.pts, aes(x=longitudes, y= latitudes, group = fly_bmk_id,size=landing_force,  colour=landing_force),
  #            alpha = 1)+
  # geom_point(data = strike.pts, aes(x=longitudes, y= latitudes, group = fly_bmk_id,size=landing_force,  colour=landing_force),
  #            alpha = 1)+
  #scale_size(range = c(2, 5))+
  scale_radius(breaks=c(5,10,15,20,25,30), limits = c(5,26),range = c(2, 4))+
  scale_colour_viridis(breaks=c(5,10,15,20,25,30),limits = c(0,26),option = "B", direction = -1)+
  
  guides(color= guide_legend(), size=guide_legend())+
  
  
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.background =element_blank(),legend.key=element_rect(fill="white"))

map_final
?guides
figpath<-"/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/"
#figpath<-"/users/kschalch/work/Figures/"
ggsave(paste0(figpath, "Figure_1.tracks_final.pdf"), plot = map_final, width = 5, height = 6, units = "in", dpi = 300) 


#########
######### figure SI foraginng trip ground truth 

#owl track
track<-fread("/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/female_M038111_170719_behaviour_GPS_meteo.csv")
# track_sub<-track %>% 
#   group_by(night_nb_real,foraging_trip_nb) %>% 
#   summarize(n=n())
# 
#   filter(night_nb_real==1) %>% 

track<-fread("/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/female_M038111_170719_behaviour_GPS_meteo.csv")
track_sub<-track %>% 
  filter(night_nb_real==1) %>% 
  filter(foraging_trip_nb==2)

track_sf<-st_as_sf(track_sub,coords = c('location-lon', 'location-lat'), crs=st_crs(4326))
track_sf<-sf::st_transform(track_sf,crs = st_crs(2056))
track_sf_df<- as.data.frame(track_sf)
coord<-data.frame(st_coordinates(track_sf))
track_sf_df$coords.x1<-coord$X
track_sf_df$coords.x2<-coord$Y
plot(track_sf)

#download raster image, crop and merge thhem together with stack to have al colour layer

img1<-raster::stack("/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/image aériene figure 1/swissimage-dop10_2020_2528-1153_0.1_2056.tif")
img2<-raster::stack("/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/image aériene figure 1/swissimage-dop10_2020_2529-1154_0.1_2056.tif")
img3<-raster::stack("/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/image aériene figure 1/swissimage-dop10_2020_2529-1153_0.1_2056.tif")
# img2<-raster::stack("/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/image aériene figure 1/swissimage-dop10_2020_2561-1189_0.1_2056.tif")
# img3<-raster::stack("/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/image aériene figure 1/swissimage-dop10_2020_2526-1157_0.1_2056.tif")
# img4<-raster::stack("/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/image aériene figure 1/swissimage-dop10_2020_2527-1157_0.1_2056.tif")

LAT1 = 1153130 ; LAT2 = 1154850
LON1 = 2528500 ; LON2 = 2529700

ext<-extent(c(LON1,LON2,LAT1,LAT2))
#plot(ext)
img1_crop<-crop(img1,ext)
img2_crop<-crop(img2,ext)
img3_crop<-crop(img3,ext)
# img4_crop<-crop(img4,ext)
#fond_de_carte=img1_crop
fond_de_carte<-merge(img1_crop,img2_crop)

#install.packages("rasterVis")
#library(rasterVis)

df <- as.data.frame(fond_de_carte, xy= TRUE)
head(df)
summary(df$layer.1)
df <- df %>% rename(Red =  layer.1,   #Rename bands
                    Green =  layer.2,
                    Blue =  layer.3) 

df<- df %>% 
  drop_na(Red)

# creat the path and percinng poitns


mov.pts<-track_sf_df%>% 
  filter(fly_bmk_category=="fly") #%>% 
#filter(GPS_id!=90318 & GPS_id !=90319)

'%ni%' <- Negate("%in%")

strike.pts<-track_sf_df%>% 
  filter(dive_bmk_category=="dive" ) %>% 
  group_by(fly_bmk_id,fly_bmk_category,dive_bmk_id) %>% 
  summarize(longitudes=median(coords.x1, na.rm = T),
            latitudes=median(coords.x2, na.rm = T),
            minutes=n()/60,
            acc_max=max(vectorial_sum_acc_max),
            hourafterdusk=round(hourafterdusk[1],2)) %>%
  filter(fly_bmk_category=="fly")%>% 
  ungroup() %>% 
  dplyr::select(-c(fly_bmk_category,dive_bmk_id)) %>% 
  arrange(fly_bmk_id)


perch.pts<-track_sf_df%>% 
  filter(fly_bmk_category=="Perching") %>% 
  group_by(fly_bmk_id) %>% 
  summarize(longitudes=median(coords.x1, na.rm = T),
            latitudes=median(coords.x2, na.rm = T),
            minutes=n()/60,
            acc_max=max(vectorial_sum_acc_max),
            hourafterdusk=round(hourafterdusk[1],2)) %>% 
  filter(fly_bmk_id %ni% c(126,128,129,132,133,136,135,134,137,138,140,141,142)) 


perch.pts.all<-bind_rows(perch.pts,strike.pts) %>% 
  arrange(fly_bmk_id)



plot(perch.pts.all$latitudes~ perch.pts.all$longitudes,col="red")
points(land_conv$Y~land_conv$X, col="blue")
points(strike.pts$latitudes~strike.pts$longitudes, col="green")
points(perch.pts.all$latitudes~perch.pts.all$longitudes, col="red")
lines(mov.pts$coords.x2~mov.pts$coords.x1)

#movement tracks ffor landing force
owl.track<-fread("/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/all_individuals_complete_land.csv")

owl.track_sub<- owl.track[owl.track$TagID=="Female_M038111_170719",] %>% 
  dplyr::filter(night_nb_roost==1) %>% 
  dplyr::filter(foraging_trip_nb==2)%>% 
  dplyr::distinct() %>% 
  #filter(row_number() <= n()-1) %>% 
  as.data.frame() %>% 
  mutate(hourafterdusk=round(hourafterdusk,2)) %>% 
  group_by(hourafterdusk) %>% 
  filter(row_number()==1) %>% 
  ungroup() %>% 
  filter(hourafterdusk %ni% c(0.93,1.56,1.64,1.83)) 



# land_conv<-st_as_sf(owl.track_sub,coords = c("location.lon","location.lat"), crs=st_crs(4326)) %>% 
#   sf::st_transform(.,crs = st_crs(2056)) %>% 
#   # as.data.frame() %>% 
#   data.frame(st_coordinates(.)) 


plot(owl.track_sub$location.lat~owl.track_sub$location.lon)
plot(land_conv$Y~land_conv$X, col="green")
points(perch.pts.all$latitudes~ perch.pts.all$longitudes)

perch.pts.all$landing_impact<-owl.track_sub$max_vecsum
perch.pts.all$landing_force<-(owl.track_sub$max_vecsum*9.81)*(322/1000)
perch.pts.all$cat<-owl.track_sub$landing_cat2


####
# final mapping
###
library(raster)
library(RColorBrewer)
black<-raster(extent(df))
black<-raster::setValues(black,1)
df_black <- as.data.frame(black, xy= TRUE)

gc()
library(ggplot2)

map_final<- ggplot(data = df, aes(x = x, y =y))+                   #plot map
  
  labs(title = "",x = "Longitude", y="Latitude")+
    geom_path(data = mov.pts, aes(x=coords.x1, y= coords.x2,  group = TagID),
            colour="grey",alpha = 0.8 ,size =0.3)+

  geom_point(data = perch.pts.all, aes(x=longitudes, y= latitudes, group = fly_bmk_id, shape=cat, colour=landing_force),
             alpha = 1, size=1)+
  
  scale_radius(breaks=c(5,10,15,20,25,30,35,40,45,50), limits = c(5,60),range = c(1, 3))+
  scale_colour_viridis(breaks=c(5,10,15,20,25,30,35,40,45,50),limits = c(0,60),option = "B", direction = -1)+
  scale_fill_viridis()+
  
  guides(color= guide_legend(), size=guide_legend())+
  
  
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.background =element_blank(),legend.key=element_rect(fill="white"))

map_final

figpath<-"/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/"
#figpath<-"/users/kschalch/work/Figures/"
ggsave(paste0(figpath, "Fig_path_SI_scale.pdf"), plot = map_final, width = 6, height = 6, units = "in", dpi = 300) 



## figure sex comp ####

male<- fread("/Volumes/aroulin1/chouette/D2c/Annotated_data_acc_GPS_meteo_preys/raw_data_new/2019/Male_M043801_140619_behaviour_GPS_meteo_preys.csv") %>% 
  dplyr::filter(night_nb_real==2) %>% 
  dplyr::filter(foraging_trip_nb %in% c(1,2,3,4,5,6,7,8,9,10,11))
male_sf<-st_as_sf(male,coords = c('location-lon', 'location-lat'), crs=st_crs(4326))
male_sf<-sf::st_transform(male_sf,crs = st_crs(2056))
male_sf_df<- as.data.frame(male_sf)
male_coord<-data.frame(st_coordinates(male_sf))
male_sf_df$coords.x1<-male_coord$X
male_sf_df$coords.x2<-male_coord$Y


female<-fread("/Volumes/aroulin1/chouette/D2c/Annotated_data_acc_GPS_meteo_preys/raw_data_new/2019/Female_M042109_140619_behaviour_GPS_meteo_preys.csv") %>% 
  dplyr::filter(night_nb_real==2) %>% 
  dplyr::filter(foraging_trip_nb %in% c(1,2,3,4))
female_sf<-st_as_sf(female,coords = c('location-lon', 'location-lat'), crs=st_crs(4326))
female_sf<-sf::st_transform(female_sf,crs = st_crs(2056))
female_sf_df<- as.data.frame(female_sf)
female_coord<-data.frame(st_coordinates(female_sf))
female_sf_df$coords.x1<-female_coord$X
female_sf_df$coords.x2<-female_coord$Y

df_tot<- rbind(female_sf_df,male_sf_df)
unique(df_tot$TagID)

## fond de carte “
img1<-raster::stack("/Users/kschalch/Downloads/swissimage-dop10_2020_2524-1157_2_2056.tif")
img2<-raster::stack("/Users/kschalch/Downloads/swissimage-dop10_2020_2524-1158_2_2056.tif")
img3<-raster::stack("/Users/kschalch/Downloads/swissimage-dop10_2020_2525-1157_2_2056.tif")
img4<-raster::stack("/Users/kschalch/Downloads/swissimage-dop10_2020_2525-1158_2_2056.tif")
img5<-raster::stack("/Users/kschalch/Downloads/swissimage-dop10_2020_2524-1156_2_2056.tif")
img6<-raster::stack("/Users/kschalch/Downloads/swissimage-dop10_2020_2525-1156_2_2056.tif")
img7<-raster::stack("/Users/kschalch/Downloads/swissimage-dop10_2020_2526-1157_2_2056.tif")
img8<-raster::stack("/Users/kschalch/Downloads/swissimage-dop10_2020_2526-1156_2_2056.tif")
img9<-raster::stack("/Users/kschalch/Downloads/swissimage-dop10_2020_2526-1158_2_2056.tif")
img10<-raster::stack("/Users/kschalch/Downloads/swissimage-dop10_2020_2525-1155_2_2056.tif")
img11<-raster::stack("/Users/kschalch/Downloads/swissimage-dop10_2020_2524-1155_2_2056.tif")
img12<-raster::stack("/Users/kschalch/Downloads/swissimage-dop10_2020_2526-1155_2_2056.tif")


LAT1 = 1155800 ; LAT2 = 1158500
LON1 = 2524600 ; LON2 = 2526200

ext<-extent(c(LON1,LON2,LAT1,LAT2))
#plot(ext)
img1_crop<-crop(img1,ext)
img2_crop<-crop(img2,ext)
img3_crop<-crop(img3,ext)
img4_crop<-crop(img4,ext)
img5_crop<-crop(img5,ext)
img6_crop<-crop(img6,ext)
img7_crop<-crop(img7,ext)
img8_crop<-crop(img8,ext)
img9_crop<-crop(img9,ext)
img10_crop<-crop(img10,ext)
img11_crop<-crop(img11,ext)
img12_crop<-crop(img12,ext)

extent(img10_crop)
rm(list = c("img1","img2","img3","img4","img5","img6","img7","img8","img9","img10","img11"))


# img4_crop<-crop(img4,ext)
#fond_de_carte=img1_crop
fond_de_carte<-merge(img1_crop,img2_crop,img3_crop,img4_crop,
                     img5_crop,img6_crop,img7_crop,img8_crop,
                     img9_crop,img10_crop,img11_crop,img12_crop)
#fond_de_carte<-merge(img1_crop,img2_crop)
head(fond_de_carte)

#install.packages("rasterVis")
#library(rasterVis)

df <- as.data.frame(fond_de_carte, xy= TRUE)
head(df)
summary(df$layer.1)
df <- df %>% rename(Red =  layer.1,   #Rename bands
                    Green =  layer.2,
                    Blue =  layer.3) 

df<- df %>% 
  drop_na(Red)


head(df)


## plot
# install.packages("RColorBrewer")
# install.packages("ggnewscale")
# install.packages("paleteer")

#library(paleteer)
library(ggnewscale)
library(RColorBrewer)
library(scales)


df_tot<-rbind(female_sf_df,male_sf_df)

dat_male <- df_tot %>% filter(TagID == "Male_M043801_140619")
dat_female <- df_tot %>% filter(TagID == "Female_M042109_140619")

# Create custom color palettes for each TagID
#
palette_male <- c("#AA6200" ,"#B87200", "#C58100", "#D19000", "#DC9D00", "#E5AA00", "#EEB600", "#F5C100", "#FACB00", "#FED300", "#FED800")
palette_female <- c("#023FA5", "#7D87B9", "#BEC1D4", "#E2E2E2")


## 
strike.pts<-df_tot%>% 
  filter(dive_bmk_category=="dive") %>% 
  group_by(TagID,foraging_trip_nb,fly_bmk_id,fly_bmk_category,dive_bmk_id) %>% 
  summarize(longitudes=median(coords.x1, na.rm = T),
            latitudes=median(coords.x2, na.rm = T),
            minutes=n()/60,
            acc_max=max(vectorial_sum_acc_max),
            hourafterdusk=round(hourafterdusk[1],2)) %>%
  filter(fly_bmk_category=="fly")%>% 
  ungroup() %>% 
  dplyr::select(-c(fly_bmk_category,dive_bmk_id)) %>% 
  arrange(fly_bmk_id)


strike.pts.succ<-df_tot%>% 
  filter(dive_bmk_category=="dive" & dive_bmk_sucess==1 ) %>% 
  group_by(TagID,foraging_trip_nb,fly_bmk_id,fly_bmk_category,dive_bmk_id) %>% 
  summarize(longitudes=median(coords.x1, na.rm = T),
            latitudes=median(coords.x2, na.rm = T),
            minutes=n()/60,
            acc_max=max(vectorial_sum_acc_max),
            hourafterdusk=round(hourafterdusk[1],2)) %>%
  filter(fly_bmk_category=="fly")%>% 
  ungroup() %>% 
  dplyr::select(-c(fly_bmk_category,dive_bmk_id)) %>% 
  arrange(fly_bmk_id)

strike.pts.unsucc<-df_tot%>% 
  filter(dive_bmk_category=="dive" & dive_bmk_sucess==0 ) %>% 
  group_by(TagID,foraging_trip_nb,fly_bmk_id,fly_bmk_category,dive_bmk_id) %>% 
  summarize(longitudes=median(coords.x1, na.rm = T),
            latitudes=median(coords.x2, na.rm = T),
            minutes=n()/60,
            acc_max=max(vectorial_sum_acc_max),
            hourafterdusk=round(hourafterdusk[1],2)) %>%
  filter(fly_bmk_category=="fly")%>% 
  ungroup() %>% 
  dplyr::select(-c(fly_bmk_category,dive_bmk_id)) %>% 
  arrange(fly_bmk_id)

strike.pts.unsucc_m<-strike.pts.unsucc %>% 
  dplyr::filter(TagID=="Male_M043801_140619")

strike.pts.unsucc_f<-strike.pts.unsucc %>% 
  dplyr::filter(TagID=="Female_M042109_140619")



strike.pts.succ_m<-strike.pts.succ %>% 
  dplyr::filter(TagID=="Male_M043801_140619")


strike.pts.succ_f<-strike.pts.succ %>% 
  dplyr::filter(TagID=="Female_M042109_140619")


`%!in%` = Negate(`%in%`)

perch.pts<-df_tot%>% 
  dplyr::filter(fly_bmk_category=="Perching") %>% 
  dplyr::filter(fly_bmk_id %!in% unique(strike.pts$fly_bmk_id+1)) %>% 
  group_by(TagID,foraging_trip_nb,fly_bmk_id) %>% 
  summarize(longitudes=median(coords.x1, na.rm = T),
            latitudes=median(coords.x2, na.rm = T),
            minutes=n()/60,
            acc_max=max(vectorial_sum_acc_max),
            hourafterdusk=round(hourafterdusk[1],2))#


perch.pts$foraging_trip_nb
perch.pts_m<-perch.pts%>% 
  dplyr::filter(TagID=="Male_M043801_140619")


perch.pts_f<-perch.pts %>% 
  dplyr::filter(TagID=="Female_M042109_140619")

library(ggmap)

map_sexcomp<-ggplot(data = df, aes(x = x, y =y))+
  geom_raster(fill = rgb(red = df$Red,
                         green = df$Green,
                         blue = df$Blue,
                         maxColorValue = 255),alpha=1,show.legend = FALSE) +
  geom_raster(alpha=0.5)+
  #geom_raster(data=df_black, aes(x = x, y =y), fill="black", color="black", alpha=0.3)+
  labs(title = "",x = "Longitude", y="Latitude")+
  
  # scale_size_continuous(name = "Minutes",breaks=c(0,1,2,4,6,8,10) ,range = c(1, 13)) +  
  new_scale_color() +
  
  geom_path(data = dat_male, aes(x = coords.x1, y = coords.x2, group = interaction(TagID, foraging_trip_nb), color = as.factor(foraging_trip_nb)), size = 1) +
  scale_color_manual(
    name = "Foraging Trip Male",
    values = palette_male,
    breaks = unique(dat_male$foraging_trip_nb),
    labels = unique(dat_male$foraging_trip_nb)
  ) +
  
  new_scale_color() +
  
  geom_path(data = dat_female, aes(x = coords.x1, y = coords.x2, group = interaction(TagID, foraging_trip_nb), color = as.factor(foraging_trip_nb)), size = 1) +
  
  scale_color_manual(
    name = "Foraging Trip Female",
    values = palette_female,
    breaks = unique(dat_female$foraging_trip_nb),
    labels = unique(dat_female$foraging_trip_nb)
  ) +
  new_scale_color() +
  
  geom_point(data = perch.pts_m, aes(x=longitudes, y= latitudes,  color = as.factor(foraging_trip_nb)),
             alpha = 1, shape = 15, size = 4)+
  geom_point(data = strike.pts.unsucc_m, aes(x=longitudes, y= latitudes, color = as.factor(foraging_trip_nb)),
             alpha = 1, size =5)+
  
  scale_color_manual(
    name = "hunting attempts male",
    values = palette_male,
    breaks = unique(dat_male$foraging_trip_nb),
    labels = unique(dat_male$foraging_trip_nb)
  ) +
  new_scale_color() +
  
  geom_point(data = perch.pts_f, aes(x=longitudes, y= latitudes, color = as.factor(foraging_trip_nb)),
             alpha = 1, shape=15, size=4)+
  
  geom_point(data = strike.pts.unsucc_f, aes(x=longitudes, y= latitudes, color = as.factor(foraging_trip_nb)),
             alpha = 1, size =5)+
  scale_color_manual(
    name = "hunting attempts female",
    values = palette_female,
    breaks = unique(dat_female$foraging_trip_nb),
    labels = unique(dat_female$foraging_trip_nb)
  ) +
  new_scale_color()+
  geom_point(data = perch.pts, aes(x = longitudes, y = latitudes,group = interaction(TagID, foraging_trip_nb)), colour="white",size = 4,shape=0,stroke = 0.5) +
  new_scale_color() +
  
  geom_point(data = strike.pts.unsucc, aes(x = longitudes, y = latitudes, group = interaction(TagID, foraging_trip_nb)), colour="white",size=5,shape=1,stroke = 0.5) +
  new_scale_color() +
  geom_point(data = strike.pts.succ, aes(x = longitudes, y = latitudes),color = "#1DCCFF", size = 5, shape = 17) +
  geom_point(data = strike.pts.succ, aes(x = longitudes, y = latitudes),color = "white", size = 5, shape = 2) +
  
  theme_minimal() +
  guides(color = guide_legend(title = "Foraging Trip"))

figpath<-"/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/"
ggsave(paste0(figpath, "Fig_sex_comp_map.pdf"), plot = map_sexcomp, width = 60, height = 30, units = "cm", dpi = 600) 


#

## sub figure sex comp ####
####### model landing force 
all_indiv_land_complete_final_data2<-all_indiv_land_complete_final_data %>% 
  drop_na(WindSpeed)

model_strike_force_divecat<- lmer(log(impact_force) ~  
                                    ###landing_cat*
                                    
                                    ###sex+
                                    #WindSpeed+
                                    #landing_cat*sex+
                                    
                                    landing_cat/(1+sex)-1+
                                    #landing_cat:WindSpeed+
                                    ##landing_cat:WindSpeed:sex+
                                    (1|TagID/night_nb_real),
                                  data = all_indiv_land_complete_final_data, REML = F)

performance(model_strike_force_divecat)
simulationOutput <- simulateResiduals(fittedModel = model_strike_force_divecat, plot = T)
plot_model(model_strike_force_divecat, type = "diag")

summary(model_strike_force_divecat)
library(emmeans)
stat_landingforce<-emmeans(model_strike_force_divecat, pairwise~ landing_cat*sex)
tbl_stat_landingforce<-data.frame(stat_landingforce$emmeans)

####### mass specific forces
boxplot(mass_spec_force_tbl$bodymass~mass_spec_force_tbl$sex)

mass_spec_force_tbl<-all_indiv_land_complete_final_data %>% 
  mutate(bodymass_kg=bodymass/1000,
         mass_spec_landing_force=impact_force/bodymass_kg) %>% 
  as_tibble()

hist(log(mass_spec_force_tbl$mass_spec_landing_force))

model_mass_spec_force<-lmer(log(mass_spec_landing_force)~
                              sex*landing_cat+
                              (1|TagID/night_nb_real),
                            data=mass_spec_force_tbl)

summary(model_mass_spec_force)
tab_model(model_mass_spec_force,  transform="exp")

plot_model(model_mass_spec_force, type="pred", terms = c("landing_cat","sex"), show.data = F, jitter = 0.1)

stat_forceunitmass<-emmeans(model_mass_spec_force, pairwise~ sex*landing_cat , nesting = NULL)
tbl_stat_forceunitmass=data.frame(stat_forceunitmass$emmeans)

all_indiv_land_complete_final_data$landing_cat




###########
###########
palette_male <- c("#AA6200" ,"#B87200", "#C58100", "#D19000", "#DC9D00", "#E5AA00", "#EEB600", "#F5C100", "#FACB00", "#FED300", "#FED800")
palette_female <- c("#023FA5", "#7D87B9", "#BEC1D4", "#E2E2E2")

col_male<-palette_male[8]
col_female<-palette_female[2]

all_indiv_land_complete_final_data<-fread("/Volumes/aroulin1/chouette/D2c/papers/data/all_indiv_land_complete_final_data.csv") 
names(all_indiv_land_complete_final_data)
all_indiv_land_complete_final_data$impact_force_mass_spec<-all_indiv_land_complete_final_data$impact_force/(all_indiv_land_complete_final_data$bodymass/1000)


#summary force per ind


landing_force_summary <- all_indiv_land_complete_final_data %>%
  group_by(TagID,sex,night_nb_real,landing_cat) %>%
  summarise(
    sd = sd(impact_force, na.rm = TRUE),
    impact_force = mean(impact_force,na.rm=T)
  )
landing_force_summary
mybreaks=log(c(5,10,20,40,80,160))
mylabels=exp(mybreaks)

plot_forcevscat <- ggplot(data = all_indiv_land_complete_final_data, aes(x = landing_cat, y = log(impact_force), fill = sex)) +
  
  geom_point(aes(color = sex), position = position_jitterdodge(dodge.width = 0.95, jitter.width = 0.3, jitter.height = 0), size = 2, alpha = 0.2) +
  
  geom_violin(color= "black",position = position_dodge(0.95),alpha=0.5) +
  geom_boxplot(color= "black",outlier.alpha = 0, position = position_dodge(0.95),alpha=0.5) +
  
  geom_pointrange(data = tbl_stat_landingforce,
                  aes(x = landing_cat, y = emmean, ymin = asymp.LCL, ymax = asymp.UCL, group = sex),
                  position = position_dodge(width = 0.95), color = "white", size=0.1) +
  
  scale_y_continuous(breaks = mybreaks, 
                     labels = mylabels)+
  # geom_point(aes(color = sex), position = position_jitterdodge(dodge.width = 0.95, jitter.width = 0.45, jitter.height = 0), size = 2, alpha = 0.5) +
  scale_colour_manual(
    name = "Sex",
    values = alpha(c(col_female, col_male), 0.5),
    labels = c( "Females","Males")
  ) +
  scale_fill_manual(
    name = "Sex",
    values = alpha(c(col_female, col_male), 0.5),
    labels = c("Females","Males")
  ) +
  facet_wrap(~landing_cat) +
  
  
  
  #facet_grid(~landing_cat) +
  theme_minimal()
?scale_y_continuous

figpath<-"/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/"
ggsave(paste0(figpath, "Fig_sex_comp_landingforce.pdf"), plot = plot_forcevscat, width = 40, height = 10, units = "cm", dpi = 600) 



#### plot vs per mass unit



mybreaks=log(c(5,10,20,40,80,160,320,640))
mylabels=exp(mybreaks)

plot_force_massvscat <- ggplot(data = all_indiv_land_complete_final_data, aes(x = landing_cat, y = log(impact_force_mass_spec), fill = sex)) +
  
  geom_point(aes(color = sex), position = position_jitterdodge(dodge.width = 0.95, jitter.width = 0.3, jitter.height = 0), size = 2, alpha = 0.2) +
  
  geom_violin(color= "black",position = position_dodge(0.95),alpha=0.5) +
  geom_boxplot(color= "black",outlier.alpha = 0, position = position_dodge(0.95),alpha=0.5) +
  
  geom_pointrange(data = tbl_stat_forceunitmass,
                  aes(x = landing_cat, y = emmean, ymin = asymp.LCL, ymax = asymp.UCL, group = sex),
                  position = position_dodge(width = 0.95), color = "white", size=0.1) +
  
  scale_y_continuous(breaks = mybreaks, 
                     labels = mylabels)+
  scale_colour_manual(
    name = "Sex",
    values = alpha(c(col_female, col_male), 0.5),
    labels = c( "Females","Males")
  ) +
  scale_fill_manual(
    name = "Sex",
    values = alpha(c(col_female, col_male), 0.5),
    labels = c("Females","Males")
  ) +
  facet_wrap(~landing_cat) +
  theme_minimal()

plot_force_massvscat

figpath<-"/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/"
ggsave(paste0(figpath, "Fig_sex_comp_landingforcemass.pdf"), plot = plot_force_massvscat, width = 40, height = 10, units = "cm", dpi = 600) 




#body mass * sex
dat_bodymass<-all_indiv_land_complete_final_data %>% 
  group_by(TagID) %>% 
  summarise(bodymass=bodymass[1],
            sex=sex[1])


dat_bodymass_summary <- dat_bodymass %>%
  group_by(sex) %>%
  summarise(
    sd = sd(bodymass, na.rm = TRUE),
    bodymass = mean(bodymass),
    nb=n())
dat_bodymass_summary


plot_bodymass<-ggplot(data=dat_bodymass, aes(x= sex, y = bodymass, fill=sex))+
  geom_point(aes(color = sex), position = position_jitterdodge(dodge.width = 0.95, jitter.width = 0.35, jitter.height = 0), size = 4, alpha = 0.5) +
  
  geom_violin( color= "black",position = position_dodge(0.95),alpha=0.5)+
  geom_boxplot( color= "black",outlier.alpha = 0, position = position_dodge(0.95),alpha=0.5, coef = 0)+
  
  geom_pointrange(
    aes(y=bodymass, ymin = bodymass-sd, ymax = bodymass+sd),
    data = dat_bodymass_summary, color="white"
  )+
  
  #facet_wrap(~landing_cat, scales = "free")+
  scale_fill_manual(
    name = "Sex",
    values = alpha(c(col_female,col_male),0.7),
    labels=c("Females","Males"))+
  geom_jitter(color="black", width=0.45,size=1, alpha=0.01) +
  scale_color_manual(
    name = "Sex",
    values = alpha(c(col_female,col_male),0.7),
    labels=c("Females","Males"))+
  geom_jitter(color="black", width=0.45,size=1, alpha=0.01) +
  theme_minimal() 
plot_bodymass


figpath<-"/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/"
ggsave(paste0(figpath, "Fig_sex_comp_bodymass.pdf"), plot = plot_bodymass, width = 40, height = 40, units = "cm", dpi = 600) 



table(all_indiv_land_complete_final_data$sex, all_indiv_land_complete_final_data$landing_cat)

#foraging activity (nb hunting attampts )


foraging_act_dat<-all_indiv_land_complete_final_data %>% 
  dplyr::filter(foraging_trip_nb!=88) %>% 
  group_by(TagID,sex, night_nb_real, landing_cat) %>% 
  summarize(nb=n()) %>% 
  dplyr::filter(nb<300)

dat_foraging_act_dat_summary <- foraging_act_dat %>%
  group_by(sex,landing_cat) %>%
  summarise(
    sd = sd(nb, na.rm = TRUE),
    mean = mean(nb,na.rm=T),
    nb = sum(nb)
  )
dat_foraging_act_dat_summary


plot_foraging_act<-ggplot(foraging_act_dat, aes(y=mean, x=landing_cat, fill=sex)) +
  geom_point(aes(color = sex), position = position_jitterdodge(dodge.width = 0.95, jitter.width = 0.35, jitter.height = 0), size = 4, alpha = 0.5) +
  geom_violin( color= "black",position = position_dodge(0.95),alpha=0.5)+
  geom_boxplot( color= "black",outlier.alpha = 0, position = position_dodge(0.95),alpha=0.5, coef = 0)+
  
  geom_pointrange(
    aes(y=nb, ymin = mean-sd, ymax = mean+sd),
    data = dat_foraging_act_dat_summary,
    position = position_dodge(width = 0.95), color = "white", size=0.4
  )+
  
  scale_fill_manual(
    name = "Sex",
    values = alpha(c(col_female,col_male),0.7),
    labels=c("Females","Males"))+
  scale_color_manual(
    name = "Sex",
    values = alpha(c(col_female,col_male),0.7),
    labels=c("Females","Males"))+
  theme_minimal()

plot_foraging_act

figpath<-"/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/"
ggsave(paste0(figpath, "Fig_sex_comp_activity.pdf"), plot = plot_foraging_act, width = 40, height = 40, units = "cm", dpi = 600) 


# foraging trip nb* sex

dat_food_supply<-all_indiv_land_complete_final_data %>% 
  group_by(TagID,sex, night_nb_real) %>% 
  mutate(last_ft_nb=last(foraging_trip_nb)) %>% 
  dplyr::filter(landing_cat=="dive") %>% 
  dplyr::filter(foraging_trip_nb!=88) %>% 
  group_by(TagID,sex, night_nb_real) %>% 
  slice_tail() %>% 
  dplyr::select(foraging_trip_nb, last_ft_nb) %>% 
  mutate(food_supply=NA) %>% 
  ungroup() %>% 
  mutate(food_supply=ifelse(last_ft_nb==88,foraging_trip_nb,foraging_trip_nb-1),
         food_supply_cor=case_when(food_supply==-1 ~ 0,
                                   food_supply==0 ~ 1,
                                   T ~ food_supply))

table(dat_food_supply$sex)

dat_food_supply_summary <- dat_food_supply %>%
  group_by(sex) %>%
  summarise(
    sd = sd(food_supply_cor, na.rm = TRUE),
    food_supply = mean(food_supply_cor,na.rm=T),
    nb = sum(food_supply_cor)
  )
dat_food_supply_summary


plot_food_supp<-ggplot(dat_food_supply, aes(y=food_supply_cor, x = sex, fill=sex))+
  geom_point(aes(color = sex), position = position_jitterdodge(dodge.width = 0.95, jitter.width = 0.35, jitter.height = 0), size = 4, alpha = 0.5) +
  geom_violin( color= "black",position = position_dodge(0.95),alpha=0.5)+
  geom_boxplot( color= "black",outlier.alpha = 0, position = position_dodge(0.95),alpha=0.5, coef = 0)+
  
  geom_pointrange(
    aes(y=food_supply, ymin = food_supply-sd, ymax = food_supply+sd),
    data = dat_food_supply_summary,
    position = position_dodge(width = 0.95), color = "white", size=0.4
  )+
  
  scale_fill_manual(
    name = "Sex",
    values = alpha(c(col_female,col_male),0.7),
    labels=c("Females", "Males"))+
  scale_color_manual(
    name = "Sex",
    values = alpha(c(col_female,col_male),0.7),
    labels=c("Females", "Males"))+
  theme_minimal()
plot_food_supp

figpath<-"/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/"
ggsave(paste0(figpath, "Fig_sex_comp_foodsupply.pdf"), plot = plot_food_supp, width = 40, height = 40, units = "cm", dpi = 600) 

# hunting strategy prop * sex

dat_hunting_strat_prop<-all_indiv_land_complete_final_data %>%
  dplyr::filter(landing_cat=="dive",
                hunting_strategy_final2!="") %>% 
  group_by(TagID, sex, night_nb_real, hunting_strategy_final2) %>%
  summarise(attempt_count = n()) %>%
  complete(hunting_strategy_final2, fill = list(attempt_count = 0)) 



dat_hunting_strat_prop <-all_indiv_land_complete_final_data %>%
  dplyr::filter(landing_cat=="dive",
                hunting_strategy_final2!="") %>% 
  dplyr::select(-sex) %>% 
  group_by(TagID, night_nb_real, hunting_strategy_final2) %>%
  summarise(attempt_count = n()) %>%
  ungroup() %>% 
  complete(TagID, night_nb_real,hunting_strategy_final2, fill = list(attempt_count = 0)) %>%
  group_by(TagID, night_nb_real) %>%
  mutate(total_attempts = sum(attempt_count)) %>%
  #ungroup() %>%
  mutate(strategy_proportion = attempt_count / total_attempts) %>%
  dplyr:: select(-total_attempts) %>% 
  dplyr::filter(!is.na(strategy_proportion)) %>% 
  mutate(sex=str_split_fixed(TagID, "_", 3)[,1]) %>% 
  group_by(TagID,sex,hunting_strategy_final2) %>% 
  summarize(strategy_proportion=mean(strategy_proportion, na.rm=T),
            nb = n()) %>% 
  dplyr::filter(TagID!="Male_M037707_260619") %>% 
  dplyr::filter(hunting_strategy_final2=="perching")


dat_hunting_strat_prop_summary <- dat_hunting_strat_prop %>%
  group_by(sex, hunting_strategy_final2) %>%
  summarise(
    sd = sd(strategy_proportion, na.rm = TRUE),
    prop = mean(strategy_proportion,na.rm=T),
    nb = n()
  )
dat_hunting_strat_prop_summary



plot_hunt_prop<-ggplot(dat_hunting_strat_prop, aes(y=strategy_proportion, x = hunting_strategy_final2, fill=sex))+
  geom_point(aes(color = sex), position = position_jitterdodge(dodge.width = 0.95, jitter.width = 0.35, jitter.height = 0), size = 4, alpha = 0.5) +
  geom_violin( color= "black",position = position_dodge(0.95),alpha=0.5)+
  geom_boxplot( color= "black",outlier.alpha = 0, position = position_dodge(0.95),alpha=0.5, coef = 0)+
  
  geom_pointrange(
    aes(y=prop, ymin = prop-sd, ymax = prop+sd),
    data = dat_hunting_strat_prop_summary,
    position = position_dodge(width = 0.95), color = "white", size=0.4
  )+
  scale_fill_manual(
    name = "Sex",
    values = alpha(c(col_female,col_male),0.7),
    labels=c("Females", "Males"))+
  scale_color_manual(
    name = "Sex",
    values = alpha(c(col_female,col_male),0.7),
    labels=c("Females", "Males")) +
  theme_minimal()

plot_hunt_prop

figpath<-"/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/"
ggsave(paste0(figpath, "Fig_sex_comp_hunt_strat_prop.pdf"), plot = plot_hunt_prop, width = 40, height = 40, units = "cm", dpi = 600) 


## test relation hunt strategy ~ time
names(all_indiv_land_complete_final_data)
all_indiv_land_complete_final_data$hunting_strategy_0s
data_hunt_time<-all_indiv_land_complete_final_data %>% 
  dplyr::select(TagID, night_nb_real,foraging_trip_nb,landing_cat, hunting_strategy_0s,timestamp) %>% 
  dplyr::filter(foraging_trip_nb!=0 & foraging_trip_nb!=88)

data_hunt_time_sumerise<-data_hunt_time %>% 
  group_by(TagID, night_nb_real, foraging_trip_nb) %>% 
  mutate(time=difftime(last(timestamp),first(timestamp), units="min")) %>% 
  dplyr::filter(landing_cat=="dive") %>% 
  group_by(TagID, night_nb_real,foraging_trip_nb,hunting_strategy_0s) %>% 
  summarise(time=time[1],
            nb=n()) %>% 
  group_by(TagID) %>% 
  complete(night_nb_real,foraging_trip_nb,hunting_strategy_0s) %>% 
  mutate_at("nb",~replace_na(., 0)) %>% 
  group_by(TagID, night_nb_real,foraging_trip_nb) %>% 
  tidyr::fill(time, .direction="downup") %>% 
  group_by(TagID, night_nb_real,foraging_trip_nb) %>% 
  mutate(nb_tot=sum(nb)) %>% 
  ungroup() %>% 
  mutate(freq=nb/nb_tot) %>% 
  dplyr::filter(hunting_strategy_0s == "perching") %>% 
  mutate(time=as.numeric(time)) %>% 
  dplyr::filter(time!=0) %>% 
  mutate(sex=str_split_fixed(TagID, "_", 3)[,1]) 



hist(data_hunt_time_sumerise$time, breaks =100)
summary(data_hunt_time_sumerise$time)



model_time<- lmer(log(time)~
                    nb_tot+
                    freq +
                    # I(freq^2)+
                    sex+
                    #freq:sex+
                    + (1|TagID/night_nb_real) , data=data_hunt_time_sumerise)

summary(model_time)
check_model(model_time)
plot_model(model_time, type = "pred", terms=c("freq[all]"), show.data = F)+
  geom_point(data=data_hunt_time_sumerise, aes(x=freq, y=time))

plot_time_strat<-plot_model(model_time, type = "pred", terms=("sex[all]"), show.data = F)

figpath<-"/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/"
ggsave(paste0(figpath, "Figure_strat_vs_time.pdf"), plot =  plot_time_strat, width = 12, height = 6, units = "in", dpi = 300) 


tab_model(model_time, transform = "exp")

tab_model(model_time,string.est = "Estimate", 
          auto.label = T,
          string.ci = "Conf. Int.", string.p = "p-value",transform = "exp" , file = "/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/trip_duration_vs_strategy.html" )

library(webshot)
webshot("/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/trip_duration_vs_strategy.html", "/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/trip_duration_vs_strategy.pdf")



taskid_moon<- data.frame(x=seq(1,31951, by=10), y=c(seq(10,31950, by =10),31953), z=seq(1,3196))



?data.frame

?seq







# flight speed*sex
mod_flightspeed<-lmer(speed_before_strike_med ~ Sex.y +
                        (1|TagID/night_nb_real),data = tb_mod_flightspeed)
flight_speed_stat=emmeans(mod_flightspeed, pairwise~ Sex.y, 
                          type="response")

table(tb_mod_flightspeed$Sex.y)

dat.flight_speed_stat<-data.frame(flight_speed_stat$emmeans)

plot_flight_speed <- ggplot(data = tb_mod_flightspeed, aes(x = Sex.y, y = speed_before_strike_mean, fill = Sex.y)) +
  
  geom_point(aes(color = Sex.y), position = position_jitterdodge(dodge.width = 0.95, jitter.width = 0.45, jitter.height = 0), size = 2, alpha = 0.2) +
  
  geom_violin(color= "black",position = position_dodge(0.95),alpha=0.5) +
  geom_boxplot(color= "black",outlier.alpha = 0, position = position_dodge(0.95),alpha=0.5, coef=0) +
  
  geom_pointrange(data = dat.flight_speed_stat,
                  aes(x = Sex.y, y = emmean, ymin = asymp.LCL, ymax = asymp.UCL, group = Sex.y),
                  position = position_dodge(width = 0.95), color = "white", size=0.1) +
  
  # geom_point(aes(color = sex), position = position_jitterdodge(dodge.width = 0.95, jitter.width = 0.45, jitter.height = 0), size = 2, alpha = 0.5) +
  scale_colour_manual(
    name = "Sex",
    values = alpha(c(col_female, col_male), 0.5),
    labels = c( "Females","Males")
  ) +
  scale_fill_manual(
    name = "Sex",
    values = alpha(c(col_female, col_male), 0.5),
    labels = c("Females","Males")
  ) +
  facet_wrap(~hunting_strategy_final) +
  theme_minimal()

plot_flight_speed

figpath<-"/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/"
ggsave(paste0(figpath, "Fig_sex_comp_hunt_flight_speed.pdf"), plot = plot_flight_speed, width = 40, height = 40, units = "cm", dpi = 600) 






#hunting success ~ foraging strategy * sex
dive_tbl_full$hunting_strategy_final


succ_sex_strat_dat<-dive_tbl_full %>% 
  group_by(TagID,sex, hunting_strategy_final, dive_bmk_sucess) %>% 
  summarize(nb=n()) %>% 
  ungroup() %>% 
  complete(TagID, hunting_strategy_final, dive_bmk_sucess, fill = list(nb = 0)) %>% 
  mutate(sex=str_split_fixed(TagID, "_", 3)[,1]) %>% 
  group_by(TagID,hunting_strategy_final) %>% 
  mutate(tot = sum(nb)) %>% 
  mutate(succ_prop = nb / tot) %>% 
  dplyr::filter(dive_bmk_sucess==1)


succ_sex_strat_datt_summary <- succ_sex_strat_dat %>%
  group_by(sex,hunting_strategy_final) %>%
  summarise(
    sd = sd(succ_prop, na.rm = TRUE),
    prop = mean(succ_prop,na.rm=T),
    nb = sum(tot)
  )
succ_sex_strat_datt_summary


try<-glmer(dive_bmk_sucess ~ 
             sex*hunting_strategy_final+
             #hunting_strategy_final/(1+sex) -1 +
             (1|TagID/night_nb_real),
           family =binomial(link = "logit"), data = dive_tbl_full)

summary(try)

emm.test<-emmeans(try, pairwise~ hunting_strategy_final*sex, 
                  type="response")

dat.emm.test<-data.frame(emm.test$emmeans)



plot_succ_sex_strat <- ggplot(data = succ_sex_strat_dat, aes(x = hunting_strategy_final, y = succ_prop, fill = sex)) +
  
  geom_point(aes(color = sex), position = position_jitterdodge(dodge.width = 0.95, jitter.width = 0.45, jitter.height = 0), size = 2, alpha = 0.2) +
  
  geom_violin(color= "black",position = position_dodge(0.95),alpha=0.5) +
  geom_boxplot(color= "black",outlier.alpha = 0, position = position_dodge(0.95),alpha=0.5, coef=0) +
  
  geom_pointrange(data = dat.emm.test,
                  aes(x = hunting_strategy_final, y = prob, ymin = asymp.LCL, ymax = asymp.UCL, group = sex),
                  position = position_dodge(width = 0.95), color = "white", size=0.1) +
  
  # geom_point(aes(color = sex), position = position_jitterdodge(dodge.width = 0.95, jitter.width = 0.45, jitter.height = 0), size = 2, alpha = 0.5) +
  scale_colour_manual(
    name = "Sex",
    values = alpha(c(col_female, col_male), 0.5),
    labels = c( "Females","Males")
  ) +
  scale_fill_manual(
    name = "Sex",
    values = alpha(c(col_female, col_male), 0.5),
    labels = c("Females","Males")
  ) +
  facet_wrap(~hunting_strategy_final) +
  theme_minimal()

plot_succ_sex_strat

figpath<-"/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/"
ggsave(paste0(figpath, "Fig_sex_comp_hunt_hunt_strat_succ.pdf"), plot = plot_succ_sex_strat, width = 30, height = 10, units = "cm", dpi = 600) 


#######
# acceleration sequence
######

acc_seq<-fread("/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/female_M038111_170719_behaviour_GPS_meteo.csv") %>% 
  filter(night_nb_real==1) %>% 
  dplyr::filter(timestamp>"2019-07-18 02:57:55" & timestamp<"2019-07-18 03:13:10") 
plot(acc_seq$vectorial_sum_acc_max, type="l")
min(acc_seq$timestamp)
max(acc_seq$timestamp)

# take the full sequences
accc_sequ_full<-fread("/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/Female_M038111_170719.csv")

#set time datetime

# set timestamp

Timestamp <- as.POSIXct(fast_strptime(accc_sequ_full$Timestamp, format="%d/%m/%Y %H:%M:%OS", tz="GMT")) ### plus vite avec fast strp time!!!
accc_sequ_full$timestamp<-Timestamp

# filter to the corresdponding time
accc_sequ_full_sub<-accc_sequ_full %>% 
  dplyr::filter(timestamp>"2019-07-18 02:57:00" & timestamp<"2019-07-18 03:13:50") 


# acc correction
# nozero<-which(accc_sequ_full_sub$X!=0)
# accc_sequ_full_sub<-accc_sequ_full_sub[-(1:nozero[1]),,drop=F]

# invert x and y and multiply acc columns by 0.62
accc_sequ_full_sub<- accc_sequ_full_sub %>%
  mutate(X=accc_sequ_full_sub$Y,Y=accc_sequ_full_sub$X) %>%
  mutate(X=(X/0.62), Y=(Y/0.62), Z=(Z/0.62))
# calculation of acceleration metrix

freq<- 50 #The Frequency of accelerometry data
secs<- 0.5 # the number of seconds over which to calculate the desired metrics.The manuscript says to use 1 second intervals,

numrows<-freq*secs # the number of rows required to calculate metrics over the chosen period.

#calculation of smoothed acceleration
library(zoo)
??rollapply
Accx_smoothed=rollapply(accc_sequ_full_sub$X,numrows,mean,fill=NA)#25 = number of cells to average across, change to desired time considering sampling frequency.
Accy_smoothed=rollapply(accc_sequ_full_sub$Y,numrows,mean,fill=NA)#25 = number of cells to average across, change to desired time considering sampling frequency.
Accz_smoothed=rollapply(accc_sequ_full_sub$Z,numrows,mean,fill=NA)#25 = number of cells to average across, change to desired time considering sampling frequency.

## calculation of ODBA and VEDBA

## first calculates the static body acceleration by applying a running mean --> correspond to the acc_smoothed

static_Accx= Accx_smoothed
static_Accy= Accy_smoothed
static_Accz= Accz_smoothed



#### Calculates DBA for each axis.

dynamic_Accx<-accc_sequ_full_sub$X-static_Accx
dynamic_Accy<-accc_sequ_full_sub$Y-static_Accy
dynamic_Accz<-accc_sequ_full_sub$Z-static_Accz
ODBA<-abs(dynamic_Accx)+abs(dynamic_Accy)+abs(dynamic_Accz)
VeDBA<-sqrt((dynamic_Accx^2)+(dynamic_Accy^2)+(dynamic_Accz^2))

### calculation of Vedba smoothed
### watch out here for VedBA we smooth over 2 secs
library(slider)
secs_VeDBA<- 2
numrows_VeDBA<-freq*secs_VeDBA
VeDBA_smoothed = slide_vec(VeDBA, mean, .before = 50,.after = 50)

### calculation vectorial sum of the raw acceleration
vec_sum_raw_acc= sqrt((accc_sequ_full_sub$X)^2 + (accc_sequ_full_sub$Y)^2 + (accc_sequ_full_sub$Z)^2)

#### ajout des metrics acceleration a la table acc and attribution des flight strategy

accc_sequ_full_sub<- accc_sequ_full_sub %>%
  as_tibble() %>% 
  mutate(static_Accx =static_Accx,
         static_Accy =static_Accy,
         static_Accz =static_Accz,
         dynamic_Accx =  dynamic_Accx,
         dynamic_Accy =  dynamic_Accy,
         dynamic_Accz =  dynamic_Accz,
         ODBA = ODBA,
         VeDBA = VeDBA,
         VeDBA_smoothed = VeDBA_smoothed,
         vectorial_sum_acc = vec_sum_raw_acc,
         force= (vec_sum_raw_acc*9.81)*0.290,
         id = row_number())

##
accc_sequ_full_sub<-accc_sequ_full_sub %>% 
  dplyr::filter(timestamp>"2019-07-18 02:57:55" & timestamp<"2019-07-18 03:13:10")
##plot
max(accc_sequ_full_sub$Z)
plot(accc_sequ_full_sub$Z, type="l", col="blue")
lines(accc_sequ_full_sub$X, type="l", col="red")
lines(accc_sequ_full_sub$Y, type="l", col="green")
lines(accc_sequ_full_sub$force, col="yellow")
lines(accc_sequ_full_sub$VeDBA, col="black")
accc_sequ_full_sub$Timestamp

names(accc_sequ_full_sub)
plot(accc_sequ_full_sub$force, type="l", col="blue")

formatter1000 <- function(x){ 
  t0=round(x[1]/50,1)
  a=(x/50)
  round(a,1)-t0
}

#
perch.pts.all




## sequence zoomed in land
#land heave

accc_sequ_full_sub_zoom_land<-accc_sequ_full_sub %>% 
  dplyr::filter(timestamp>"2019-07-18 03:09:13" & timestamp<"2019-07-18 03:09:35") 

plot_sequence_acc_zoom_land_Heave<-ggplot(accc_sequ_full_sub_zoom_land)+
  geom_line(aes(x=id, y=Z,group = 1),color="black", size=0.2)+
  #geom_line(aes(x=id,y=X,group = 1), color="red", size=0.2)+
  #geom_line(aes(x=id,y=Y,group = 1), color="green",size=0.2)+
  #geom_line(aes(x=id,y=vectorial_sum_acc,group = 1), color="purple",size=0.5)+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(limits = c(36650, 37550)
                     ,breaks = seq(36650, 37550, by = 150),labels = formatter1000)+
  scale_y_continuous(limits = c(-4,10))+
  labs(x="Time (sec)", y="Heave acceleration (g)")


plot_sequence_acc_zoom_land_Heave
figpath<-"/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/"
ggsave(paste0(figpath, "Figure_1.acc_seq_zoom_land_Heave.pdf"), plot =plot_sequence_acc_zoom_land_Heave, width = 6, height = 6, units = "in", dpi = 360) 


# land force

df_land<-accc_sequ_full_sub_zoom_land %>% 
  as_tibble() %>% 
  dplyr::select(id,force) 
land_app<-approx(df_land$id,df_land$force, n=10000)
force_approx_land<-data.frame(id=land_app[1], force=land_app[2])
names(force_approx_land)<-c("id","force")


plot_sequence_acc_zoom_land_force<-ggplot(force_approx_land)+
  #geom_line(aes(x=id, y=Z,group = 1),color="blue", size=0.2)+
  #geom_line(aes(x=id,y=X,group = 1), color="red", size=0.2)+
  #geom_line(aes(x=id,y=Y,group = 1), color="green",size=0.2)+
  geom_line(aes(x=id,y=force,group = 1, colour=force), size=0.5)+
  scale_x_continuous(limits = c(36650, 37550)
                     ,breaks = seq(36650, 37550, by = 150),labels = formatter1000)+
  scale_y_continuous(limits = c(-4,30))+
  scale_colour_viridis(limits = c(0,26),option = "B", direction = -1)+
  labs(x="Time (sec)", y="Force (N)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")


plot_sequence_acc_zoom_land_force
figpath<-"/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/"
ggsave(paste0(figpath, "Figure_1.acc_seq_zoom_land_force.pdf"), plot =plot_sequence_acc_zoom_land_force, width = 6, height = 6, units = "in", dpi = 360) 



## sequence zoomed in strike
# strike heave

acc_sequ_full_sub_zoom_strike<-accc_sequ_full_sub %>% 
  dplyr::filter(timestamp>"2019-07-18 03:12:50" & timestamp<"2019-07-18 03:13:10")

plot_sequence_acc_zoom_strike_Heave<-ggplot(acc_sequ_full_sub_zoom_strike)+
  geom_line(aes(x=id, y=Z,group = 1),color="black", size=0.2)+
  #geom_line(aes(x=id,y=X,group = 1), color="red", size=0.2)+
  #geom_line(aes(x=id,y=Y,group = 1), color="green",size=0.2)+
  #geom_line(aes(x=id,y=vectorial_sum_acc,group = 1), color="purple",size=0.5)+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(limits = c(47900, 48400)
                     ,breaks = seq(47900, 48400, by = 150),labels = formatter1000
  )+
  scale_y_continuous(limits = c(-4,10))+
  labs(x="Time (min)", y="Heave acceleration (g)")

plot_sequence_acc_zoom_strike_Heave
figpath<-"/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/"
ggsave(paste0(figpath, "Figure_1.acc_seq_zoom_strike_Heave.pdf"), plot =plot_sequence_acc_zoom_strike_Heave, width = 6, height = 6, units = "in", dpi = 360) 


# strike force

df_strike<-acc_sequ_full_sub_zoom_strike %>% 
  as_tibble() %>% 
  dplyr::select(id,force) 
strike_app<-approx(df_strike$id,df_strike$force, n=10000)
force_approx_strike<-data.frame(id=strike_app[1], force=strike_app[2])
names(force_approx_strike)<-c("id","force")

plot_sequence_acc_zoom_strike_force<-ggplot(force_approx)+
  #geom_line(aes(x=id, y=Z,group = 1),color="black", size=0.2)+
  #geom_line(aes(x=id,y=X,group = 1), color="red", size=0.2)+
  #geom_line(aes(x=id,y=Y,group = 1), color="green",size=0.2)+
  geom_line(aes(x=id,y=force,group = 1,colour=force),size=0.5)+
  scale_x_continuous(limits = c(47900, 48400)
                     ,breaks = seq(47900, 48400, by = 150),labels = formatter1000
  )+
  scale_y_continuous(limits = c(-4,30))+
  scale_colour_viridis_c(limits = c(0,26),option = "B", direction = -1)+
  labs(x="Time (min)", y="Force (N)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")

plot_sequence_acc_zoom_strike_force
figpath<-"/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/"
ggsave(paste0(figpath, "Figure_1.acc_seq_zoom_strike_force.pdf"), plot =plot_sequence_acc_zoom_strike_force, width = 6, height = 6, units = "in", dpi = 360) 


#require(gridExtra)
min(acc_sequ_full_sub_zoom_strike$force)

grid.arrange(plot_sequence_acc_zoom_land_Heave,
             plot_sequence_acc_zoom_strike_Heave,
             plot_sequence_acc_zoom_land_force,
             plot_sequence_acc_zoom_strike_force, ncol=2)




test<-acc_sequ_full_sub_zoom_strike %>% 
  as_tibble() %>% 
  dplyr::select(id,force) 

app<-approx(test$id,test$force, n=10000)

force_approx<-df<-data.frame(id=app[1], force=app[2])
names(force_approx)<-c("id","force")


plot_sequence_acc_zoom_strike_force_test<-ggplot(force_approx)+
 
  geom_line(aes(x=id,y=force,group = 1,colour=force),size=0.5)+
  scale_x_continuous(limits = c(47900, 48400)
                     ,breaks = seq(47900, 48400, by = 150),labels = formatter1000
  )+
  scale_y_continuous(limits = c(-4,30))+
  scale_colour_viridis_c(limits = c(0,26),option = "B", direction = -1)+
  labs(x="Time (min)", y="Force (N)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")





## tag deployment folder 

Tag_deployment<-fread("/Users/kschalch/Documents/Tag_Deployment_230421.csv") 


all_indiv_land_complete_final_data<-fread("/Volumes/aroulin1/chouette/D2c/papers/data/all_indiv_land_complete_final_data.csv") 
indiv_per_year<-all_indiv_land_complete_final_data %>% 
  group_by(TagID) %>% 
  slice(1)

depl_ind<- indiv_per_year %>% 
  left_join(Tag_deployment,by = c("TagID" = "Data_Directory")) %>% 
  dplyr::select(Date_attached,Weight_attached,Date_recapture,Weight_recapture,GPS_TotalWeight) %>% 
  mutate(Date_attached_new = as.POSIXct(Date_attached,format="%d.%m.%y"),
         Date_recapture_new = as.POSIXct(Date_recapture,format="%d.%m.%y"),
         diff_attach_recapt= difftime(Date_recapture_new,Date_attached_new, units="days"),
         date_record_end = Date_recapture_new + diff_attach_recapt)


deployment_length_mean<- mean(as.numeric(depl_ind$diff_attach_recapt[which(depl_ind$diff_attach_recapt<30)]), na.rm=T)
deployment_length_sd<-sd(as.numeric(depl_ind$diff_attach_recapt[which(depl_ind$diff_attach_recapt<30)]), na.rm=T)
hist(as.numeric(depl_ind$diff_attach_recapt[-11]), breaks=160)
hist(as.numeric(depl_ind$diff_attach_recapt), breaks=160)


  mutate(Year= year(Date_attached_new),
         Month= month(Date_attached_new),
         sex= strsplit(TagID, "[_]")[[1]][1]) %>% 
  group_by(Year,sex) %>% 
  summarize(n=n())


min_month<-depl_ind %>% 
  mutate(Year= year(Date_attached_new),
         Month= month(Date_attached_new)) %>% 
  group_by(Year) %>% 
  summarise_at(vars(Month),
               list(max = min))


max_month<-depl_ind %>% 
  mutate(Year= year(Date_attached_new),
         Month= month(Date_attached_new)) %>% 
  group_by(Year) %>% 
  summarise_at(vars(Month),
               list(max = max))


### nb of tracked pairs, nb of nest boxes! ####
names(Tag_deployment)

indiv_per_year<-all_indiv_land_complete_final_data %>% 
  group_by(TagID) %>% 
  slice(1)

depl_ind2<- indiv_per_year %>% 
  left_join(Tag_deployment,by = c("TagID" = "Data_Directory")) %>% 
  dplyr::select(Date_attached,Weight_attached,Date_recapture,Weight_recapture, Site_name, Installation_Year, SiteID, Sex ) %>% 
  mutate(Date_attached_new = as.POSIXct(Date_attached,format="%d.%m.%y"),
         Date_recapture_new = as.POSIXct(Date_recapture,format="%d.%m.%y"),
         diff_attach_recapt= difftime(Date_recapture_new,Date_attached_new, units="days"),
         date_record_end = Date_recapture_new + diff_attach_recapt) %>% 
  arrange(Installation_Year,Site_name) %>% 
  group_by(Installation_Year, Site_name, Sex) %>% 
  summarise(n=n()) %>% 
  group_by(Installation_Year, n) %>% 
  summarise(nb=n())



### tot nb perch vs hunt per ind per night ####
names(all_indiv_land_complete_final_data)
all_indiv_land_complete_final_data$landing_cat
all_indiv_land_complete_final_data %>% group_by(TagID,night_nb_roost,landing_cat) %>% 
  filter(night_nb_roost!=0) %>% 
  summarise( nb_attempts = n()) %>% 
  #group_by(TagID, landing_cat) %>% 
  group_by(TagID) %>% 
  filter(night_nb_roost!=max(night_nb_roost)) %>% 
  group_by(TagID, landing_cat) %>% 
  #group_by(sex,night_nb_roost,landing_cat) %>% 
  summarise(mean_attempts=mean(nb_attempts), sd_attempts=sd(nb_attempts)) %>% 
  rowwise() %>% 
  mutate(sex= strsplit(TagID, "[_]")[[1]][1]) %>% 
  group_by(sex,landing_cat) %>% 
  summarise(mean_tot_attempts=mean(mean_attempts), sd_tot_attempts=mean(sd_attempts, na.rm=T))


### nb of foraging trips per ind and per night
all_indiv_land_complete_final_data$foraging_trip_nb

all_indiv_land_complete_final_data %>% group_by(TagID,night_nb_roost,foraging_trip_nb) %>% 
  slice(1) %>% 
  filter(night_nb_roost!=0) %>% 
  summarise(ft_nb=n()) %>% 
  group_by(TagID) %>% 
  filter(night_nb_roost!=max(night_nb_roost)) %>% 
  group_by(TagID, night_nb_roost) %>% 
  summarise(ft_nb=n()) %>% 
  group_by(TagID) %>% 
  summarise(mean_ft_nb=mean(ft_nb),
            sd_ft_nb=sd(ft_nb)) %>% 
  rowwise() %>% 
  mutate(sex= strsplit(TagID, "[_]")[[1]][1]) %>% 
  group_by(sex) %>% 
  summarise(mean_ft=mean(mean_ft_nb),
            sd_ft=mean(sd_ft_nb))





##nb perching events per pole types ####

lands2<-fread("/Volumes/aroulin1/chouette/D2c/Papers/data/table_gam_lands_complete_final.csv") %>% 
  dplyr::mutate(perch_type=as.factor(perch_type),
                TagID=as.factor(TagID),
                Tag_night = as.factor(paste0(TagID,night_nb_real))) #%>% 
drop_na(night_nb_real)


perch_prop<-lands2 %>% 
  group_by(TagID, perch_type) %>% 
  summarize(n=n()) %>% 
  group_by(TagID) %>% 
  mutate(freq=n/sum(n)) %>% 
  group_by(perch_type) %>%
  summarise(mean_freq=mean(freq),
            sd_freq=sd(freq))


## reviewer 2 ####

### limit of 5% bodymass ####

depl_ind<- indiv_per_year %>% 
  left_join(Tag_deployment,by = c("TagID" = "Data_Directory")) %>% 
  dplyr::select(Date_attached,Weight_attached,Date_recapture,Weight_recapture,GPS_TotalWeight) %>% 
  mutate(Date_attached_new = as.POSIXct(Date_attached,format="%d.%m.%y"),
         Date_recapture_new = as.POSIXct(Date_recapture,format="%d.%m.%y"),
         diff_attach_recapt= difftime(Date_recapture_new,Date_attached_new, units="days"),
         date_record_end = Date_recapture_new + diff_attach_recapt) %>% 
  mutate(percent_bodymass = GPS_TotalWeight/Weight_attached)


mean(depl_ind$GPS_TotalWeight, na.rm=T)
sd(depl_ind$GPS_TotalWeight, na.rm=T)

mean(depl_ind$percent_bodymass, na.rm=T)
min(depl_ind$percent_bodymass, na.rm=T)
max(depl_ind$percent_bodymass, na.rm=T)

### R2 problem ####
mod_fit_dredge<-glm_test<-glmer(dive_bmk_sucess ~ 
                                  #scale(landing_force_previous_landing)*
                                  ##scale(I(landing_force_previous_landing^2))+
                                  #hunting_strategy_final*
                                  sex+
                                  #scale(time_to_hunt_from_previous_landing)+
                                  #scale(dist_from_previous_landing)+ # check savoir si je l'ai lassé a la fin ou pas
                                  ##scale(I(dist_from_previous_landing^2))+
                                  #scale(WindSpeed)+
                                  ##scale(impact_force)+
                                  hunting_strategy_final/(1+scale(landing_force_previous_landing))-1+
                                  #hunting_strategy_5s_0_1/(1+scale(landing_force_previous_landing))-1+ # this gives same result but better effect
                                  #perch_type_previous_landing+
                                  ##scale(dist_from_previous_landing):scale(landing_force_previous_landing)+
                                  (1|TagID/night_nb_real),
                                family =binomial(link = "logit"), data = dive_tbl_no_na)
plot_model(mod_fit_dredge, terms=c("landing_force_previous_landing","hunting_strategy_final"), type="eff")

summary(mod_fit_dredge)

library(MuMIn)
options(na.action = "na.fail")
options(na.action = "na.omit")
dredge<- dredge(mod_fit_dredge,extra = "margR^2")






# ANALYSES ####

## force vs landing cat ####
all_indiv_land_complete_final_data<-fread("/Volumes/aroulin1/chouette/D2c/papers/data/all_indiv_land_complete_final_data.csv") 
model_strike_force_divecat<- lmer(log(impact_force) ~  
                                    ###landing_cat*
                                    
                                    ###sex+
                                    #WindSpeed+
                                    #landing_cat*sex+
                                    
                                    landing_cat/(1+sex)-1+
                                    #landing_cat:WindSpeed+
                                    ##landing_cat:WindSpeed:sex+
                                    (1|TagID/night_nb_real),
                                  data = all_indiv_land_complete_final_data, REML = F)

## strike force model ####
table_dives<-all_indiv_land_complete_final_data  %>% 
  ungroup() %>% 
  mutate(hunting_strategy_20s= case_when(hunt_strat_flyprop_20sec>0.1 ~ "flying",
                                         T ~ "perching"),
         hunting_strategy_0s= case_when(hunt_strat_flyprop_0sec==1 ~ "flying",
                                        T ~ "perching")) %>% 
  mutate(hunting_strategy_20s= as.factor(hunting_strategy_20s),
         hunting_strategy_0s= as.factor(hunting_strategy_0s)) %>% 
  mutate(hunting_strategy_final=case_when(hunt_strat_flyprop_5sec==1~"flying",
                                          hunt_strat_flyprop_0sec==0~"perching",
                                          T ~ NA_character_)) %>% 
  mutate(hunting_strategy_5s_0_1=case_when(hunt_strat_flyprop_5sec==1~"flying",
                                           hunt_strat_flyprop_5sec==0~"perching",
                                           T ~NA_character_)) %>% 
  filter(landing_cat=="dive") %>% 
  arrange(TagID,timestamp) %>% 
  group_by(TagID) 

model_strike_force_divesucc<- lmer(log(impact_force) ~ landing_cat2+
                                     sex+
                                     hunting_strategy_final+
                                     
                                     #hunting_strategy_5s_0_1+
                                     #WindSpeed+
                                     # landing_cat2:sex + 
                                     #sex:hunting_strategy_final + 
                                     #sex:landing_cat2+
                                     
                                     landing_cat2:hunting_strategy_final+
                                     #landing_cat2:sex+
                                     #hunting_strategy_final:sex+
                                     #landing_cat2:hunting_strategy_5s_0_1+
                                     (1|TagID/night_nb_real),
                                   data = table_dives, REML = F)



# sex comp huntin success ####

setwd("/Volumes/aroulin1/chouette/D2c/Papers/land_study_robj")
load("/Volumes/aroulin1/chouette/D2c/Papers/land_study_robj/data_mod_gam.RData" )
tbl_timetohunt<-lands2 %>% 
  dplyr::select(TagID,landing_cat,landing_id,time_to_hunt_corrected,dist_to_hunt,perch_type)

test<- tbl_timetohunt %>% 
  as_tibble() %>% 
  right_join(all_indiv_land_complete_final_data) %>% 
  group_by(TagID, landing_cat, landing_id) %>% 
  arrange(TagID,landing_id)

library(geosphere)

landing_force_previous_landing<-rep(NA, length(test$TagID)) 
time_to_hunt_from_previous_landing<-rep(NA, length(test$TagID))
dist_from_previous_landing<- rep(NA,length(test$TagID))
perch_type_previous_landing<-rep(NA,length(test$TagID))

i=1

for(i in 2:nrow(test)){
  
  if(test$TagID[i]!=test$TagID[i-1]){
    next
  } else if(test$landing_cat[i]=="landing"){
    next
  }else if( test$landing_cat[i]=="dive"){
    
    if (test$landing_cat[i-1]=="dive"){
      next
    }
    
    landing_force_previous_landing[i]= test$impact_force[i-1]
    time_to_hunt_from_previous_landing[i]=test$time_to_hunt_corrected[i-1]
    dist_from_previous_landing[i]=  abs(distm(c(test$location.lon[i],test$location.lat[i]),
                                              c(test$location.lon[i-1],test$location.lat[i-1])))
    perch_type_previous_landing[i] = as.character(test$perch_type[i-1])
  }
}
test$landing_force_previous_landing<-landing_force_previous_landing
test$time_to_hunt_from_previous_landing<-time_to_hunt_from_previous_landing
test$dist_from_previous_landing<-dist_from_previous_landing
test$perch_type_previous_landing<-perch_type_previous_landing

dive_tbl<-test %>% 
  as_tibble() %>% 
  filter(landing_cat2!="landing") %>% 
  #filter(max_vecsum<9) %>% ### warning check here
  #filter(max_vecsum>1.5) %>% ### warning chheck here
  #filter(hunting_strategy_20s=="perching") %>% 
  #filter(landing_force_previous_landing<30) %>% ## tripple fixed to 40
  mutate(dive_bmk_sucess= as.factor(ceiling(as.numeric(as.character(dive_bmk_sucess))))) %>% 
  mutate(sex=as.factor(sex)) %>% 
  mutate(TagID=as.factor(TagID)) %>% 
  #filter(time_to_hunt_from_previous_landing>0) %>% 
  mutate(noise_pollution2=ifelse(is.na(noise_pollution),0,noise_pollution)) %>% 
  #mutate(hunting_strategy_20s=as.factor(hunting_strategy_20s)) %>% 
  mutate(hunting_strategy_20s_05=case_when(hunt_strat_flyprop_20sec>=0.5~"flying",
                                           T ~"perching")) %>% 
  mutate(hunting_strategy_20s_0_1=case_when(hunt_strat_flyprop_20sec==1~"flying",
                                            hunt_strat_flyprop_20sec==0~"perching",
                                            T ~NA_character_)) %>% 
  mutate(hunting_strategy_5s_0_1=case_when(hunt_strat_flyprop_5sec==1~"flying",
                                           hunt_strat_flyprop_5sec==0~"perching",
                                           T ~NA_character_)) %>% 
  mutate(hunting_strategy_5s_05=case_when(hunt_strat_flyprop_5sec==1~"flying",
                                          hunt_strat_flyprop_0sec==1.~"perching",
                                          T ~ NA_character_)) %>% 
  mutate(hunting_strategy_final=case_when(hunt_strat_flyprop_5sec==1~"flying",
                                          hunt_strat_flyprop_0sec==0~"perching",
                                          T ~ NA_character_)) %>% 
  mutate(hunting_strategy_final2=case_when(hunt_strat_flyprop_5sec==1~"flying",
                                           hunt_strat_flyprop_5sec==0~"perching",
                                           T ~ "hybrid")) %>% 
  mutate(hunting_strategy_20s_05= as.factor(hunting_strategy_20s_05)) %>% 
  mutate(hunting_strategy_20s_0_1= as.factor(hunting_strategy_20s_0_1)) %>% 
  mutate(hunting_strategy_final= as.character(hunting_strategy_final)) %>% 
  #filter(hunting_strategy_final=="perching") %>% 
  filter(time_to_hunt_from_previous_landing<90) %>%  # %>% when tripple interaction, fixed to 500 other wise select only time to hunt within 60 sec
  mutate(Genus= ifelse(Genus == "", NA, Genus)) %>% 
  drop_na(WindSpeed)


#filter(sex=="Male")


#### dive_tbl_all
dive_tbl_full<-test %>% 
  as_tibble() %>% 
  filter(landing_cat2!="landing") %>% 
  #filter(max_vecsum<9) %>% ### warning check here
  #filter(max_vecsum>1.5) %>% ### warning chheck here
  #filter(hunting_strategy_20s=="perching") %>% 
  #filter(landing_force_previous_landing<30) %>% ## tripple fixed to 40
  mutate(dive_bmk_sucess= as.factor(ceiling(as.numeric(as.character(dive_bmk_sucess))))) %>% 
  mutate(sex=as.factor(sex)) %>% 
  mutate(TagID=as.factor(TagID)) %>% 
  #filter(time_to_hunt_from_previous_landing>0) %>% 
  mutate(noise_pollution2=ifelse(is.na(noise_pollution),0,noise_pollution)) %>% 
  #mutate(hunting_strategy_20s=as.factor(hunting_strategy_20s)) %>% 
  mutate(hunting_strategy_20s_05=case_when(hunt_strat_flyprop_20sec>=0.5~"flying",
                                           T ~"perching")) %>% 
  mutate(hunting_strategy_20s_0_1=case_when(hunt_strat_flyprop_20sec==1~"flying",
                                            hunt_strat_flyprop_20sec==0~"perching",
                                            T ~NA_character_)) %>% 
  mutate(hunting_strategy_5s_0_1=case_when(hunt_strat_flyprop_5sec==1~"flying",
                                           hunt_strat_flyprop_5sec==0~"perching",
                                           T ~NA_character_)) %>% 
  mutate(hunting_strategy_5s_05=case_when(hunt_strat_flyprop_5sec==1~"flying",
                                          hunt_strat_flyprop_0sec==1.~"perching",
                                          T ~ NA_character_)) %>% 
  mutate(hunting_strategy_final=case_when(hunt_strat_flyprop_5sec==1~"flying",
                                          hunt_strat_flyprop_0sec==0~"perching",
                                          T ~ NA_character_)) %>% 
  mutate(hunting_strategy_final2=case_when(hunt_strat_flyprop_5sec==1~"flying",
                                           hunt_strat_flyprop_5sec==0~"perching",
                                           T ~ "hybrid")) %>% 
  mutate(hunting_strategy_20s_05= as.factor(hunting_strategy_20s_05)) %>% 
  mutate(hunting_strategy_20s_0_1= as.factor(hunting_strategy_20s_0_1)) %>% 
  mutate(hunting_strategy_final= as.character(hunting_strategy_final)) %>% 
  #filter(hunting_strategy_final=="perching") %>% 
  #filter(time_to_hunt_from_previous_landing<90) %>%  # %>% when tripple interaction, fixed to 500 other wise select only time to hunt within 60 sec
  mutate(Genus= ifelse(Genus == "", NA, Genus))



try<-glmer(dive_bmk_sucess ~ 
             #sex*hunting_strategy_final+
             sex/(1+hunting_strategy_final) -1 +
             (1|TagID/night_nb_real),
           family =binomial(link = "logit"), data = dive_tbl_full)
check_model(try)
summary(try)

tab_model(try,transform = "exp")
tab_model(try,string.est = "Estimate", 
          auto.label = T,
          string.ci = "Conf. Int.", string.p = "p-value",transform = "exp" , file = "/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/table_overall_hunting_success.html" )

library(webshot)
webshot("/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/table_overall_hunting_success.html", "/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/Land_study/Figures/table_overall_hunting_success.html.pdf")

??webshot

dredge(try)


library(sjPlot)
plot_model(glm_test, terms = c("sex", "hunting_strategy_final"), type="pred")
plot_model(glm_test, terms = c("landing_force_previous_landing[all]", "hunting_strategy_final"), type="eff")

plot_model(try, terms = c("sex", "hunting_strategy_final"), type="pred")

library(emmeans)
emm.test<-emmeans(try, pairwise~ hunting_strategy_final*sex, 
                  type="response")

dat.emm.test<-emm.test$emmeans

# landing force vs hunting success







# GAM
gamtest_allsex0<- gam(log(impact_force) ~ 
                        perch_type+
                        WindSpeed+
                        sex+
                        #perch_type/(1+WindSpeed)-1+
                        #sex+
                        s(time_to_hunt_corrected,by =perch_type,k=9)+
                        #s(Tag_night,bs="re") ,
                        s(TagID,bs="re"),## ajouter time to hunt dans "re" pour avoir un random slope
                      data = lands2,
                      #method = "REML",
                      family=gaussian(link="identity"))

summary(gamtest_allsex0)
plot(gamtest_allsex0)
table(lands2$perch_type)

tab_model(gamtest_allsex0, transform = "exp")

gam.check(gamtest_allsex0)


# fliht speed model ####
flight_final<-readRDS("/Volumes/aroulin1/chouette/D2c/script_flight_energetics/final_tables/flight_final_new.rds")# %>%

flight_final<-flight_final %>%
  mutate(pressure_smooth_3sec=slide_vec(Pressure, mean, .before = 1,.after = 1),
         pressure_smooth_5sec=slide_vec(Pressure, mean, .before = 2,.after = 2),
         pressure_smooth_7sec=slide_vec(Pressure, mean, .before = 3,.after = 3))

flight_speed<-flight_final %>%   
  ungroup() %>% 
  dplyr::select(TagID,timestamp,GPS_id,vectorial_sum_acc_max,
                Pressure,location.lat,location.lon,hunt_strat_flyprop_0sec,
                hunt_strat_flyprop_5sec,hourafterdusk,dive_bmk_sucess,dive_bmk_id,
                dive_bmk_category,fly_bmk_id,fly_bmk_category,has_prey,foraging_trip_nb,night_nb_real,
                location.lat_to,location.lon_to,distance_covered,trip_section2,WindSpeed) %>% 
  mutate(time_next=c(timestamp[-1],NA),
         diff_time_s=abs(as.numeric(difftime(timestamp,time_next))),
         flight_speed_ms=(distance_covered/diff_time_s)) %>% 
  group_by(TagID,night_nb_real,foraging_trip_nb,trip_section2,fly_bmk_id,fly_bmk_category)# %>% 
# dplyr::filter(fly_bmk_category=="fly",
#               !is.na(trip_section2),
#               WindSpeed==0)


flight_speed$TagID[flight_speed$TagID=="Male_M031195_250619"]="Male_M041195_250619"



dive_tb_tot<-fread("/Users/kschalch/Desktop/UNIL/Projets/GPS:AXY/moon_project/attack_angle/dive_tbl_angle.csv")
setwd("//Volumes/aroulin1/chouette/D2c/tables_Kerrian")

dive_tb_tot_colour<-readRDS("//Volumes/aroulin1/chouette/D2c/tables_Kerrian/dive_tb_tot_colour.RDS")


dive_tb_tot_speed<-dive_tb_tot_colour %>% 
  mutate(speed_before_strike_mean=as.numeric(NA),
         speed_before_strike_mean_minus2 = as.numeric(NA),
         speed_before_strike_med=as.numeric(NA),
         speed_before_strike_med_minus2 = as.numeric(NA))

#ind ="Female_M038289_140520"

for(ind in unique(dive_tb_tot_speed$TagID)) {
  print(ind)
  
  # owl_data<- owl_data_tot %>% 
  #   dplyr::filter(TagID==ind)
  
  owl_data_speed <- flight_speed[which(flight_speed$TagID==ind),]
  
  # dive_tb<-dive_tb_tot %>% 
  #   dplyr::filter(TagID==ind)
  
  dive_tb_speed <- dive_tb_tot_speed[which(dive_tb_tot_speed$TagID==ind),]
  
  
  # ici on prends la moyenne des flight speed sur les 20 secs avabt chaqze dive 
  for(i in 1:nrow(dive_tb_speed)){  
    
    if(is.na(dive_tb_speed$hunting_strategy_final[i]) | dive_tb_speed$hunting_strategy_final[i]=="perching"){
      
      dive_tb_speed$speed_before_strike_mean[i]=NA
      dive_tb_speed$speed_before_strike_mean_minus2[i]=NA 
      dive_tb_speed$speed_before_strike_med[i]=NA
      dive_tb_speed$speed_before_strike_med_minus2[i]=NA
      
    }
    # ici, on s'assure que les 20 derniers points sont sur le même fly id
    # si fly cat du dive == perch
    
    else if(dive_tb_speed$fly_bmk_category[i]=="fly"){ 
      
      if(length(owl_data_speed$flight_speed_ms[owl_data_speed$fly_bmk_id==dive_tb_speed$fly_bmk_id[i] & owl_data_speed$fly_bmk_category == "fly"])> 20) {
        
        dive_tb_speed$speed_before_strike_mean[i]= mean(owl_data_speed$flight_speed_ms[which(owl_data_speed$GPS_id %in% 
                                                                                               c((dive_tb_speed$GPS_id[i]-20) : dive_tb_speed$GPS_id[i]) & 
                                                                                               owl_data_speed$fly_bmk_category == "fly" &
                                                                                               owl_data_speed$fly_bmk_id== dive_tb_speed$fly_bmk_id[i])], na.rm=T)
        
        dive_tb_speed$speed_before_strike_mean_minus2[i]= mean(owl_data_speed$flight_speed_ms[which(owl_data_speed$GPS_id %in% 
                                                                                                      c((dive_tb_speed$GPS_id[i]-20) : (dive_tb_speed$GPS_id[i]-2)) & 
                                                                                                      owl_data_speed$fly_bmk_category == "fly" &
                                                                                                      owl_data_speed$fly_bmk_id== dive_tb_speed$fly_bmk_id[i])], na.rm=T)
        
        
        dive_tb_speed$speed_before_strike_med[i]= median(owl_data_speed$flight_speed_ms[which(owl_data_speed$GPS_id %in% 
                                                                                                c((dive_tb_speed$GPS_id[i]-20) : dive_tb_speed$GPS_id[i]) & 
                                                                                                owl_data_speed$fly_bmk_category == "fly" &
                                                                                                owl_data_speed$fly_bmk_id== dive_tb_speed$fly_bmk_id[i])], na.rm=T)
        
        dive_tb_speed$speed_before_strike_med_minus2[i]= median(owl_data_speed$flight_speed_ms[which(owl_data_speed$GPS_id %in% 
                                                                                                       c((dive_tb_speed$GPS_id[i]-20) : (dive_tb_speed$GPS_id[i]-2)) & 
                                                                                                       owl_data_speed$fly_bmk_category == "fly" &
                                                                                                       owl_data_speed$fly_bmk_id== dive_tb_speed$fly_bmk_id[i])], na.rm=T)
        
        
        
        
        
      }else {
        
        dive_tb_speed$speed_before_strike_mean[i] = mean(owl_data_speed$flight_speed_ms[which(owl_data_speed$fly_bmk_id==dive_tb_speed$fly_bmk_id[i] & 
                                                                                                owl_data_speed$fly_bmk_category == "fly")], na.rm=T)
        
        dive_tb_speed$speed_before_strike_mean_minus2[i] = mean(owl_data_speed$flight_speed_ms[which(owl_data_speed$fly_bmk_id==dive_tb_speed$fly_bmk_id[i] & 
                                                                                                       owl_data_speed$fly_bmk_category == "fly")][-c((length(which(owl_data_speed$fly_bmk_id==dive_tb_speed$fly_bmk_id[i] & owl_data_speed$fly_bmk_category == "fly"))-1):
                                                                                                                                                       length(which(owl_data_speed$fly_bmk_id==dive_tb_speed$fly_bmk_id[i] & owl_data_speed$fly_bmk_category == "fly")))], na.rm=T)
        
        dive_tb_speed$speed_before_strike_med[i] = median(owl_data_speed$flight_speed_ms[which(owl_data_speed$fly_bmk_id==dive_tb_speed$fly_bmk_id[i] & 
                                                                                                 owl_data_speed$fly_bmk_category == "fly")], na.rm=T)
        
        dive_tb_speed$speed_before_strike_med_minus2[i] = median(owl_data_speed$flight_speed_ms[which(owl_data_speed$fly_bmk_id==dive_tb_speed$fly_bmk_id[i] & 
                                                                                                        owl_data_speed$fly_bmk_category == "fly")][-c((length(which(owl_data_speed$fly_bmk_id==dive_tb_speed$fly_bmk_id[i] & owl_data_speed$fly_bmk_category == "fly"))-1):
                                                                                                                                                        length(which(owl_data_speed$fly_bmk_id==dive_tb_speed$fly_bmk_id[i] & owl_data_speed$fly_bmk_category == "fly")))], na.rm=T)
        
        
        
      }
    } else {
      
      if(length(owl_data_speed$flight_speed_ms[owl_data_speed$fly_bmk_id==dive_tb_speed$fly_bmk_id[i]-1 & owl_data_speed$fly_bmk_category == "fly"])> 20) {
        
        
        
        dive_tb_speed$speed_before_strike_mean[i]= mean(owl_data_speed$flight_speed_ms[which(owl_data_speed$GPS_id %in% 
                                                                                               c((dive_tb_speed$GPS_id[i]-20) : dive_tb_speed$GPS_id[i]) & 
                                                                                               owl_data_speed$fly_bmk_category == "fly" &
                                                                                               owl_data_speed$fly_bmk_id== dive_tb_speed$fly_bmk_id[i]-1)], na.rm=T)
        
        dive_tb_speed$speed_before_strike_mean_minus2[i]= mean(owl_data_speed$flight_speed_ms[which(owl_data_speed$GPS_id %in% 
                                                                                                      c((dive_tb_speed$GPS_id[i]-20) : (dive_tb_speed$GPS_id[i]-2)) & 
                                                                                                      owl_data_speed$fly_bmk_category == "fly" &
                                                                                                      owl_data_speed$fly_bmk_id== dive_tb_speed$fly_bmk_id[i]-1)], na.rm=T)
        
        
        dive_tb_speed$speed_before_strike_med[i]= median(owl_data_speed$flight_speed_ms[which(owl_data_speed$GPS_id %in% 
                                                                                                c((dive_tb_speed$GPS_id[i]-20) : dive_tb_speed$GPS_id[i]) & 
                                                                                                owl_data_speed$fly_bmk_category == "fly" &
                                                                                                owl_data_speed$fly_bmk_id== dive_tb_speed$fly_bmk_id[i]-1)], na.rm=T)
        
        dive_tb_speed$speed_before_strike_med_minus2[i]= median(owl_data_speed$flight_speed_ms[which(owl_data_speed$GPS_id %in% 
                                                                                                       c((dive_tb_speed$GPS_id[i]-20) : (dive_tb_speed$GPS_id[i]-2)) & 
                                                                                                       owl_data_speed$fly_bmk_category == "fly" &
                                                                                                       owl_data_speed$fly_bmk_id== dive_tb_speed$fly_bmk_id[i]-1)], na.rm=T)
        
        
      }else {
        
        
        dive_tb_speed$speed_before_strike_mean[i] = mean(owl_data_speed$flight_speed_ms[which(owl_data_speed$fly_bmk_id==dive_tb_speed$fly_bmk_id[i] -1 & 
                                                                                                owl_data_speed$fly_bmk_category == "fly")], na.rm=T)
        
        dive_tb_speed$speed_before_strike_mean_minus2[i] = mean(owl_data_speed$flight_speed_ms[which(owl_data_speed$fly_bmk_id==dive_tb_speed$fly_bmk_id[i] -1 & 
                                                                                                       owl_data_speed$fly_bmk_category == "fly")][-c((length(which(owl_data_speed$fly_bmk_id==dive_tb_speed$fly_bmk_id[i] & owl_data_speed$fly_bmk_category == "fly"))-1):
                                                                                                                                                       length(which(owl_data_speed$fly_bmk_id==dive_tb_speed$fly_bmk_id[i] & owl_data_speed$fly_bmk_category == "fly")))], na.rm=T)
        
        dive_tb_speed$speed_before_strike_med[i] = median(owl_data_speed$flight_speed_ms[which(owl_data_speed$fly_bmk_id==dive_tb_speed$fly_bmk_id[i]-1 & 
                                                                                                 owl_data_speed$fly_bmk_category == "fly")], na.rm=T)
        
        dive_tb_speed$speed_before_strike_med_minus2[i] = median(owl_data_speed$flight_speed_ms[which(owl_data_speed$fly_bmk_id==dive_tb_speed$fly_bmk_id[i]-1 & 
                                                                                                        owl_data_speed$fly_bmk_category == "fly")][-c((length(which(owl_data_speed$fly_bmk_id==dive_tb_speed$fly_bmk_id[i] & owl_data_speed$fly_bmk_category == "fly"))-1):
                                                                                                                                                        length(which(owl_data_speed$fly_bmk_id==dive_tb_speed$fly_bmk_id[i] & owl_data_speed$fly_bmk_category == "fly")))], na.rm=T)
        
        
      }
    }
  }
  
  dive_tb_tot_speed[dive_tb_tot_speed$TagID==ind,c("speed_before_strike_mean",
                                                   "speed_before_strike_mean_minus2",
                                                   "speed_before_strike_med",
                                                   "speed_before_strike_med_minus2")] <- dive_tb_speed[,c("speed_before_strike_mean",
                                                                                                          "speed_before_strike_mean_minus2",
                                                                                                          "speed_before_strike_med",
                                                                                                          "speed_before_strike_med_minus2")] #
}


boxplot(dive_tb_tot_speed$speed_before_strike_med~dive_tb_tot_speed$Sex.y, ylim=c(5,6))
hist(dive_tb_tot_speed$speed_before_strike[which(dive_tb_tot_speed$speed_before_strike<20)])

##########

measurement<- fread("/Users/kschalch/Downloads/BarnOwls_Legacy_20231114121253_BirdMeasurement.csv") %>% 
  dplyr::filter(GrowthStage=="Adult") %>% 
  group_by(RingId) %>% 
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

head(measurement)

dive_tb_tot_speed$hunting_strategy_final

tb_mod_flightspeed<-dive_tb_tot_speed %>% 
  dplyr::filter(speed_before_strike_mean<20) %>% 
  left_join(measurement)

mod_flightspeed<-lmer(speed_before_strike_med ~ Sex.y +
                        (1|TagID/night_nb_real),data = tb_mod_flightspeed)
summary(mod_flightspeed)
tab_model(mod_flightspeed)
