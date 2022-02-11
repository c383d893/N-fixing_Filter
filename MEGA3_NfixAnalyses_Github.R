############################
###### LOAD PACKAGES #######
############################

library(lme4);library(lmerTest);library(lsmeans);library(ggeffects);library(spdep);library(ggplot2);library(ncf);library(ape);library(sjPlot);library(gridExtra);library(MuMIn);library(tidyverse)
options(na.action = "na.fail") 

#############################
#############################
######### ALL SPECIES #######
#############################
#############################

#############################
######### MODEL S1 ##########
########### BROAD ###########
#############################

gdat<-read.table("data/Nfixdrivers_speciesall_10.27.21.csv") %>%
  filter(!entity_class== "undetermined") %>%                                         
  filter(!nfix==0) %>%                                                               
  filter(!nfixno==0) %>%                                                            
  select(c("entity_ID","entity_class","nfix","nfixno","latitude","longitude","geology", "area",  
           "elev_range","CHELSA_annual_mean_Temp","CHELSA_annual_Prec","dist")) %>%    
  mutate(entity_class2 = case_when(geology=="dev" ~ "oceanic",                    
                                   geology=="nondev" ~ "nonoceanic",
                                   entity_class =="Mainland" ~ "mainland")) %>%
  mutate(elev_range = ifelse(elev_range== 0,1,elev_range)) %>%                                        
  mutate(elev_range = ifelse(is.na(elev_range),1,elev_range)) %>%                    
  select(-geology) %>%                                                              
  drop_na()   

gdat$area <-scale(log10((gdat$area)+.01))                   
gdat$elev_range <- scale(log10(gdat$elev_range+.01))        
gdat$CHELSA_annual_mean_Temp<-scale(gdat$CHELSA_annual_mean_Temp)
gdat$CHELSA_annual_Prec<-scale(gdat$CHELSA_annual_Prec)
gdat$abslatitude<-scale(abs(gdat$latitude))
gdat$squaredlat<-scale(abs(gdat$latitude)^2)

mod.S1<- glmer(cbind(nfix,nfixno)~entity_class2+abslatitude+squaredlat+ area+elev_range+CHELSA_annual_Prec+CHELSA_annual_mean_Temp+(1|entity_class2:entity_ID) , data =gdat, family=binomial(link ="logit"),  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(mod.S1)

#############################
########### PLOT ############
#############################

colScale <- scale_colour_manual(values=c("darkseagreen3","deepskyblue4"))
fillScale <- scale_fill_manual(values=c("darkseagreen3","deepskyblue4"))

ref<-lsmeans(mod.S1,pairwise~entity_class2, data= gdat)
ref.table<-as.data.frame(ref$lsmeans)
ref.table$clsmeans<-exp(ref.table$lsmean)
ref.table$plus.se<-(exp(ref.table$lsmean+ref.table$SE))-ref.table$clsmeans
ref.table$minus.se<-ref.table$clsmeans-(exp(ref.table$lsmean-ref.table$SE))

ref.table<-ref.table[c(1,3),] #remove non oceanic row
ref.table$LandType<-c("Mainland","Island")

gdat2cat<- gdat %>% mutate(propnfix=nfix/nfixno) %>% filter(!entity_class2=="nonoceanic")

fig.1A<- 
  ggplot(ref.table, aes(x=entity_class2, y=clsmeans,fill=entity_class2,color=entity_class2))  + 
  geom_point(position=position_dodge(1), size =2) + 
  geom_errorbar(aes(ymin=clsmeans-minus.se, ymax=clsmeans+plus.se), width=0.4,size=1, position=position_dodge(1)) + 
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  guides(fill=guide_legend(title="Land Type")) +
  geom_bar(stat="identity", alpha=0.3)+
  geom_point(data=gdat2cat,aes(x=entity_class2,y=propnfix,color=entity_class2),position=position_jitter(width =.2),size=2,alpha=0.3)+
  coord_cartesian(ylim=c(0,.2))+
  colScale+
  fillScale+
  labs(y="Proportion N-fixing Plant Species")+
  theme_classic(base_size = 20) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position = "none")+
  aes(x = fct_inorder(c("mainland","oceanic"))) +
  annotate("text", x = .8, y = 0.18, label = "A.", cex=7)

#############################
########### PLOT ############
#############################

colScale <- scale_colour_manual(values=c("darkseagreen3","cyan4","deepskyblue4"))
fillScale <- scale_fill_manual(values=c("darkseagreen3","cyan4","deepskyblue4"))

ref<-lsmeans(mod.S1,pairwise~entity_class2, data= gdat)
ref.table<-as.data.frame(ref$lsmeans)
ref.table$clsmeans<-exp(ref.table$lsmean)
ref.table$plus.se<-(exp(ref.table$lsmean+ref.table$SE))-ref.table$clsmeans
ref.table$minus.se<-ref.table$clsmeans-(exp(ref.table$lsmean-ref.table$SE))

ref.table$LandType<-c("Mainland","Oceanic Island","Non-Oceanic Island")

gdat<- gdat %>% mutate(propnfix=nfix/nfixno) 

fig.S1<- 
  ggplot(ref.table, aes(x=entity_class2, y=clsmeans,fill=entity_class2,color=entity_class2))  + 
  geom_point(position=position_dodge(1), size =2) + 
  geom_errorbar(aes(ymin=clsmeans-minus.se, ymax=clsmeans+plus.se), width=0.4,size=1, position=position_dodge(1)) + 
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  guides(fill=guide_legend(title="Land Type")) +
  geom_bar(stat="identity", alpha=0.3)+
  geom_point(data=gdat,aes(x=entity_class2,y=propnfix,color=entity_class2),position=position_jitter(width =.2),size=2,alpha=0.3)+
  coord_cartesian(ylim=c(0,.2))+
  colScale+
  fillScale+
  labs(y="Proportion N-fixing Plant Species")+
  theme_classic(base_size = 20) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position = "none")+
  aes(x = fct_inorder(c("mainland","nonoceanic","oceanic"))) +
  annotate("text", x = .8, y = 0.18, label = "A.", cex=7)

#############################
######### MODEL S2 ##########
###### ISLAND PRESENCE ######
#############################

gdat<-read.table("data/Nfixdrivers_speciesall_10.27.21.csv") %>%
  filter(!entity_class== "undetermined") %>%                                          
  select(c("entity_ID","entity_class","nfix","nfixno","latitude","longitude","geology", "area",  
           "elev_range","CHELSA_annual_mean_Temp","CHELSA_annual_Prec","dist")) %>%    
  mutate(entity_class2 = case_when(geology=="dev" ~ "oceanic",                        
                                   geology=="nondev" ~ "nonoceanic",
                                   entity_class =="Mainland" ~ "mainland")) %>%
  mutate(elev_range = ifelse(elev_range==0,1,elev_range)) %>%                                          
  mutate(elev_range = ifelse(is.na(elev_range),1,elev_range)) %>%                    
  select(-geology) %>%                                                                
  drop_na() %>%                                                                     
  mutate(hasnfix= ifelse(nfix>0, 1, 0)) %>%                                          
  filter(entity_class2=="oceanic")                                                 

#for model
gdat$area <-scale(log10(gdat$area + .01))                         
gdat$dist<-scale(gdat$dist)                                       
gdat$elev_range <- scale(log10(gdat$elev_range+.01))         
gdat$CHELSA_annual_mean_Temp<-scale(gdat$CHELSA_annual_mean_Temp)
gdat$CHELSA_annual_Prec<-scale(gdat$CHELSA_annual_Prec)
gdat$abslatitude<-scale(abs(gdat$latitude))
gdat$squaredlat<-scale(abs(gdat$latitude)^2)

#for figs
gdat$area <-log10(gdat$area + .01) 
gdat$dist<-gdat$dist
gdat$elev_range <- log10(gdat$elev_range+.01)
gdat$CHELSA_annual_mean_Temp<-gdat$CHELSA_annual_mean_Temp
gdat$CHELSA_annual_Prec<-gdat$CHELSA_annual_Prec
gdat$abslatitude<-abs(gdat$latitude)
gdat$squaredlat<-abs(gdat$latitude)^2

mod.S2<- glm(hasnfix ~  area+dist+CHELSA_annual_Prec +abslatitude+squaredlat+
                CHELSA_annual_mean_Temp +elev_range, data = gdat, family = binomial(link="logit"))
summary(mod.S2)

N <- nrow(gdat)
p <- length(coef(mod.S2))
E1 <- resid(mod.S2, type = "pearson")
Dispersion <- sum(E1^2)/ (N-p)

q.mod.S2 <- glm(hasnfix ~   area+dist+ CHELSA_annual_Prec+abslatitude+squaredlat+
                   CHELSA_annual_mean_Temp +elev_range, data = gdat, family = quasibinomial(link="logit"))
summary(q.mod.S2)

coords <- cbind(gdat$longitude, gdat$latitude)
matrix.dist = as.matrix(dist(cbind(gdat$longitude, gdat$latitude)))
matrix.dist[1:10, 1:10]
matrix.dist.inv <- 1/matrix.dist
matrix.dist.inv[1:10, 1:10]
diag(matrix.dist.inv) <- 0
matrix.dist.inv[1:10, 1:10]
nonsp_moranI = Moran.I(resid(mod.S2), matrix.dist.inv)
myDist = 1000
rac.a1 <- autocov_dist(resid(q.mod.S2), coords, nbs = myDist, type = "inverse", zero.policy = TRUE, style = "W", longlat=T) 
sp.mod.S2<- glm(hasnfix~  area +dist+ CHELSA_annual_Prec+abslatitude+squaredlat+
                    CHELSA_annual_mean_Temp +elev_range+ rac.a1, data = gdat, family = binomial(link="logit"))
summary(sp.mod.S2)
sp_moranI = Moran.I(resid(sp.mod.S2), matrix.dist.inv)

#############################
########### PLOT ############
#############################

gdat$rac.a1 <- rac.a1
new.data <- as.data.frame(t(apply(gdat[,c("area" , "dist", "CHELSA_annual_Prec", "abslatitude", "squaredlat","CHELSA_annual_mean_Temp", "elev_range", "rac.a1","hasnfix")], 2,function(x) mean(x, na.rm=TRUE))))
new.data.dist     <- new.data[rep(1,100),]
dist.range<- seq(from=quantile(gdat$dist, 0.001),to=quantile(gdat$dist, 0.999),length.out=100) # set for pres below too!
new.data.dist$dist <- dist.range
predi.dist.pres    <-predict(sp.mod.S2,type="response",se.fit=TRUE,newdata=new.data.dist) %>% 
  as.data.frame() %>% mutate(dist=dist.range)

fig2.A<-
  ggplot(data = predi.dist.pres, mapping = aes(x = dist, y = fit))+
  geom_line(color= "cyan4") +
  geom_ribbon(aes(ymin = fit-se.fit, ymax = fit+se.fit),fill= "cyan4", alpha = .6)+
  labs(x="Distance (km)", cex=5)+
  labs(y="Probability N-fixing Plant Species")+
  theme_classic(base_size = 20) +
  theme(legend.position = "none")+
  annotate("text", x = 500, y = 0.97, label = "A.", cex=7)

#############################
######### MODEL S3 ##########
##### ISLAND PROPORTION #####
#############################

gdat<-read.table("data/Nfixdrivers_speciesall_10.27.21.csv") %>%
  filter(!entity_class== "undetermined") %>%                                        
  filter(!nfix==0) %>%                                                               
  filter(!nfixno==0) %>%                                                             
  select(c("entity_ID","entity_class","nfix","nfixno","latitude","longitude","geology", "area",   
           "elev_range","CHELSA_annual_mean_Temp","CHELSA_annual_Prec","dist")) %>%    
  mutate(entity_class2 = case_when(geology=="dev" ~ "oceanic",                      
                                   geology=="nondev" ~ "nonoceanic",
                                   entity_class =="Mainland" ~ "mainland")) %>%
  mutate(elev_range = ifelse(elev_range==0,1,elev_range)) %>%                                          
  mutate(elev_range = ifelse(is.na(elev_range),1,elev_range)) %>%                    
  select(-geology) %>%                                                              
  drop_na() %>%                                                                     
  filter(entity_class2=="oceanic")                                                  

#for model
gdat$area <-scale(log10(gdat$area + .01)) #log of area
gdat$dist<-scale(gdat$dist)
gdat$elev_range <- scale(log10(gdat$elev_range+.01))
gdat$CHELSA_annual_mean_Temp<-scale(gdat$CHELSA_annual_mean_Temp)
gdat$CHELSA_annual_Prec<-scale(gdat$CHELSA_annual_Prec)
gdat$abslatitude<-scale(abs(gdat$latitude))
gdat$squaredlat<-scale(abs(gdat$latitude)^2)

#for figs
gdat$area <-log10(gdat$area + .01) #log of area
gdat$dist<-gdat$dist
gdat$elev_range <- log10(gdat$elev_range+.01)
gdat$CHELSA_annual_mean_Temp<-gdat$CHELSA_annual_mean_Temp
gdat$CHELSA_annual_Prec<-gdat$CHELSA_annual_Prec
gdat$abslatitude<-abs(gdat$latitude)
gdat$squaredlat<-abs(gdat$latitude)^2

mod.S3<- nfixgdat.noage <- glm(cbind(gdat$nfix,gdat$nfixno) ~  area*dist+CHELSA_annual_Prec +abslatitude+squaredlat+
                                   CHELSA_annual_mean_Temp +elev_range, data = gdat, family = binomial(link="logit"))
summary(mod.S3) 

N <- nrow(gdat)
p <- length(coef(mod.S3))
E1 <- resid(mod.S3, type = "pearson")
Dispersion <- sum(E1^2)/ (N-p)

q.mod.S3 <- glm(cbind(gdat$nfix,gdat$nfixno) ~   area*dist+ CHELSA_annual_Prec+abslatitude+squaredlat+
                   CHELSA_annual_mean_Temp +elev_range, data = gdat, family = quasibinomial(link="logit"))
summary(q.mod.S3) 

coords <- cbind(gdat$longitude, gdat$latitude)
matrix.dist = as.matrix(dist(cbind(gdat$longitude, gdat$latitude)))
matrix.dist[1:10, 1:10]
matrix.dist.inv <- 1/matrix.dist
matrix.dist.inv[1:10, 1:10]
diag(matrix.dist.inv) <- 0
matrix.dist.inv[1:10, 1:10]
nonsp_moranI = Moran.I(resid(mod.S3), matrix.dist.inv)
myDist = 1000
rac.a1 <- autocov_dist(resid(q.mod.S3), coords, nbs = myDist, type = "inverse", zero.policy = TRUE, style = "W", longlat=T) 
sp.mod.S3 <- glm(cbind(gdat$nfix,gdat$nfixno)~  area*dist+ CHELSA_annual_Prec+abslatitude+squaredlat+
                    CHELSA_annual_mean_Temp +elev_range+ rac.a1, data = gdat, family = binomial(link="logit"))
summary(sp.mod.S3)
sp_moranI = Moran.I(resid(sp.mod.S3), matrix.dist.inv)

#############################
########### PLOT ############
#############################

datplot<- ggeffect(sp.mod.S3, terms=c("dist","area")) %>%
  rename(dist=x,islandsize=group)%>%
  mutate(islandsize= as.numeric(as.character(islandsize)))

datplot$islandsize[datplot$islandsize==0.53]<-"Small"
datplot$islandsize[datplot$islandsize==1.57]<-"Medium"
datplot$islandsize[datplot$islandsize==2.62]<-"Large"

gdat<-gdat %>% mutate(propnfix=(nfix/nfixno)) %>% 
  mutate(islandsize = case_when(area < 0.53 ~ "Small",  
                                area > 2.62 ~ "Large")) %>%
  mutate(islandsize = replace_na(islandsize, "Medium"))

#Create a custom color scale
colScale <- scale_colour_manual(values= c("darkcyan","darkslategray3", "lightblue1"))
fillScale <- scale_fill_manual(values= c("darkcyan","darkslategray3", "lightblue1"))

fig2.B<-
  ggplot(data = datplot,aes(x=dist, y=predicted, group=islandsize,color=islandsize)) +
  geom_line(show.legend=FALSE)+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = islandsize), alpha = .6,color =NA)+
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  coord_cartesian(ylim=c(0,.2))+
  labs(x="Distance (km)", cex=5)+
  labs(y="Proportion N-fixing Plant Species")+
  theme_classic(base_size = 20) +
  guides(fill=guide_legend(title="Land Type"))+
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  colScale +
  fillScale + 
  guides(fill=guide_legend(title="Island Area"))+
  annotate("text", x = 500, y = .175, label = "B.", cex=7)

#############################
######### MODEL S4 ##########
######### SYNERGISM #########
#############################

simp<-read.table("data/Nfixbroad_genusMYC_speciesall_12.01.21.csv") %>%
  filter(!entity_class== "undetermined") %>%    # Remove undetermined
  filter(!count == 0) %>%                       # Remove zero counts
  select(c("count","entity_ID","entity_class","nfixmyc","latitude","geology", "area","dist",  # Select columns
           "elev_range","CHELSA_annual_mean_Temp","CHELSA_annual_Prec")) %>%
  mutate(entity_class2 = case_when(geology=="dev" ~ "oceanic",                      # Create entity_class2
                                   geology=="nondev" ~ "nonoceanic",
                                   entity_class =="Mainland" ~ "mainland")) %>%
  mutate(elev_range = ifelse(elev_range==0,1,elev_range)) %>%                       # Make 0 elevations 1                  
  mutate(elev_range = ifelse(is.na(elev_range),1,elev_range)) %>%                   # Make unknown elevations 1                  
  mutate(logcounts1 = log(count+0.01))  %>%                                         # Log counts nfix
  separate(nfixmyc, c("myc", "nfix"), sep="_", remove=FALSE) %>%                    # Make myc and nfix cols
  select(-geology) %>%                                                              # Remove geology
  drop_na()                                                                         # Remove rows with NA

simp<- simp %>%
  filter(!entity_class2=="nonoceanic")

simp$area <-scale(log10((simp$area)+.01))                   
simp$elev_range <- scale(log10(simp$elev_range+.01))        
simp$CHELSA_annual_mean_Temp<-scale(simp$CHELSA_annual_mean_Temp)
simp$CHELSA_annual_Prec<-scale(simp$CHELSA_annual_Prec)
simp$abslatitude<-scale(abs(simp$latitude))
simp$squaredlat<-scale(abs(simp$latitude)^2)

mod.S4 <- glmer(count~nfix*myc*entity_class2 +(1|entity_class2:entity_ID) +(1|entity_class2:entity_ID:nfix)+(1|entity_class2:entity_ID:myc)+(1|entity_class2:entity_ID:nfix:myc),weights = logcounts1,family = poisson, data = simp,  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) 
summary(mod.S4 )

ref.1<-lsmeans(mod.S4,pairwise~nfix*myc*entity_class2, data=simp)
ref.table.1<-as.data.frame(ref.1$lsmeans)
ref.table.1$clsmeans<-exp(ref.table.1$lsmean)
ref.table.1$plus.se<-(exp(ref.table.1$lsmean+ref.table.1$SE))-ref.table.1$clsmeans
ref.table.1$minus.se<-ref.table.1$clsmeans-(exp(ref.table.1$lsmean-ref.table.1$SE))

ref.table.1 <- ref.table.1 %>% 
  group_by(nfix,myc) %>%
  mutate(clsmeans.rel = clsmeans/clsmeans[entity_class2 == "mainland"]) %>%
  mutate(clsmeans.rel.diff= -(clsmeans.rel))

Fig.S2<-
  ggplot(ref.table.1 %>% filter(entity_class2=="oceanic"|entity_class2=="nonoceanic"), aes(x=entity_class2, y=clsmeans.rel.diff,fill=entity_class2)) +
  geom_bar(stat="identity", alpha=0.8, fill="darkgrey")+
  coord_cartesian(ylim=c(-0.1,0))+
  labs(y="Relative Change in Proportion from Mainland to Island")+
  theme_classic(base_size = 20) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.spacing.x = unit(0, "null"),
        strip.background = element_blank())+
  geom_abline(intercept = 0, slope = 0, linetype="solid")+
  theme(legend.position = "none")+
  facet_grid(~nfix*myc)

#############################
######### MODEL S5 ##########
#### ENDEMISM PROPORTION ####
#############################

gdat<-read.table("data/Nfixdrivers_ENDspeciesall_10.27.21.csv") %>%
  filter(!entity_class== "undetermined") %>%                                         
  filter(!endem==0) %>%                                                              
  filter(!noendem==0) %>%                                                             
  select(c("entity_ID","nfix","latitude","longitude","geology", "area","elev_range","dist",
           "entity_class","CHELSA_annual_Prec","CHELSA_annual_mean_Temp","endem","noendem")) %>%    
  mutate(entity_class2 = case_when(geology=="dev" ~ "oceanic",                       
                                   geology=="nondev" ~ "nonoceanic",
                                   entity_class =="Mainland" ~ "mainland")) %>%
  mutate(elev_range = ifelse(elev_range==0,1,elev_range)) %>%                                      
  mutate(elev_range = ifelse(is.na(elev_range),1,elev_range)) %>%                     
  select(-geology) %>%                                                               
  drop_na() %>%                                                                      
  filter(entity_class2=="oceanic")  

#for models
gdat$area <-scale(log10(gdat$area + .01))
gdat$dist<-scale(gdat$dist)
gdat$elev_range <- scale(log10(gdat$elev_range+.01))
gdat$CHELSA_annual_mean_Temp<-scale(gdat$CHELSA_annual_mean_Temp)
gdat$CHELSA_annual_Prec<-scale(gdat$CHELSA_annual_Prec)
gdat$abslatitude<-scale(abs(gdat$latitude))
gdat$squaredlat<-scale(abs(gdat$latitude)^2)

#for figs
gdat$area <-gdat$area
gdat$dist<-gdat$dist
gdat$elev_range <- gdat$elev_range
gdat$CHELSA_annual_mean_Temp<-gdat$CHELSA_annual_mean_Temp
gdat$CHELSA_annual_Prec<-gdat$CHELSA_annual_Prec
gdat$abslatitude<-abs(gdat$latitude)
gdat$squaredlat<-abs(gdat$latitude)^2

mod.S5 <- nfixgdat.noage <- glm(cbind(endem,noendem) ~  nfix*dist+ nfix*area +CHELSA_annual_Prec +abslatitude+squaredlat+
                                     CHELSA_annual_mean_Temp +elev_range, data = gdat, family = binomial(link="logit"))
summary(mod.S5) 

N <- nrow(gdat)
p <- length(coef(mod.S5))
E1 <- resid(mod.S5, type = "pearson")
Dispersion <- sum(E1^2)/ (N-p)

q.mod.S5 <- glm(cbind(endem,noendem) ~ nfix*dist+ nfix*area+CHELSA_annual_Prec+abslatitude+squaredlat+
                     CHELSA_annual_mean_Temp +elev_range, data = gdat, family = quasibinomial(link="logit"))
summary(q.mod.S5) 

coords <- cbind(gdat$longitude, gdat$latitude)+matrix(runif(2*nrow(gdat), 0, 0.00001), nrow = nrow(gdat), ncol = 2)
matrix.dist = as.matrix(dist(coords))
matrix.dist[1:10, 1:10]
matrix.dist.inv <- 1/matrix.dist
matrix.dist.inv[1:10, 1:10]
diag(matrix.dist.inv) <- 0
matrix.dist.inv[1:10, 1:10]
nonsp_moranI = Moran.I(resid(mod.S5), matrix.dist.inv)
myDist = 1000
rac.a1 <- autocov_dist(resid(q.mod.S5), coords, nbs = myDist, type = "inverse", zero.policy = TRUE, style = "W", longlat=T)
sp.mod.S5 <- glm(cbind(endem,noendem)~  nfix*dist+ nfix*area + CHELSA_annual_Prec+abslatitude+squaredlat+
                      CHELSA_annual_mean_Temp +elev_range+ rac.a1, data = gdat, family = binomial(link="logit"))
summary(sp.mod.S5)

#############################
########### PLOT ############
#############################

datplot<- ggeffect(sp.mod.S5, terms=c("area","nfix")) %>%
  mutate(nfix = case_when(group=="nfix" ~ "N fixing",                       
                          group== "nfixno" ~ "Non N-fixing")) %>%
  rename(area=x)

colScale <- scale_colour_manual(values= c("mediumpurple1","brown3"))
fillScale <- scale_fill_manual(values= c("mediumpurple1","brown3"))

fig.S3A<-
  ggplot(data = datplot,aes(x=area, y=predicted, group=nfix, colour=nfix)) +
  geom_line()+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = nfix), alpha = .6,color =NA)+
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  coord_cartesian(ylim=c(0,1.3))+
  labs(x="Area (squared km)", cex=5)+
  labs(y="Proportion Endemic Plant Species")+
  theme_classic(base_size = 20) +
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  colScale +
  fillScale + 
  guides(fill=guide_legend(title="N-fixing status"))+
  theme(legend.position = "none")+
  annotate("text", x = 2500, y = 1.2, label = "A.", cex=7)

##########################
##########################
######## LEGUMES #########
##########################
##########################

#############################
######### MODEL S6 ##########
########### BROAD ###########
#############################

gdat<-read.table("data/Nfixdrivers_Legume.speciesall_10.27.21.csv") %>%
  filter(!entity_class== "undetermined") %>%                                         
  filter(!nfix==0) %>%                                                              
  filter(!nfixno==0) %>%                                                            
  select(c("entity_ID","entity_class","nfix","nfixno","latitude","longitude","geology", "area",   
           "elev_range","CHELSA_annual_mean_Temp","CHELSA_annual_Prec","dist")) %>%    
  mutate(entity_class2 = case_when(geology=="dev" ~ "oceanic",                     
                                   geology=="nondev" ~ "nonoceanic",
                                   entity_class =="Mainland" ~ "mainland")) %>%
  mutate(elev_range = ifelse(elev_range==0,1,elev_range)) %>%                                         
  mutate(elev_range = ifelse(is.na(elev_range),1,elev_range)) %>%                    
  select(-geology) %>%                                                               
  drop_na()                                                                          

gdat$area <-scale(log10((gdat$area)+.01))                   
gdat$elev_range <- scale(log10(gdat$elev_range+.01))        
gdat$CHELSA_annual_mean_Temp<-scale(gdat$CHELSA_annual_mean_Temp)
gdat$CHELSA_annual_Prec<-scale(gdat$CHELSA_annual_Prec)
gdat$abslatitude<-scale(abs(gdat$latitude))
gdat$squaredlat<-scale(abs(gdat$latitude)^2)

mod.S6<- glmer(cbind(gdat$nfix,gdat$nfixno)~entity_class2+abslatitude +  
                     squaredlat + area + elev_range + CHELSA_annual_Prec + CHELSA_annual_mean_Temp +(1|entity_class2:entity_ID) , data =gdat, family=binomial(link ="logit"))
summary(mod.S6)

#############################
########### PLOT ############
#############################

colScale <- scale_colour_manual(values=c("darkseagreen3","cyan4","deepskyblue4"))
fillScale <- scale_fill_manual(values=c("darkseagreen3","cyan4","deepskyblue4"))

ref<-lsmeans(mod.S6,pairwise~entity_class2, data= gdat)
ref.table<-as.data.frame(ref$lsmeans)
ref.table$clsmeans<-exp(ref.table$lsmean)
ref.table$plus.se<-(exp(ref.table$lsmean+ref.table$SE))-ref.table$clsmeans
ref.table$minus.se<-ref.table$clsmeans-(exp(ref.table$lsmean-ref.table$SE))

ref.table$LandType<-c("Mainland","Oceanic Island","Non-Oceanic Island")

gdat<- gdat %>% mutate(propnfix=nfix/nfixno)

fig.S1B<- 
  ggplot(ref.table, aes(x=entity_class2, y=clsmeans,fill=entity_class2,color=entity_class2))  + 
  geom_point(position=position_dodge(1), size =2) + 
  geom_errorbar(aes(ymin=clsmeans-minus.se, ymax=clsmeans+plus.se), width=0.4,size=1, position=position_dodge(1)) + 
  geom_bar(stat="identity", alpha=0.3)+
  geom_point(data=gdat,aes(x=entity_class2,y=propnfix,color=entity_class2),position=position_jitter(width =.2),size=2,alpha=0.3)+
  coord_cartesian(ylim=c(0,30))+
  colScale+
  fillScale+
  labs(y="Proportion N-fixing Legumes")+
  theme_classic(base_size = 20) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  aes(x = fct_inorder(c("mainland","nonoceanic","oceanic"))) +
  theme(legend.position = "none")+
  annotate("text", x = .8, y = 28, label = "B.", cex=7)

#############################
######### MODEL S7 ##########
##### ISLAND PROPORTION #####
#############################

gdat<-read.table("data/Nfixdrivers_Legume.speciesall_10.27.21.csv") %>%
  filter(!entity_class== "undetermined") %>%                                         
  filter(!nfix==0) %>%                                                               
  filter(!nfixno==0) %>%                                                             
  select(c("entity_ID","entity_class","nfix","nfixno","latitude","longitude","geology", "area",   
           "elev_range","CHELSA_annual_mean_Temp","CHELSA_annual_Prec","dist")) %>%    
  mutate(entity_class2 = case_when(geology=="dev" ~ "oceanic",                       
                                   geology=="nondev" ~ "nonoceanic",
                                   entity_class =="Mainland" ~ "mainland")) %>%
  mutate(elev_range = ifelse(elev_range==0,1,elev_range)) %>%                                        
  mutate(elev_range = ifelse(is.na(elev_range),1,elev_range)) %>%                    
  select(-geology) %>%                                                               
  drop_na() %>%                                                                   
  filter(entity_class2=="oceanic")                                                   

gdat$area <-scale(log10(gdat$area + .01)) 
gdat$dist<-scale(gdat$dist)
gdat$elev_range <- scale(log10(gdat$elev_range+.01))
gdat$CHELSA_annual_mean_Temp<-scale(gdat$CHELSA_annual_mean_Temp)
gdat$CHELSA_annual_Prec<-scale(gdat$CHELSA_annual_Prec)
gdat$abslatitude<-scale(abs(gdat$latitude))
gdat$squaredlat<-scale(abs(gdat$latitude)^2)

mod.S7 <- nfixgdat.noage <- glm(cbind(gdat$nfix,gdat$nfixno) ~  area*dist+CHELSA_annual_Prec +abslatitude+squaredlat+
                                   CHELSA_annual_mean_Temp +elev_range, data = gdat, family = binomial(link="logit"))
summary(mod.S7) 

N <- nrow(gdat)
p <- length(coef(mod.S7))
E1 <- resid(mod.S7, type = "pearson")
Dispersion <- sum(E1^2)/ (N-p)

q.mod.S7 <- glm(cbind(gdat$nfix,gdat$nfixno) ~   area*dist+ CHELSA_annual_Prec+abslatitude+squaredlat+
                   CHELSA_annual_mean_Temp +elev_range, data = gdat, family = quasibinomial(link="logit"))
summary(q.mod.S7) 

coords <- cbind(gdat$longitude, gdat$latitude)
matrix.dist = as.matrix(dist(cbind(gdat$longitude, gdat$latitude)))
matrix.dist[1:10, 1:10]
matrix.dist.inv <- 1/matrix.dist
matrix.dist.inv[1:10, 1:10]
diag(matrix.dist.inv) <- 0
matrix.dist.inv[1:10, 1:10]
nonsp_moranI = Moran.I(resid(mod.S7), matrix.dist.inv)
myDist = 1000
rac.a1 <- autocov_dist(resid(q.mod.S7), coords, nbs = myDist, type = "inverse", zero.policy = TRUE, style = "W", longlat=T)
sp.mod.S7 <- glm(cbind(gdat$nfix,gdat$nfixno)~  area*dist+ CHELSA_annual_Prec+abslatitude+squaredlat+
                    CHELSA_annual_mean_Temp +elev_range+ rac.a1, data = gdat, family = binomial(link="logit"))
summary(sp.mod.S7)
sp_moranI = Moran.I(resid(sp.mod.S7), matrix.dist.inv)

#############################
######### MODEL S8 ##########
#### ENDEMISM PROPORTION ####
#############################

gdat<-read.table("data/Nfixdrivers_ENDLegumespeciesall_10.27.21.csv") %>%
  filter(!entity_class== "undetermined") %>%                                         
  filter(!endem==0) %>%                                                               
  filter(!noendem==0) %>%                                                             
  select(c("entity_ID","nfix","latitude","longitude","geology", "area","elev_range","dist",
           "entity_class","CHELSA_annual_Prec","CHELSA_annual_mean_Temp","endem","noendem")) %>%    
  mutate(entity_class2 = case_when(geology=="dev" ~ "oceanic",                     
                                   geology=="nondev" ~ "nonoceanic",
                                   entity_class =="Mainland" ~ "mainland")) %>%
  mutate(elev_range = ifelse(elev_range==0,1,elev_range)) %>%                                         
  mutate(elev_range = ifelse(is.na(elev_range),1,elev_range)) %>%                    
  select(-geology) %>%                                                             
  drop_na() %>%                                                                      
  filter(entity_class2=="oceanic")  

#for models
gdat$area <-scale(log10(gdat$area + .01))
gdat$dist<-scale(gdat$dist)
gdat$elev_range <- scale(log10(gdat$elev_range+.01))
gdat$CHELSA_annual_mean_Temp<-scale(gdat$CHELSA_annual_mean_Temp)
gdat$CHELSA_annual_Prec<-scale(gdat$CHELSA_annual_Prec)
gdat$abslatitude<-scale(abs(gdat$latitude))
gdat$squaredlat<-scale(abs(gdat$latitude)^2)

#for figs
gdat$area <-gdat$area + .01 
gdat$dist<-gdat$dist
gdat$elev_range <- log10(gdat$elev_range+.01)
gdat$CHELSA_annual_mean_Temp<-gdat$CHELSA_annual_mean_Temp
gdat$CHELSA_annual_Prec<-gdat$CHELSA_annual_Prec
gdat$abslatitude<-abs(gdat$latitude)
gdat$squaredlat<-abs(gdat$latitude)^2

mod.S8 <- nfixgdat.noage <- glm(cbind(endem,noendem) ~nfix*area + nfix*dist+CHELSA_annual_Prec +abslatitude+squaredlat+
                                     CHELSA_annual_mean_Temp +elev_range, data = gdat, family = binomial(link="logit"))
summary(mod.S8) 

N <- nrow(gdat)
p <- length(coef(mod.S8))
E1 <- resid(mod.S8, type = "pearson")
Dispersion <- sum(E1^2)/ (N-p)

q.mod.S8 <- glm(cbind(endem,noendem) ~    nfix*area + nfix*dist+ CHELSA_annual_Prec+abslatitude+squaredlat+
                     CHELSA_annual_mean_Temp +elev_range, data = gdat, family = quasibinomial(link="logit"))
summary(q.mod.S8) 

coords <- cbind(gdat$longitude, gdat$latitude)+matrix(runif(2*nrow(gdat), 0, 0.00001), nrow = nrow(gdat), ncol = 2)
matrix.dist = as.matrix(dist(cbind(gdat$longitude, gdat$latitude)))
matrix.dist[1:10, 1:10]
matrix.dist.inv <- 1/matrix.dist
matrix.dist.inv[1:10, 1:10]
matrix.dist.inv[is.infinite(matrix.dist.inv)] <- 0
matrix.dist.inv[1:10, 1:10]
nonsp_moranI = Moran.I(resid(mod.S8), matrix.dist.inv)
myDist = 1000
rac.a1 <- autocov_dist(resid(q.mod.S8), coords, nbs = myDist, type = "inverse", zero.policy = TRUE, style = "W", longlat=T) 
sp.mod.S8 <- glm(cbind(endem,noendem)~  nfix*dist+ nfix*area + CHELSA_annual_Prec+abslatitude+squaredlat+
                      CHELSA_annual_mean_Temp +elev_range+ rac.a1, data = gdat, family = binomial(link="logit"))
summary(sp.mod.S8)

#############################
########### PLOT ############
#############################

datplot<- ggeffect(sp.mod.S8, terms=c("area","nfix")) %>%
  mutate(nfix = case_when(group=="nfix" ~ "N fixing",                       
                          group== "nfixno" ~ "Non N-fixing")) %>%
  rename(area=x)

colScale <- scale_colour_manual(values= c("mediumpurple1","brown3"))
fillScale <- scale_fill_manual(values= c("mediumpurple1","brown3"))

fig.S3B<-
  ggplot(data = datplot,aes(x=area, y=predicted, group=nfix, colour=nfix)) +
  geom_line(show.legend=FALSE)+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = nfix), alpha = .6,color =NA)+
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  coord_cartesian(ylim=c(0,1.3))+
  labs(x="Area (squared km)", cex=5)+
  labs(y="Proportion Endemic Legumes")+
  theme_classic(base_size = 20) +
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  colScale +
  fillScale + 
  guides(fill=guide_legend(title="N-fixing status"),
         line)+
  annotate("text", x = 2500, y = 1.2, label = "B.", cex=7)
