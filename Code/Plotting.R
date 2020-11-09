###########################################################
###  Code for Rassweiler et al. (2021) Ecol. Appl. 
###  Authors: D.K. Okamoto, A. Rassweiler
###  Last Updated: May, 2018 (commented Nov. 2020)  
###  
###  
###  **CAUTION: code is not designed for efficiency or for
###  other applications, but for this paper alone!**
###########################################################

###  Code provided here is for plotting simulation output only 

### load packages 
library("ggplot2")
library("plyr")
library("dplyr")
library("gdata")
library("reshape2")
library("RColorBrewer")
library("lme4")
library("boot")
library("gridExtra")
library("MBA")
library("fields")
library("cowplot")

### read in raw output and summaries
### download the large file from google drive link below
## https://drive.google.com/file/d/1-FIVbGKsh1y4IUK_0BH-t_gydH1n1oWJ/view?usp=sharing
out <- read.csv("runs_output_full.csv") 

### summarize data across groups
out_sum <- out%>%
  group_by(species,incl,severity,domain)%>%
  dplyr::summarize(log_est= mean(leffect2,na.rm=T),
            log_est2= mean(leffect3,na.rm=T),
            sd_log_est2 = sd(exp(leffect2),na.rm=T),
            sd_log_est3 = sd(exp(leffect3),na.rm=T),
            psc = mean(coef2_p<0.025&leffect2<0&!is.na(baci_p2)),
            ppc = mean(coef2_p<0.025&leffect2>=0&!is.na(baci_p2)),
            ps = mean(baci_p2<0.025&leffect2<0&!is.na(baci_p2)),
            pp = mean(baci_p2<0.025&leffect2>=0&!is.na(baci_p2)),
            ns = mean(baci_p2>=0.025&!is.na(leffect2)&!is.na(baci_p2)),
            nsc = mean(coef2_p>=0.025&!is.na(leffect2)&!is.na(coef2_p)),
            ps2 = mean(baci_p3<0.025&leffect3<0&!is.na(baci_p3)),
            pp2 = mean(baci_p3<0.025&leffect3>=0&!is.na(baci_p3)),
            ns2 = mean(baci_p3>=0.025&!is.na(leffect3)&!is.na(baci_p3)),
            ps2 = mean(baci_p3<0.025&leffect3<0&!is.na(baci_p3)),
            ## all non-na results
            nodata = mean(is.na(ctl_acf)),
            noconv =mean(is.na(baci_p2)&!is.na(ctl_acf)),
            noconv2 =mean(is.na(baci_p3)&!is.na(ctl_acf)),
            n_pres= mean(!is.na(leffect2)),
            n_pres2= mean(!is.na(leffect3)),
            ctl_acf= mean(ctl_acf,na.rm=T),
            imp_acf= mean(imp_acf,na.rm=T),
            ctl_cv= mean(ctl_cv,na.rm=T),
            imp_cv= mean(imp_cv,na.rm=T),
            sv_pre =mean(sv_pre,na.rm=T),
            ctl_mu= mean(ctl_mu,na.rm=T),
            imp_mu= mean(imp_mu,na.rm=T),
            ctl_n= mean(ctl_n,na.rm=T),
            imp_n= mean(imp_n,na.rm=T),
            ctl_n_years= mean(ctl_ny+ctlp_ny,na.rm=T),
            imp_n_years= mean(imp_ny+impp_ny,na.rm=T),
            effect_cv= abs(mean(leffect2_se/leffect2,na.rm=T)),
            effect2_cv= abs(mean(leffect3_se/leffect3,na.rm=T)),
            n= length(baci_p2),
            n2= length(baci_p3)
          )%>%
  mutate(group= ifelse(species%in%c(16:24),"pc","count"))%>%
  data.frame()

### read and merge species codes
codes <- read.csv("Data/Species_Codes.csv")
out_sum <- left_join(out_sum,codes,by= "species")
out_sum <- out_sum%>%
  drop.levels()

### subset by quarter hit
out_full <- out%>%
  filter(domain=="quart")%>%
  group_by(species,incl,severity,domain)%>%
  filter(ctl_n>1&imp_n>1)%>%
  summarize(
         N= length(coef2_p),
         ctl_cv= mean(ctl_cv,na.rm=T),
         ctl_n= mean(ctl_n,na.rm=T),
         imp_n= mean(imp_n,na.rm=T),
         log_est=mean(leffect2,na.rm=T),
         log_est_AR1= mean(leffect3,na.rm=T),
         sd_log_est = sd(exp(leffect2),na.rm=T),
         sd_log_est_AR1 = sd(exp(leffect3),na.rm=T),
         sig_neg_z = mean(coef2_p<0.05&leffect2<0&!is.na(baci_p2),na.rm=T),
         sig_pos_z = mean(coef2_p<0.05&leffect2>=0&!is.na(baci_p2),na.rm=T),
         sig_neg = mean(baci_p2<0.05&leffect2<0&!is.na(baci_p2),na.rm=T),
         sig_pos = mean(baci_p2<0.05&leffect2>=0&!is.na(baci_p2),na.rm=T),
         sig_neg_z_AR1 = mean(coef2_p<0.05&leffect3<0&!is.na(baci_p3),na.rm=T),
         sig_pos_z_AR1 = mean(coef2_p<0.05&leffect3>=0&!is.na(baci_p3),na.rm=T),
         sig_neg_AR1 = mean(baci_p3<0.05&leffect3<0&!is.na(baci_p3),na.rm=T),
         sig_pos_AR1 = mean(baci_p3<0.05&leffect3>=0&!is.na(baci_p3),na.rm=T),
         sig= mean(baci_p2<0.05&!is.na(baci_p2),na.rm=T),
         sig_AR1= mean(baci_p3<0.05&!is.na(baci_p3),na.rm=T),
         pre_synchrony = mean(sv_pre,na.rm=T),
         ctl_autocorr = mean(ctl_acf,na.rm=T))%>%
  left_join(codes[,1:3],by= "species")

### subset for surfperch and macro
dat <- out_sum%>%
  filter(taxa_name%in%c("Black Surfperch","Macrocystis")&
         domain%in%unique(out_sum$domain)[1:24])%>%
  mutate(response= ps+pp)%>%
  select(severity,domain,response,taxa_name)%>%
  mutate(severity= log10(ifelse((1-severity)>1,1-severity,ifelse((1-severity)<1,1-severity,1))),
         domain = as.numeric(as.character(domain)))

radius_n <- out_sum%>%
  filter(taxa_name%in%c("Black Surfperch","Macrocystis")&
         domain%in%c(25,50,75,100))%>%
 group_by(domain)%>%
 summarize(mean_n= mean(imp_n/23))

for(i in 1:2){
  interp_fun <- function(x){  mba.int<- mba.surf(x[,1:3], 50,50, extend=T)$xyz.est
    expand.grid(severity=mba.int[[1]],domain=mba.int[[2]])%>%
    mutate(response= melt(mba.int[[3]])[,3])}
  if(i==1){
    interp <- dat%>%
      filter(taxa_name==c("Black Surfperch","Macrocystis")[i])%>%
      interp_fun()%>%
      mutate(taxa_name= c("Black Surfperch","Macrocystis")[i])
  } else {
    temp <- dat%>%
      filter(taxa_name==c("Black Surfperch","Macrocystis")[i])%>%
      interp_fun()%>%
      mutate(taxa_name= c("Black Surfperch","Macrocystis")[i])
    interp<- rbind(interp,temp)
  }
}

rgb2 <- c("white",tim.colors(9))
rgb.palette <- colorRampPalette(rgb2, space = "rgb",bias= 2)
  
heatmap <- ggplot(aes(y=10^(severity),x=domain,
                      fill= ifelse(response<0,0,response),
                      z= ifelse(response<0,0,response)),
                  data=subset(interp,severity<=0))+
  geom_raster()+
  facet_wrap(~taxa_name)+
  scale_fill_gradientn(colours= rgb.palette(10),
                       name= "Power to detect impact",
                       limits= c(0,0.8))+
  scale_colour_gradient(high= "black",low="black")+
  geom_point(size= 0.1,data= subset(dat,severity<=0))+
  scale_y_continuous(trans= "log10",breaks= c(1,1/2,1/5),labels= c("0","-50%","-80%"),expand= c(0,0))+
  scale_x_continuous(breaks= c(25,50,75,100),expand= c(0,0),
                     labels= c("25\n(9%)",
                               "50\n(20%)",
                               "75\n(55%)",
                               "100\n(79%)"))+
  xlab("Radius of Impact [km]\n (Average % of Sites Impacted)")+
  ylab("Severity of Impact")+
  guides(fill= guide_colourbar(title.position= "top",
                               title.hjust=0.5),
         colour= "none")+
  annotation_logticks(side= "l",base= 10)+
  theme(axis.text.y= element_text(hjust=1),
        legend.position= "top",
        strip.text= element_blank(),
        legend.background = element_blank(),
        legend.key.height= unit(0.25,"cm"),
        legend.key.width= unit(0.8,"cm"),
        strip.background=element_blank(),
        axis.text.x= element_blank(),
        legend.direction= "horizontal",
        axis.title.x= element_blank(),
        panel.background= element_blank(),
        plot.background= element_blank(),
        panel.border= element_rect(colour= "black",fill=NA))+
  geom_hline(yintercept= 1)

#pdf(width=6,height=4.5,file="4_Figures/heatmaps.pdf")
heatmap
#dev.off()

### relative error for BSP and Macro
error_dat <-out%>%
  filter(leffect2>(-5)&leffect2<5&severity==0.8&species%in%c(1,27))%>%
  mutate(species=ifelse(species==1,"Macrocystis","Black Surfperch"))

errorplot <- ggplot(aes(x= as.numeric(as.character(domain)),
  y=(leffect2-log(1-severity))/(log(1-severity))),
  data= error_dat)+
  stat_summary(
    fun.ymin = function(z) { quantile(z,0.95) },
    fun.ymax = function(z) { quantile(z,0.05) },
    fun.y = median, geom= "ribbon",fill= "grey60")+
   stat_summary(
    fun.ymin = function(z) { quantile(z,0.75) },
    fun.ymax = function(z) { quantile(z,0.25) },
    fun.y = median, geom= "ribbon",fill= "grey30")+
   stat_summary(
    fun.y = median, geom= "line",linetype= "dotted")+
  coord_cartesian(ylim= c(-1.5,1.5),expand= c(0,0))+
   scale_x_continuous(breaks= c(25,50,75,100),expand= c(0,0),
                     labels= c("25\n(9%)",
                               "50\n(20%)",
                               "75\n(55%)",
                               "100\n(79%)"))+
  facet_grid(~species)+
  geom_abline(slope=0)+
  ylab("Relative error in \n estimate of effect size")+
  xlab("Impact Radius (km)")+
   theme(axis.text.y= element_text(hjust=1),
          legend.position= "bottom",
          strip.text= element_blank(),
          legend.background = element_blank(),
          legend.key.width= unit(2,"cm"),
          strip.background=element_blank(),
          panel.background= element_blank(),
          plot.background= element_blank(),
          panel.border= element_rect(colour= "black",fill=NA))

#pdf(width=5,height= 6,file="4_Figures/Figure_2.pdf")
plot_grid(heatmap,errorplot,ncol= 1,align= "v",axis= "lr")
#dev.off()

plot_groups <- read.csv("~/Dropbox/DanO_DOI_Partnership/DOI_Impacts/3_Output/runs_output_groups.csv")
plot_groups <- plot_groups%>%
  mutate(incl= 0.15)%>%
  filter(severity!= "4x impact")

plot_dat <- out_sum%>%
  filter(species!=32)%>%
  filter(severity%in%c(-4,0,0.8))%>%
  filter(!(taxa_name%in%c("Blue Rockfish","Blackeye Goby","Striped Surfperch",
    "Pachythyone","Haliotis corrugata","Haliotis rufescens","Rock wrasse")))%>%
  mutate(domain= factor(domain))%>%
  dplyr::select(ps,pp,ns,noconv,nodata,severity,species,taxa_name,domain,incl,group,taxa_group,num_2)%>%
  melt(id.vars= c("severity","species","taxa_name","domain","incl","group","taxa_group","num_2"))%>%
  group_by(severity,domain,variable,taxa_group)%>%
  mutate(mean_val= mean(value))

plot_dat$domain= factor(plot_dat$domain, levels= c(levels(plot_dat$domain)[-c(25:28)][order(as.numeric(as.character(levels(plot_dat$domain)[-c(25:28)])))],levels(plot_dat$domain)[c(25:28)]))
levels(plot_dat$variable) <- c("significant & negative","significant & positive",
  "not significant","no convergence","insufficient data")

plot_dat$severity <- factor(plot_dat$severity, labels= c("+4x impact","no impact","-4x impact"))
plot_dat$severity <- factor(plot_dat$severity, levels= c("-4x impact","no impact","+4x impact"))
taxa_levels <- unique(plot_dat$taxa_name)[order((plot_dat%>%filter(domain=="quart"&
  severity=="no impact"&variable== "significant & negative"))$value)]
plot_dat$group<- factor(plot_dat$group, labels= c("count","point contact"))
plot_dat$taxa_name <- as.character(plot_dat$taxa_name)
plot_dat <- join(plot_dat%>%data.frame(),plot_groups,type= "full")

plot_dat$taxa_name <- factor(
  plot_dat$taxa_name,
  levels= c(c("all algae","all mobile Inverts",
  "all fish","all sessile inverts"),
  as.character(taxa_levels))
)

### read in informedness rank
latin_binom <-   read.csv("Output/latin_binomial_table.csv")
inf_rank <- read.csv("Output/informedness_table.csv")

### join and sort by informedness
plot_dat <- join(plot_dat,inf_rank)%>%
  mutate(taxa_name= as.character(taxa_name),
         binomial= as.character(binomial))%>%
  mutate(taxa_name= ifelse(taxa_name%in%c("all algae","all mobile Inverts","all fish","all sessile inverts"),taxa_name,binomial))

levels(plo)
orders <- unique(plot_dat[,c("taxa_name","rank")])%>%
         mutate(rank= ifelse(is.na(rank),row_number(),rank))%>%
  arrange(rank)

plot_dat$variable = factor(plot_dat$variable, levels= c("significant & negative","significant & positive", 
                                                        "not significant",  "no convergence","insufficient data"))    
plot_dat$taxa_name = factor(plot_dat$taxa_name,levels= orders$taxa_name)


##### plot detections for half and quart (significant and negative)
## remove 150
## plot spatial variation time series
## highlight outliers
## 
plot_dat$severity <- as.factor(plot_dat$severity)
levels(plot_dat$severity)[1] <- "-80% impact"
impact_by_spp <- ggplot(aes(taxa_name,-value),
       data= subset(plot_dat,domain%in%c("quart")&incl==0.15&severity!="+4x impact"))+
  geom_bar(position = "fill", stat = "identity",aes(fill= variable))+
  facet_grid(taxa_group~severity,scales= "free_y",space="free_y")+
  coord_equal()+
  coord_flip()+
  theme_bw()+
  ylab("")+
  guides(fill= guide_legend(ncol= 1),shape=guide_legend(ncol=1))+
  scale_y_continuous(breaks= seq(-1,0,by = 0.1),
                     labels= c(0,"",0.2,"",0.4,"",0.6,"",0.8,"",""),expand= c(0,0))+
  ylab("outcome frequency")+
  scale_shape_manual(values= c(9,22))+
  scale_fill_manual(values= c(rev(brewer.pal(5,"RdYlBu"))))+
  theme(strip.background= element_blank(),
        legend.title= element_blank(),
        legend.position= "right",
        axis.title.y= element_blank(),
        legend.background= element_blank(),
        plot.background= element_blank(),
        legend.key= element_blank())
#pdf(width=7,height= 4.75,file="4_Figures/Figure_4.pdf")
impact_by_spp
#dev.off()

# write.csv("Table_probababilites_quart.csv")

#### another way of looking at impact radius 
summaries <- out%>%
  filter(species!=32&!(domain%in%c("half","half_shuffle","quart","quart_shuffle")))%>%
  group_by(domain,severity,incl)%>%
  dplyr::summarize(frac= mean(imp_n/(imp_n+ctl_n)),
            ps = mean(baci_p2<0.05&leffect2<0&!is.na(baci_p2)),
            pp = mean(baci_p2<0.05&leffect2>=0&!is.na(baci_p2)),
            ns = mean(baci_p2>=0.05&!is.na(leffect2)&!is.na(baci_p2)),
            ps = mean(baci_p2<0.05&leffect2<0&!is.na(baci_p2)),
            ## all non-na results
            nodata = mean(is.na(ctl_acf)),
            noconv =mean((err=="insufficient data"&!is.na(ctl_acf))|(err=="sufficient data"&is.na(baci_p2))))

plot1 <- ggplot(aes(value,x=as.numeric(as.character(domain))),
                data= subset(plot_dat,!(domain%in%c("half","half_shuffle","quart","quart_shuffle"))&
                             variable=="significant & negative"&severity=="-80% impact"))+
  geom_line(aes(group= species),size= 0.25,alpha= 0.5)+
  stat_summary(col= "blue", geom= "line",fun= mean,size= 1)+
  geom_hline(yintercept= 0)+
  theme_bw()+
  scale_y_continuous(limits= c(0,1),breaks= seq(0,1,by = 0.2),expand= c(0,0))+
  scale_x_continuous(breaks= seq(0,125,by= 25))+
  ylab("probability of detection")+
  xlab("impact radius (km)")+
  theme(panel.grid= element_blank(),
        panel.background= element_blank(),
        strip.background= element_blank(),
        legend.title= element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        legend.background= element_blank(),
        plot.background= element_blank());plot1

plot2 <- ggplot(aes(value,x=as.numeric(as.character(domain))),
                data= subset(plot_dat,!(domain%in%c("half","half_shuffle","quart","quart_shuffle"))&
                            variable=="insufficient data"&severity=="-80% impact"))+
  geom_line(aes(group= species),size= 0.25,alpha= 0.5)+
  stat_summary(col= "blue", geom= "line",fun= mean,size= 1)+
  geom_hline(yintercept= 0)+
  theme_bw()+
  scale_y_continuous(limits= c(0,1),breaks= seq(0,1,by = 0.2),expand= c(0,0))+
  scale_x_continuous(breaks= seq(0,125,by= 25))+
  ylab("data insufficiency")+
  xlab("impact radius (km)")+
  theme(panel.grid= element_blank(),
        panel.background= element_blank(),
        strip.background= element_blank(),
        legend.title= element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        legend.background= element_blank(),
        plot.background= element_blank(),
        legend.key= element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank());plot2

plot3 <- ggplot(aes(as.numeric(as.character(domain)),frac*100),data= summaries)+
  stat_summary(fun= mean,geom= "line")+
  ylab("% of sites impacted")+
  xlab("impact radius")+
  theme_bw()+
  scale_y_continuous(limits= c(0,100),breaks= seq(0,100,by = 20),expand= c(0,0))+
  scale_x_continuous(expand= c(0,0),breaks= seq(0,150,by= 25))+
  theme(panel.grid= element_blank(),
        strip.background=element_blank(),
        strip.text= element_blank(),
        panel.background= element_rect(fill= NA),
        legend.title= element_blank(),
        legend.position= c(0.1,0.25),
        axis.text.x=element_blank(),
        axis.title.x=element_blank());plot3

#pdf(width= 5,height= 6.5,file="4_Figures/Figure_3.pdf")
plot_grid(plot3,plot2,plot1,align= "v",ncol=1,rel_heights = c(0.5,0.8,1))
#dev.off()

#### impact radius versus impacted sites and data sufficiency
plot_dat2 <-out_sum%>%
  filter(severity==0&domain%in%c("quart","quart_shuffle"))%>%
  mutate(domain= ifelse(domain=="quart","aggregated","isolated"))

plot_dat3 <-out_sum%>%
  filter(severity%in%c(0.8,0)&domain%in%c("quart","quart_shuffle"))%>%
  mutate(domain= ifelse(domain=="quart","aggregated","isolated"))
           
binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}

plot_AR <- ggplot(aes(ctl_acf,((1-exp(-log_est))-(1-severity))/(1-severity)),data= plot_dat3)+
  geom_point(aes(y= ((1-exp(-log_est))-(1-severity))/(1-severity),colour= "AR(1) error model"))+
  geom_point(aes(y= ((1-exp(-log_est2))-(1-severity))/(1-severity),colour="white noise error model"))+
  stat_smooth(se=F,aes(colour="white noise error model"))+
  stat_smooth(aes(y=((1-exp(-log_est2))-(1-severity))/(1-severity),colour= "AR(1) error model"),se= F)+
  xlab("partial autocorrelation function")+
  ylab("false positives")+
  facet_wrap(severity~domain)+
  scale_colour_manual(values= c("grey","black"))+
  theme(legend.title=element_blank(),
        legend.position= c(0.2,0.85),
        panel.grid= element_blank(),
        plot.background= element_blank(),
        panel.background= element_rect(colour= "black",fill= "white"),
        legend.key= element_blank(),
        strip.background=element_blank(),
        legend.background= element_blank())+
  scale_y_continuous(limits= c(0,0.6),expand= c(0,0))

plot_AR <- ggplot(aes(ctl_acf,ps2+pp2),data= plot_dat2)+
  geom_point(aes(y= ps+pp,colour= "AR(1) error model"))+
  geom_point(aes(y= ps2+pp2,colour="white noise error model"))+
  binomial_smooth(se=F,aes(colour="white noise error model",weight=n))+
  binomial_smooth(aes(y=ps+pp,colour= "AR(1) error model",weight=n),se= F)+
  xlab("partial autocorrelation function")+
  ylab("false impact detections")+
  facet_wrap(~domain)+
  scale_colour_manual(values= c("grey","black"))+
  theme(legend.title=element_blank(),
        legend.position= c(0.2,0.85),
        panel.grid= element_blank(),
        plot.background= element_blank(),
        panel.background= element_rect(colour= "black",fill= "white"),
        legend.key= element_blank(),
        strip.background=element_blank(),
        legend.background= element_blank())+
  scale_y_continuous(limits= c(0,0.6),expand= c(0,0))

#pdf(width= 6,height= 3.75,file="4_Figures/Figure_5.pdf")
plot_AR
#dev.off()

#### Type II error against control cv and pre-impact spatial variation
out_sum <- out_sum%>%
  mutate(n_prop= I((ctl_n+imp_n)/23))
sub <- subset(out_sum, incl==0.15&domain=="quart"&severity=="0.8")
sub2 <- subset(out_sum, incl==0.15&domain=="quart"&severity=="0")
fit <- glm((ps)~ctl_cv+n_prop+sv_pre+ctl_acf,data= sub%>%data.frame(),family = binomial(link= 'logit'), weights= n)
fit2 <- glm((ps)~ctl_cv+n_prop+sv_pre+ctl_acf,data= sub2%>%data.frame(),family = binomial(link= 'logit'), weights= n)

sub1a <- sub%>%
  data.frame()%>%
  mutate(n_prop= mean(n_prop,na.rm=T),
         sv_pre= mean(sv_pre,na.rm=T),
         group2= "CV",
         ctl_resp = ctl_cv)  
sub1a$pred <- predict(fit,newdata= sub1a,type= "response")

sub1b <- sub%>%
  data.frame()%>%
  mutate(ctl_cv= mean(ctl_cv,na.rm=T),
         sv_pre= mean(sv_pre,na.rm=T),
         group2= "n",
         ctl_resp = n_prop)  
sub1b$pred <- predict(fit,newdata= sub1b,type= "response")

sub1c <- sub%>%
  data.frame()%>%
  mutate(ctl_cv= mean(ctl_cv,na.rm=T),
         n_prop= mean(n_prop,na.rm=T),
         group2= "sv",
         ctl_resp = sv_pre)  
sub1c$pred <- predict(fit,newdata= sub1c,type= "response")

sub5 <- bind_rows(list(sub1a,sub1b,sub1c))

sub1a <- sub2%>%
  data.frame()%>%
  mutate(n_prop= mean(n_prop,na.rm=T),
         sv_pre= mean(sv_pre,na.rm=T),
         group2= "CV",
         ctl_resp = ctl_cv)  
sub1a$pred <- predict(fit2,newdata= sub1a,type= "response")

sub1b <- sub2%>%
  data.frame()%>%
  mutate(ctl_cv= mean(ctl_cv,na.rm=T),
         sv_pre= mean(sv_pre,na.rm=T),
         group2= "n",
         ctl_resp = n_prop)  
sub1b$pred <- predict(fit2,newdata= sub1b,type= "response")

sub1c <- sub2%>%
  data.frame()%>%
  mutate(ctl_cv= mean(ctl_cv,na.rm=T),
         n_prop= mean(n_prop,na.rm=T),
         group2= "sv",
         ctl_resp = sv_pre)  
sub1c$pred <- predict(fit2,newdata= sub1c,type= "response")


sub6 <- bind_rows(list(sub5,sub1a,sub1b,sub1c))

binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}

table_data <- subset(out_sum, incl==0.15&domain=="quart"&severity%in%c("0","0.8"))%>%
  select(species,ps,ctl_cv,n_prop,sv_pre,ctl_acf,severity)%>%
  left_join(codes,by= "species")

#write.csv(table_data,"~/Dropbox/DanO_DOI_partnership/DOI_Impacts/4_Figures/Figure_5_Table.csv",row.names= F)

sub6$group2 <- factor(sub6$group2,levels= c("n","CV","sv"))
levels(sub6$taxa_group)[3] <- "mobile invert"
p1 <- ggplot(aes(ctl_resp,y= ps),data=sub6)+
      binomial_smooth(aes(weight= n),colour= "black",se= F)+
      facet_grid(severity~group2,scales= "free")+
      geom_point(aes(y= ps,shape= taxa_group,fill= taxa_group))+
      ylab("Power to detect impact                  False impact detections")+
      xlab("Prop Sites Sufficient            Mean Control CV               Spatial Synchrony")+
      theme_bw()+
      theme(panel.grid= element_blank(),
            strip.background=element_blank(),
            strip.text= element_blank(),
            legend.title= element_blank(),
            legend.position= c(0.15,0.9),
            legend.key.height = unit(0.07,"cm"),
            legend.background= element_blank(),
            plot.background= element_blank(),
            legend.key= element_blank())+
      scale_fill_manual(values= c("#1b9e77","#7570b3","black","#d95f02"))+
      scale_shape_manual(values= c(21,25,8,23));p1
    
#pdf(width= 6,height= 4,file="4_Figures/Figure_6.pdf")
p1
#dev.off()
    
sub <- subset(out_sum, incl==0.15&domain=="quart"&severity=="0.8")
fit <- glm((ps)~ctl_cv+sv_pre+ctl_acf,data= sub%>%data.frame(),family = binomial(link= 'logit'), weights= n)

sub2 <- sub%>%
  data.frame()%>%
  mutate(ctl_cv = mean(ctl_cv,na.rm=T),
         ctl_acf= mean(ctl_acf,na.rm=T),
         group2= "SV",
         ctl_resp =sv_pre)

sub2$pred <- predict(fit,newdata= sub2,type= "response")

sub3 <- sub%>%
  data.frame()%>%
  mutate(sv_pre= mean(sv_pre,na.rm=T),
         ctl_acf= mean(ctl_acf,na.rm=T),
         group2= "CV",
         ctl_resp = ctl_cv)  
sub3$pred <- predict(fit,newdata= sub3,type= "response")

sub4 <- sub%>%
  data.frame()%>%
  mutate(sv_pre= mean(sv_pre,na.rm=T),
         ctl_cv= mean(ctl_cv,na.rm=T),
         group2= "ACF",
         ctl_resp = ctl_acf)  
sub4$pred <- predict(fit,newdata= sub4,type= "response")

sub5 <- rbind_all(list(sub2,sub3,sub4))

p2 <- ggplot(aes(ctl_resp,ps+pp),data=sub5)+
      binomial_smooth(aes(weight= n),colour= "black")+
      facet_wrap(~group2,scales= "free_x")+
      geom_point(aes(y= (ps+pp),shape= group,fill= group))+
      ylab("Type II Error")+
      xlab("Mean Pre-Impact Autocorrelation          Mean Control Coefficient of Variation          Mean Pre-Impact    Spatial Variation")+
      theme_bw()+
      theme(panel.grid= element_blank(),
            strip.background=element_blank(),
            strip.text= element_blank(),
            panel.background= element_rect(fill= NA),
            legend.title= element_blank(),
            legend.position= "top")+
      #scale_y_continuous(limits= c(0,0.25))+
      scale_fill_manual(values= c("grey60","black"))+
      scale_shape_manual(values= c(21,24));p2

#### type II vs type I error rates

out2 <- read.csv("Output/output_summary_full.csv")

a <-  out2%>%
  filter(severity%in%c(0.8)&domain=="quart")%>%
  data.frame()%>%
  mutate(tp=(sig_negative))%>%
  select(severity,species,tp)

b <-  out2%>%
  filter(severity%in%c(0)&domain=="quart")%>%
  data.frame()%>%
  mutate(tn=1-(sig_negative))%>%
  select(species,tn)

test <- join(a,b)%>%
  mutate(Code=species)%>%
  join(inf_rank)%>%
  mutate(Inf2 = tp/(tp+1-tn)+tn/(tn+1-tp)-1)

plot(Inf2~Informedness,data= test)

Pwr_vs_FP <- ggplot(data= test,aes(1-tn,tp))+
  geom_point(aes(colour=Inf2))+
  geom_text(data= subset(test,Code%in%c(1,27)),aes(x= (1-tn)+0.015,y= tp+0.01,label =c("Macrocystis pyrifera","Embiotoca jacksoni")),fontface= "italic",hjust=.5)+
  geom_segment(data= subset(test,Code%in%c(1,27)),aes(xend= (1-tn)+0.01,yend= tp-0.01))+
  xlab("False Positives")+
  ylab("Power")+
  coord_flip()+
  theme_bw()+
  scale_x_continuous(limits= c(0,0.2),expand= c(0,0))+
  scale_y_continuous(limits= c(0,01),expand= c(0,0))+
  geom_vline(xintercept= 0.05,linetype= "dotted")+
  guides(colour=guide_colourbar(title.position = "top"))+
  theme(panel.grid= element_blank(),
        legend.direction= "horizontal",
        legend.position= c(0.3,0.8),
        legend.key.width= unit(0.8,"cm"))+
  scale_colour_gradientn(colors= rainbow(20),name= "Informedness");Pwr_vs_FP

#pdf(width= 5,height= 4.5,file="4_Figures/Figure_7.pdf")
Pwr_vs_FP
#dev.off()

### Informedness Plot ###
inf_data <- read.csv("Output/cum_informedness.csv")

head(inf_data)

inf_plot <- ggplot(aes(rank,YJ,colour= order),data= inf_data)+
  geom_point()+
  geom_path()+
  scale_colour_manual(values=c("blue","red"),name= "Individual Informedness")+
  stat_smooth(se=F)+
  ylab("Informedness (Youden's J)")+
  xlab("Number of taxa analyzed")+
  theme(axis.text.y= element_text(hjust=1),
        legend.background = element_blank(),
        legend.key.height= unit(0.2,"cm"),
        legend.key.width= unit(0.4,"cm"),
        legend.key=element_blank(),
        legend.position= c(0.7,0.3),
        strip.background=element_blank(),
        panel.background= element_blank(),
        plot.background= element_blank(),
        panel.border= element_rect(colour= "black",fill=NA));inf_plot

#pdf(width= 5,height= 5,file="4_Figures/Figure_8.pdf")
#inf_plot
#dev.off()

