###########################################################
###  Code for Rassweiler et al. (2021) Ecol. Appl. 
###  Authors: D.K. Okamoto, A. Rassweiler
###  Last Updated: May, 2018 (commented Nov. 2020)  
###  
###  
###  **CAUTION: code is not designed for efficiency or for
###  other applications, but for this paper alone!**
###########################################################

###  Code provided here is for *ungrouped* analysis (e.g. a single taxa)

### load necessary packages

library("ggplot2")
library("VGAM")
library("car")
library("DHARMa")
library("lme4")
library("plyr")
library("dplyr")
library("lme4")
library("reshape2")
library("boot")
library("parallel")
library("msm")
library("tidyr")
library("glmmTMB")

### load necessary functions
source("2_Code/Functions.R")

########################################
### scenarios to run 
### sets of impacts: 0,50,80
### all species
### quarter hit, quarter hit shuffle 
### test cases - quarter hit shuffle 
### with zero severity - FALSE POSITIVES
########################################

#read in ecological data
allData = read.csv("Data/Species_Data.csv")%>%
  mutate(density=count/totalArea)

### species of interest
species_list= c(1:12,16:24,25:39)

### min proportion of years with positive counts in each site
min_props= c(0.15,0.5)

### number of samples per species
sample_size = 1000

### severities to consider (a subset for some, a grid for others)
severities = c(0.8,0.5,0,-.5,-1)
# severities = c(rev(seq(0,0.75,by= 0.075)),1/(1-seq(.075,0.75,by= 0.075))-1)

### cores to use
n.cores= 20

### KEY OUTPUT: 
### significance of impact 
### effect size for impact 
### diagnostics (overdispersion/zero-inflation fit)
### other error codes
### no impact or control sites with observations
### power given sample size (?)

### dataframe of species, severity, and index
hit_df <- expand.grid(severity= severities,
                      min_prop= min_props)

hit_df <- hit_df[rep(1:nrow(hit_df),each= sample_size),]

rows.to.list <- function( df ) {
  ll<-apply(df,1,list)
  ll<-lapply(ll,unlist)
}

zeropois <- function(lambda) ifelse(lambda==0,0,rpois(1,lambda))
#Species codes we need to evaluate 
#1-12, 16-24 (point contact), 25-39

## warnings 
withWarnings <- function (expr) 
{ 
  warnings <- character() 
  retval <- withCallingHandlers(expr, warning = function(ex) { 
    warnings <<- c(warnings, conditionMessage(ex)) 
    invokeRestart("muffleWarning") 
  }) 
  list(Value = retval, Warnings = warnings) 
} 

pred_data <- data.frame(totalArea= 1,density= 1,count= 1,postImpact= factor(c(0,1)),Impact=factor(c(1),levels= c(0,1)),BACI= factor(c(0,1)),new_count= 1,INDID=1)

hit_fun <- function(x, power= FALSE,plot= FALSE,pred_data= pred_data){
  x <- id_list[[1]]
  id_data <- data.frame(t(x))
  IY= sample(1990:2005,1)
  w.list=NA
  e.list= NA
  
  HL = hit_list%>%
    dplyr::filter(indexNo==id_data$index)%>%
    dplyr::select(site)
  
  newData <- allData%>%
    filter(NEW_CODE%in%id_data$species&YEAR>=(IY-10)&YEAR<=(IY+5))%>%
    mutate(postImpact = factor(ifelse(YEAR>=IY,1,0),levels= c(0,1)),
           Impact = factor(ifelse(SITE%in%HL$site,1,0),levels= c(0,1)),
           BACI = factor(ifelse(SITE%in%HL$site&YEAR>=IY,1,0),levels= c(0,1)))
  
  predData <- data.frame(totalArea= 1,density= 1,count= 1,post_impact= 1,Impact=c(0,1),BACI= c(0,1),new_count= 1,INDID=1)
  
  sum_data <- newData%>%
    group_by(SITE)%>%
    dplyr::summarize(lgz = mean.default(count>0))%>%
    data.frame()
  
  newData <- subset(newData,SITE%in%subset(sum_data,lgz>id_data$min_prop)$SITE)
  
  sync <- function(y){sum(y)/sum(sqrt(diag(y)))^2}
  all_sync = data.frame(expand.grid(Impact=c(1,0),postImpact= c(1,0),sync=NA))
  for(i in 1:length(unique(factor(newData$Impact):factor(newData$postImpact)))){
    all_sync[factor(all_sync$Impact):factor(all_sync$postImpact)==unique(factor(newData$Impact):factor(newData$postImpact))[i],"sync"] <- newData%>%
      filter(factor(newData$Impact):factor(newData$postImpact)==unique(factor(newData$Impact):factor(newData$postImpact))[i])%>%
      dplyr::select(SITE,YEAR,density)%>%
      dcast(YEAR~SITE,value.var= "density")%>%
      dplyr::select(-one_of("YEAR"))%>%
      filter(complete.cases(.))%>%
      var(.)%>%
      sync(.)
  }
  
  pre_sync <- newData%>%
    filter(newData$postImpact==0)%>%
    dplyr::select(SITE,YEAR,density)%>%
    dcast(YEAR~SITE,value.var= "density")%>%
    dplyr::select(-one_of("YEAR"))%>%
    filter(complete.cases(.))%>%
    var(.)%>%
    sync(.)
  
  
  mean_acf <- newData%>%
    group_by(SITE,Impact)%>%
    dplyr::summarize(acf= pacf(density,plot= F,lag.max=1)[[1]][[1]],
                     mu= mean.default(density),
                     site_sd= sd(density,na.rm=T))%>%
    group_by(Impact)%>%
    dplyr::summarize(macf=mean.default(acf,na.rm=TRUE),
                     cv= mean.default(site_sd)/mean.default(mu,na.rm=TRUE),
                     mu= mean.default(mu,na.rm=TRUE))%>%
    dplyr::select(Impact,macf,cv,mu)%>%
    melt()
  
  ctl_n=length(unique(subset(newData,Impact==0)$SITE))
  imp_n=length(unique(subset(newData,Impact==1)$SITE))
  
  mean_years = newData%>%
    group_by(Impact,postImpact,SITE)%>%
    dplyr::summarize(n_year=length(unique(YEAR)))%>%
    group_by(Impact,postImpact)%>%
    dplyr::summarize(mean_years = mean(n_year,na.rm=T))%>%
    complete(Impact,postImpact)
  
  if(length(unique(newData$Impact))==2&length(unique(newData$postImpact))==2){ 
    newprob <- subset(newData,postImpact==1&Impact==1)$count*(1-id_data$severity)
    impdata <- floor(newprob)+sapply(newprob - floor(newprob ),function(x) rbinom(1,1,x))
    
    newData <-  newData%>%
      mutate(new_count=ifelse(postImpact==1&Impact==1,impdata,count))%>%
      data.frame()
    
    newData$INDID <- 1:nrow(newData)
    newData$time <- factor(newData$YEAR)
    newData$SITE <- factor(newData$SITE)
    newData <<- newData
    if(id_data$species %in% c(16:24)){
      ### fit a glmm with overdispersion
      # fit <- tryCatch(
      #   withWarnings(glmer(new_count/totalArea~postImpact*Impact+
      #                                             (1|SITE)+(1|INDID),
      #                                             family= binomial(link= "logit"),
      #                                             data= newData,weights= totalArea,
      #                           control=glmerControl(optimizer="bobyqa",
      #                                                optCtrl=list(maxfun=1e5)))),
      #   error = function(error) {return(list(Value= NULL,Warnings= error))})
      
      fit3 <- tryCatch(
        withWarnings(glmmTMB(new_count/totalArea~postImpact*Impact+
                               (1|SITE)+ar1(time+0|SITE)+(1|INDID),
                             weights= totalArea,
                             family= binomial(link= "logit"),
                             data= newData)),
        error = function(error) {return(list(Value= NULL,Warnings= error))})
      
      fit4 <- tryCatch(
        withWarnings(glmmTMB(new_count/totalArea~postImpact*Impact+
                               (1|SITE)+(1|INDID),
                             weights= totalArea,
                             family= binomial(link= "logit"),
                             data= newData)),
        error = function(error) {return(list(Value= NULL,Warnings= error))})
    } else {
      # fit <-   tryCatch(
      #   withWarnings(glmer.nb(new_count~postImpact*Impact+(1|SITE),
      #                                         offset= log(totalArea),
      #                                          data= newData,
      #                                          control=glmerControl(optimizer="bobyqa",
      #                                          optCtrl=list(maxfun=1e5)))),
      #   error = function(error) {return(list(Value= NULL,Warnings= error))})
      
      fit3 <- tryCatch(
        withWarnings(glmmTMB(new_count~postImpact*Impact+
                               (1|SITE)+ar1(time+0|SITE)+(1|INDID),
                             offset= log(totalArea),
                             family= poisson(link= "log"),
                             data= newData)),
        error = function(error) {return(list(Value= NULL,Warnings= error))})
      
      fit4 <- tryCatch(
        withWarnings(glmmTMB(new_count~postImpact*Impact+
                               (1|SITE)+(1|INDID),
                             offset= log(totalArea),
                             family= poisson(link= "log"),
                             data= newData)),
        error = function(error) {return(list(Value= NULL,Warnings= error))})
    }
    if(!(is.null(fit3))){
      # fit2 <- tryCatch(
      #   withWarnings(update(fit$Value,.~.-postImpact:Impact)),
      # error = function(error) {return(list(Value= NULL,Warnings= error))})
      fit3.b <- tryCatch(
        withWarnings(update(fit3$Value,.~.-postImpact:Impact)),
        error = function(error) {return(list(Value= NULL,Warnings= error))})
      fit4.b <- tryCatch(
        withWarnings(update(fit4$Value,.~.-postImpact:Impact)),
        error = function(error) {return(list(Value= NULL,Warnings= error))})
      
      if(!(is.null(fit3.b$Value))){
        
        # ### look at residuals and diagnostics (see how it improves things!)
        # simulationOutput <- simulateResiduals(fittedModel = fit$Value, n = 250)
        # #plotSimulatedResiduals(simulationOutput = simulationOutput)
        # 
        # ### look at diagnostics for assumption of uniform residuals and assess zero-inflation
        # od <- testUniformity(simulationOutput = simulationOutput)$p.value
        # zi <- testZeroInflation(simulationOutput = simulationOutput,plot=F)$p.value
        # cf <- withWarnings(summary(fit$Value)$coefficients)
        cf2 <- withWarnings(summary(fit3$Value)$coefficients$cond)
        cf3 <- withWarnings(summary(fit4$Value)$coefficients$cond)
        
        # leffect <- ifelse(id_data$species %in% c(16:24),
        #                  log(inv.logit(sum(cf$Value[,1]))/inv.logit(sum(cf$Value[c(1,3),1]))),
        #                  log(exp(sum(cf$Value[,1]))/exp(sum(cf$Value[c(1,3),1]))))
        # 
        # se_leffect <- ifelse(id_data$species %in% c(16:24),
        #                     deltamethod(~log(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))/((exp(x1+x3))/(1+exp(x1+x3)))),cf$Value[,1],vcov(fit$Value)),
        #                     deltamethod(~log(exp(x1+x2+x3+x4)/(exp(x1+x3))),cf$Value[,1],vcov(fit$Value)))
        
        leffect2 <- ifelse(id_data$species %in% c(16:24),
                           log(inv.logit(sum(cf2$Value[,1]))/inv.logit(sum(cf2$Value[c(1,3),1]))),
                           log(exp(sum(cf2$Value[,1]))/exp(sum(cf2$Value[c(1,3),1]))))
        
        se_leffect2 <- ifelse(id_data$species %in% c(16:24),
                              deltamethod(~log(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))/((exp(x1+x3))/(1+exp(x1+x3)))),cf2$Value[,1],vcov(fit3$Value)$cond),
                              deltamethod(~log(exp(x1+x2+x3+x4)/(exp(x1+x3))),cf2$Value[,1],vcov(fit3$Value)$cond))
        
        leffect3 <- ifelse(id_data$species %in% c(16:24),
                           log(inv.logit(sum(cf3$Value[,1]))/inv.logit(sum(cf3$Value[c(1,3),1]))),
                           log(exp(sum(cf3$Value[,1]))/exp(sum(cf3$Value[c(1,3),1]))))
        
        se_leffect3 <- ifelse(id_data$species %in% c(16:24),
                              deltamethod(~log(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))/((exp(x1+x3))/(1+exp(x1+x3)))),cf3$Value[,1],vcov(fit4$Value)$cond),
                              deltamethod(~log(exp(x1+x2+x3+x4)/(exp(x1+x3))),cf3$Value[,1],vcov(fit4$Value)$cond))
        
        # baci <- anova(fit$Value,fit2$Value)$'Pr(>Chisq)'[2]
        baci2 <- anova(fit3$Value,fit3.b$Value)$'Pr(>Chisq)'[2]
        baci3 <- anova(fit4$Value,fit4.b$Value)$'Pr(>Chisq)'[2]
        
        df <- data.frame(list(#od_p=signif(od,5),
          #                       zi_p=signif(zi,5),
          #                       coef=signif(cf$Value[4,1],5),
          #                       coef_se=signif(cf$Value[4,2],5),
          #                       coef_p=signif(cf$Value[4,4],5),
          coef2_p=signif(cf2$Value[4,4],5),
          coef3_p=signif(cf3$Value[4,4],5),
          # leffect=signif(leffect,5),
          # leffect_se=signif(se_leffect,5),
          leffect2=signif(leffect2,5),
          leffect2_se=signif(se_leffect2,5),
          leffect3=signif(leffect3,5),
          leffect3_se=signif(se_leffect3,5),
          # baci_p = baci,
          baci_p2 = baci2,
          baci_p3 = baci3,
          species= id_data$species,
          severity = id_data$severity,
          index = id_data$index,
          incl= id_data$min_prop,
          ctl_mu= subset(mean_acf,Impact==0&variable== "mu")$value,
          imp_mu= subset(mean_acf,Impact==1&variable== "mu")$value,
          ctl_acf= subset(mean_acf,Impact==0&variable== "macf")$value,
          imp_acf= subset(mean_acf,Impact==1&variable== "macf")$value,
          ctl_cv= subset(mean_acf,Impact==0&variable== "cv")$value,
          imp_cv= subset(mean_acf,Impact==1&variable== "cv")$value,
          ctl_sv_pre= subset(all_sync,Impact==0&postImpact==0)$sync,
          imp_sv_pre= subset(all_sync,Impact==1&postImpact==0)$sync,
          ctl_sv_post= subset(all_sync,Impact==0&postImpact==1)$sync,
          imp_sv_post= subset(all_sync,Impact==1&postImpact==1)$sync,
          sv_pre = pre_sync,
          ctl_n=ctl_n,
          imp_n=imp_n,
          ctl_ny = mean_years$mean_years[1],
          ctlp_ny = mean_years$mean_years[2],
          imp_ny = mean_years$mean_years[3],
          impp_ny = mean_years$mean_years[4],
          IY= IY,
          err = NA,
          # warn1 = ifelse(length(fit$Warnings)>0,paste0(fit$Warnings,sep=" "),NaN),
          warn2 = ifelse(length(fit3$Warnings)>0,paste0(fit3$Warnings,sep=" "),NaN),
          warn3 = ifelse(length(fit4$Warnings)>0,paste0(fit4$Warnings,sep=" "),NaN)))
      } else {
        df <- data.frame(list(   coef2_p=NA,
                                 coef3_p=NA,
                                 # leffect=signif(leffect,5),
                                 # leffect_se=signif(se_leffect,5),
                                 leffect2=NA,
                                 leffect2_NA,
                                 leffect3=NA,
                                 leffect3_NA,
                                 # baci_p = baci,
                                 baci_p2 = NA,
                                 baci_p3 = NA,
                                 species= id_data$species,
                                 severity = id_data$severity,
                                 index = id_data$index,
                                 incl= id_data$min_prop,
                                 ctl_mu= subset(mean_acf,Impact==0&variable== "mu")$value,
                                 imp_mu= subset(mean_acf,Impact==1&variable== "mu")$value,
                                 ctl_acf= subset(mean_acf,Impact==0&variable== "macf")$value,
                                 imp_acf= subset(mean_acf,Impact==1&variable== "macf")$value,
                                 ctl_cv= subset(mean_acf,Impact==0&variable== "cv")$value,
                                 imp_cv= subset(mean_acf,Impact==1&variable== "cv")$value,
                                 ctl_sv_pre= subset(all_sync,Impact==0&postImpact==0)$sync,
                                 imp_sv_pre= subset(all_sync,Impact==1&postImpact==0)$sync,
                                 ctl_sv_post= subset(all_sync,Impact==0&postImpact==1)$sync,
                                 imp_sv_post= subset(all_sync,Impact==1&postImpact==1)$sync,
                                 sv_pre = pre_sync,
                                 ctl_n=ctl_n,
                                 imp_n=imp_n,
                                 ctl_ny = mean_years$mean_years[1],
                                 ctlp_ny = mean_years$mean_years[2],
                                 imp_ny = mean_years$mean_years[3],
                                 impp_ny = mean_years$mean_years[4],
                                 IY= IY,
                                 err = NA,
                                 warn = "no convergence"))
      }
    } else {
      df <- data.frame(list(  coef2_p=NA,
                              coef3_p=NA,
                              # leffect=signif(leffect,5),
                              # leffect_se=signif(se_leffect,5),
                              leffect2=NA,
                              leffect2_NA,
                              leffect3=NA,
                              leffect3_NA,
                              # baci_p = baci,
                              baci_p2 = NA,
                              baci_p3 = NA,
                              species= id_data$species,
                              severity = id_data$severity,
                              index = id_data$index,
                              incl= id_data$min_prop,
                              ctl_mu= subset(mean_acf,Impact==0&variable== "mu")$value,
                              imp_mu= subset(mean_acf,Impact==1&variable== "mu")$value,
                              ctl_acf= subset(mean_acf,Impact==0&variable== "macf")$value,
                              imp_acf= subset(mean_acf,Impact==1&variable== "macf")$value,
                              ctl_cv= subset(mean_acf,Impact==0&variable== "cv")$value,
                              imp_cv= subset(mean_acf,Impact==1&variable== "cv")$value,
                              ctl_sv_pre= subset(all_sync,Impact==0&postImpact==0)$sync,
                              imp_sv_pre= subset(all_sync,Impact==1&postImpact==0)$sync,
                              ctl_sv_post= subset(all_sync,Impact==0&postImpact==1)$sync,
                              imp_sv_post= subset(all_sync,Impact==1&postImpact==1)$sync,
                              sv_pre = pre_sync,
                              ctl_n=ctl_n,
                              imp_n=imp_n,
                              ctl_ny = mean_years$mean_years[1],
                              ctlp_ny = mean_years$mean_years[2],
                              imp_ny = mean_years$mean_years[3],
                              impp_ny = mean_years$mean_years[4],
                              IY= IY,
                              err= "insufficient data",
                              warn = "NULL"))
    }
  } else {
    df <- data.frame(list(  coef2_p=NA,
                            coef3_p=NA,
                            # leffect=signif(leffect,5),
                            # leffect_se=signif(se_leffect,5),
                            leffect2=NA,
                            leffect2_NA,
                            leffect3=NA,
                            leffect3_NA,
                            # baci_p = baci,
                            baci_p2 = NA,
                            baci_p3 = NA,
                            species= id_data$species,
                            severity = id_data$severity,
                            index = id_data$index,
                            incl= id_data$min_prop,
                            ctl_mu= NA,
                            imp_mu= NA,
                            ctl_acf= NA,
                            imp_acf= NA,
                            ctl_cv= NA,
                            imp_cv= NA,
                            ctl_sv= NA,
                            imp_sv= NA,
                            ctl_n=ctl_n,
                            imp_n=imp_n,
                            ctl_ny = mean_years$mean_years[1],
                            ctlp_ny = mean_years$mean_years[2],
                            imp_ny = mean_years$mean_years[3],
                            impp_ny = mean_years$mean_years[4],
                            IY= IY,
                            err= "insufficient data",
                            warn = "NULL"))
  }
  # df
  # ### plotting code for testing
  if(plot==TRUE){
    plot <- ggplot(aes(YEAR,new_count/totalArea),data= newData)+
      geom_ribbon(aes(ymin= -Inf,ymax= Inf),data= subset(newData,Impact==1&postImpact==1))+
      geom_line(aes(colour= factor(SITE)))+
      facet_wrap(~Impact)+
      theme_bw()+#geom_line(aes(y=pred/totalArea,colour= factor(SITE)))+
      ylab("proportional cover")+
      theme(legend.position= "none",
            strip.background = element_blank())
    return(list(df=df,plot= plot))
  } else {
    return(cbind(df,id_data))
  }
}

#read in hitList for impacted sites (there are several of these that we will use but this is probably the main one)
files = paste0(list.files(pattern= "domain",path= "Data"))[c(1:4)]

### test function
hit_list = read.csv(paste0("Data/",files[1]))[,-3]	
hit_df$index <- as.vector(replicate(nrow(hit_df)/sample_size,sample(unique(hit_list$indexNo),sample_size)))
id_list <- rows.to.list(hit_df)

### profile function for efficiency
Rprof(tmp <- tempfile())
a <- system.time(hit_fun(id_list[[10000]],plot=FALSE))
Rprof()
summaryRprof(tmp)

### run if you don't want to overwrite
read_file <- function(x) {
  a<- gsub("hit_analysis_|_domain_|R[[:digit:]]+|[a-z]+|_|.csv.rds","",x)
  return(a)
}
#file_list <- list.files(pattern="\\.rds",path= "Output")
file_list= NULL
file_nums <- sapply(file_list,read_file)

### 
for(i in 1:length(files)){
  hit_list = read.csv(paste0("Data/",files[i]))[,-3]	
  hit_df$index <- as.vector(replicate(nrow(hit_df)/sample_size,sample(1:9999,sample_size)))
  id_list <- rows.to.list(hit_df)
  groups <- split(1:length(id_list), cut_number(1:length(id_list), n=400))
  ranges <- sapply(groups,function(x) paste(range(x),collapse="-"))
  groups <- groups[!(ranges%in%file_nums)]
  ap_fun <- function(x){
    out <- ldply(id_list[x],hit_fun)
    saveRDS(out,file= paste0("Output/hit_analysis_",paste(range(x),collapse="-"),"_",files[i],".rds"))
    return(out)
  }
  out <- mclapply(groups,ap_fun,mc.cores= n.cores,mc.silent = TRUE)
}
#

file_list <- list.files(pattern="\\.rds",path= "Output")
read_RDS <- function(x) {
  a <- readRDS(paste0("Output/",x))
  a$domain <- gsub("hit_analysis_[[:digit:]]+-[[:digit:]]+||_domain_|R|-|.csv.rds","",x)
  return(a)
}

file_list <- list.files(pattern="\\.rds",path= "Output")

###
b <- ldply(file_list,read_RDS)
write.csv(b,"runs_output_10.10.17.csv",row.names= F)


b <- read.csv("runs_output_1")
# ggplot(aes(factor(species),abs(exp(leffect)-(1-severity))/(1-severity),
#            colour= factor(ifelse(species%in%c(16:24),1,0))),alpha= 0.5,data=b)+
#   geom_boxplot()+
#   facet_grid(~severity)+
#   geom_hline(yintercept= 0.4)+
#   coord_cartesian(ylim= c(0,1))+
#   theme_bw()
# 
# ggplot(aes(imp_cv,as.numeric(baci_xs_p<0.05)),data=b)+
#   geom_point()+
#   facet_grid(~severity)+
#   geom_hline(yintercept= 0.4)+
#   ylim(0,1)+
#   theme_bw()+
#   stat_smooth(method= "gam",
#                   method.args= list(family= "binomial"))