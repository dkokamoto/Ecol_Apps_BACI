
### convert rows to list
rows.to.list <- function( df ) {
  ll<-apply(df,1,list)
  ll<-lapply(ll,unlist)
}

### simulate a single Poisson random variable unless lambda = 0
zeropois <- function(lambda) ifelse(lambda==0,0,rpois(1,lambda))

### warnings 
withWarnings <- function (expr) 
{ 
  warnings <- character() 
  retval <- withCallingHandlers(expr, warning = function(ex) { 
    warnings <<- c(warnings, conditionMessage(ex)) 
    invokeRestart("muffleWarning") 
  }) 
  list(Value = retval, Warnings = warnings) 
} 

### Impact and statistical model functions ---
### Returns output with the following variables: 
# coef2_p = p-value associated with BACI coefficient w/ AR(1)
# coef3_p = p-value associated with BACI coefficient w/o AR(1)
# leffect2 = BACI effect size (AR(1))
# leffect2_se = se of BACI effect size (AR(1))
# leffect2 = BACI effect size (w/o AR(1))
# leffect2_se = se of BACI effect size (w/o AR(1))
# baci_p2 =  p-value associated with BACI LR test (AR(1))
# baci_p3 = p-value associated with BACI LR test (w/o AR(1))
# species= species ID
# severity = severity of the impact imposed
# index = index ID
# incl= minimum prorportion of sites with positive values
# ctl_mu= mean control value
# imp_mu= mean impact value
# ctl_acf= mean acf of the control sites
# imp_acf= mean acf of the impact sites
# ctl_cv= mean CV of the control sites
# imp_cv= mean CV of the impact sites
# ctl_sv_pre= spatial variation of the control sites pre-impact
# imp_sv_pre= spatial variation of the impact sites pre-impact
# ctl_sv_post= spatial variation of the control sites post-impact
# imp_sv_post= spatial variation of the impact sites post-impact
# sv_pre = total spatial variation pre-impact
# ctl_n = number of control sites with sufficient data
# imp_n = number of impact sites with sufficient data
# ctl_ny = mean pre-impact time series length of control sites 
# ctlp_ny = mean post-impact time series length of control sites 
# imp_ny = mean pre-impact time series length of impact sites 
# impp_ny = mean post-impact time series length of impact sites 
# IY= impact year
# err = any resulting errors
# warn2 = warnings for the model with AR(1)
# warn3 = warnings for the model w/o AR(1) model


### ungrouped function 
hit_fun <- function(x, power= FALSE,plot= FALSE,pred_data= pred_data){
  id_data <- data.frame(t(x))
  IY= sample(1990:2005,1)
  w.list=NA
  e.list= NA
  
  HL = hit_list%>%
    dplyr::filter(indexNo==id_data$index)%>%
    dplyr::select(site)
  
  newData <- allData%>%
    filter(NEW_CODE==id_data$species&YEAR>=(IY-10)&YEAR<=(IY+5))%>%
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

### grouped function
hit_fun_groups <- function(x, power= FALSE,plot= FALSE,pred_data= pred_data){
  x =id_list[[10000]]
  id_data <- data.frame(t(x))
  IY= sample(1990:2005,1)
  w.list=NA
  e.list= NA
  
  HL = hit_list%>%
    dplyr::filter(indexNo==id_data$index)%>%
    dplyr::select(site)
  
  newData <- allData%>%
    filter(group%in%id_data$group&YEAR>=(IY-10)&YEAR<=(IY+5))%>%
    mutate(postImpact = factor(ifelse(YEAR>=IY,1,0),levels= c(0,1)),
           Impact = factor(ifelse(SITE%in%HL$site,1,0),levels= c(0,1)),
           BACI = factor(ifelse(SITE%in%HL$site&YEAR>=IY,1,0),levels= c(0,1)))
  
  predData <- data.frame(totalArea= 1,density= 1,count= 1,post_impact= 1,Impact=c(0,1),BACI= c(0,1),new_count= 1,INDID=1)
  
  newData <- newData%>%
    group_by(SITE,NEW_CODE)%>%
    dplyr::summarize(lgz = mean.default(count>0))%>%
    data.frame()%>%
    join(newData)%>%
    filter(lgz>id_data$min_prop)
  
  
  ctl_n=length(unique(subset(newData,Impact==0)$SITE))
  imp_n=length(unique(subset(newData,Impact==1)$SITE))
  
  if(length(unique(newData$Impact))==2&length(unique(newData$postImpact))==2){ 
    
    mean_years = newData%>%
      group_by(Impact,postImpact,SITE)%>%
      dplyr::summarize(n_year=length(unique(YEAR)))%>%
      group_by(Impact,postImpact)%>%
      dplyr::summarize(mean_years = mean(n_year,na.rm=T))%>%
      complete(Impact,postImpact)
    
    newprob <- subset(newData,postImpact==1&Impact==1)$count*(1-id_data$severity)
    impdata <- floor(newprob)+sapply(newprob - floor(newprob ),function(x) rbinom(1,1,x))
    
    newData <-  newData%>%
      mutate(new_count=ifelse(postImpact==1&Impact==1,impdata,count))%>%
      data.frame()
    
    newData$INDID <- 1:nrow(newData)
    newData$time <- factor(newData$YEAR)
    newData$SITE <- factor(newData$SITE)
    newData$NEW_CODE <- factor(newData$NEW_CODE)
    
    fit3 <- tryCatch(
      withWarnings(glmmTMB(new_count~postImpact*Impact+
                             (0+postImpact*Impact|NEW_CODE)+(1|NEW_CODE:SITE)+ar1(time+0|NEW_CODE%in%SITE)+(1|INDID),
                           offset= log(totalArea),
                           family= poisson(link= "log"),
                           data= newData)),
      error = function(error) {return(list(Value= NULL,Warnings= error))})
    
    fit4 <- tryCatch(
      withWarnings(glmmTMB(new_count~  postImpact*Impact+ (0+postImpact*Impact|NEW_CODE)+(1|NEW_CODE:SITE)+(1|INDID),
                           offset= log(totalArea),
                           family= poisson(link= "log"),
                           data= newData)),
      error = function(error) {return(list(Value= NULL,Warnings= error))})
    if(!(is.null(fit3))){
      cf2 <- withWarnings(summary(fit3$Value)$coefficients$cond)
      if(!is.null(fit4)){
        cf3 <- withWarnings(summary(fit4$Value)$coefficients$cond)
        leffect3 <- log(exp(sum(cf3$Value[,1]))/exp(sum(cf3$Value[c(1,3),1])))
        se_leffect3 <- deltamethod(~log(exp(x1+x2+x3+x4)/(exp(x1+x3))),cf3$Value[,1],vcov(fit4$Value)$cond)
      } else {
        cf3 <- NA
        leffect3 <- NA
        se_leffect3 <- NA
      }
      
      leffect2 <- log(exp(sum(cf2$Value[,1]))/exp(sum(cf2$Value[c(1,3),1])))
      se_leffect2 <- deltamethod(~log(exp(x1+x2+x3+x4)/(exp(x1+x3))),cf2$Value[,1],vcov(fit3$Value)$cond)
      
      df <- data.frame(list(coef2_p=signif(cf2$Value[4,4],5),
                            coef3_p=signif(cf3$Value[4,4],5),
                            leffect2=signif(leffect2,5),
                            leffect2_se=signif(se_leffect2,5),
                            leffect3=signif(leffect3,5),
                            leffect3_se=signif(se_leffect3,5),
                            group= id_data$group,
                            severity = id_data$severity,
                            index = id_data$index,
                            incl= id_data$min_prop,
                            ctl_n=ctl_n,
                            imp_n=imp_n,
                            ctl_ny = mean_years$mean_years[1],
                            ctlp_ny = mean_years$mean_years[2],
                            imp_ny = mean_years$mean_years[3],
                            impp_ny = mean_years$mean_years[4],
                            IY= IY,
                            err = NA,
                            warn2 = ifelse(length(fit3$Warnings)>0,paste0(fit3$Warnings,sep=" "),NaN),
                            warn3 = ifelse(length(fit4$Warnings)>0,paste0(fit4$Warnings,sep=" "),NaN)))
    } else {
      df <- data.frame(list(   coef2_p=NA,
                               coef3_p=NA,
                               leffect2=NA,
                               leffect2_se=NA,
                               leffect3=NA,
                               leffect3_se=NA,
                               species= id_data$species,
                               severity = id_data$severity,
                               index = id_data$index,
                               incl= id_data$min_prop,
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
                            leffect2=NA,
                            leffect2_se=NA,
                            leffect3=NA,
                            leffect3_se=NA,
                            species= id_data$species,
                            severity = id_data$severity,
                            index = id_data$index,
                            incl= id_data$min_prop,
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
  return(cbind(df,id_data))
}
