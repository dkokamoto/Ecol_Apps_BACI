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
source("Code/Functions.R")

########################################
### scenarios to run 
### sets of impacts: 0,50,80
### all species
### quarter hit, quarter hit shuffle 
### test cases - quarter hit shuffle 
### with zero severity - FALSE POSITIVES
########################################


#read in ecological data
allData = read.csv("1_Data/Species_Data.csv")%>%
  mutate(density=count/totalArea)

### species of interest
species_list= c(1:12,16:24,25:39)

### min proportion of years with positive counts in each site
min_props= c(0.15,0.5)

### number of samples per species
sample_size = 1000

### severities to consider
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
hit_df <- expand.grid(species= species_list,
                      severity= severities,
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
#file_list <- list.files(pattern="\\.rds",path= "3_Output")
file_list= NULL
file_nums <- sapply(file_list,read_file)

### 
for(i in 1:length(files)){
  hit_list = read.csv(paste0("1_Data/",files[i]))[,-3]	
  hit_df$index <- as.vector(replicate(nrow(hit_df)/sample_size,sample(1:9999,sample_size)))
  id_list <- rows.to.list(hit_df)
  groups <- split(1:length(id_list), cut_number(1:length(id_list), n=400))
  ranges <- sapply(groups,function(x) paste(range(x),collapse="-"))
  groups <- groups[!(ranges%in%file_nums)]
  ap_fun <- function(x){
    out <- ldply(id_list[x],hit_fun)
    saveRDS(out,file= paste0("3_Output/hit_analysis_",paste(range(x),collapse="-"),"_",files[i],".rds"))
    return(out)
  }
  out <- mclapply(groups,ap_fun,mc.cores= n.cores,mc.silent = TRUE)
}
#

file_list <- list.files(pattern="\\.rds",path= "3_Output")
read_RDS <- function(x) {
  a <- readRDS(paste0("3_Output/",x))
  a$domain <- gsub("hit_analysis_[[:digit:]]+-[[:digit:]]+||_domain_|R|-|.csv.rds","",x)
  return(a)
}

file_list <- list.files(pattern="\\.rds",path= "3_Output")

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