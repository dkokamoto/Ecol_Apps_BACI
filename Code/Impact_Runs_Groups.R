###########################################################
###  Code for Rassweiler et al. (2021) Ecol. Appl. 
###  Authors: D.K. Okamoto, A. Rassweiler
###  Last Updated: May, 2018 (commented Nov. 2020)  
###  
###  
###  **CAUTION: code is not designed for efficiency or for
###  other applications, but for this paper alone!**
###########################################################

###  Code provided here is for *grouped* analysis (e.g. by taxonomic
###  group or by individual informedness)

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

### read in ecological data and define groups
### groups can be defined as taxonomic groupings (as below) or as "informedness" rank
allData = read.csv("Data/Species_Data.csv")%>%
  mutate(density=count/totalArea,
         group=ifelse(NEW_CODE%in%c(1,4,9,15),1,
                      ifelse(NEW_CODE%in%c(2,3,5,10,11,13),2,
                             ifelse(NEW_CODE%in%c(8,12,14),3,
                                    ifelse(NEW_CODE%in%c(25:39),4,5)))))

### species of interest
group_list= c(1:4)

### min proportion of years with positive counts in each site
min_props= c(0.15)

### number of samples per species
sample_size = 1000

### severities to consider
severities = c(0.8,0.5,0)

### cores to use
n.cores= 2

### dataframe of species, severity, and index
hit_df <- expand.grid(group= group_list,
                      severity= severities,
                      min_prop= min_props)

### replicate rows of the dataframe by sample size
hit_df <- hit_df[rep(1:nrow(hit_df),each= sample_size),]

pred_data <- data.frame(totalArea= 1,density= 1,count= 1,postImpact= factor(c(0,1)),
                        Impact=factor(c(1),levels= c(0,1)),BACI= factor(c(0,1)),new_count= 1,INDID=1)

### read in hitList for impacted sites
files = paste0(list.files(pattern= "domain_quart",path= "Data"))

### test function
hit_list = read.csv("Data/domain_quart.csv")[,-3]	
hit_df$index <- as.vector(replicate(nrow(hit_df)/sample_size,sample(unique(hit_list$indexNo),sample_size)))
id_list <- rows.to.list(hit_df)

### profile function for efficiency
Rprof(tmp <- tempfile())
system.time(test <- hit_fun_groups(id_list[[10000]],plot=FALSE))
Rprof()
summaryRprof(tmp)
test


### run if you don't want to overwrite existing files because of an interruption
#read_file <- function(x) {
#  a<- gsub("hit_analysis_|_domain_|R[[:digit:]]+|[a-z]+|_|.csv.rds","",x)
#  return(a)
#}

### create groups for saving snippets of 400 at a time as the simulation runs
### we are only focused on impacting a quarter of the sites for these groups
groups <- split(1:length(id_list), cut_number(1:length(id_list), n=400))
ranges <- sapply(groups,function(x) paste(range(x),collapse="-"))
file_list <- list.files(pattern="hit_analysis_group",path= "Output")
file_nums <- sapply(file_list,read_file)
groups <- groups[!(ranges%in%file_nums)]

### create a wrapper to run the function and save the collection 
ap_fun <- function(x){
  out <- ldply(id_list[x],hit_fun_groups)
  saveRDS(out,file= paste0("Output/hit_analysis_group_",paste(range(x),collapse="-"),".rds"))
  return(out)
}

### run the function in parallel - will save each snippet as it runs in case the process hangs
out <- mclapply(groups,ap_fun,mc.cores= n.cores,mc.silent = TRUE,mc.preschedule= FALSE)

### read the snippets and compile
file_list <- list.files(pattern="\\hit_analysis_group_",path= "Output")
read_RDS <- function(x) {
  runs <- readRDS(paste0("Output/",x))
  runs$domain <- gsub("hit_analysis_group_[[:digit:]]+-[[:digit:]]+||_domain_|R|-|.csv.rds","",x)
  return(runs)
}
runs_full <- ldply(file_list,read_RDS)