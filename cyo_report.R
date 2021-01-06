################################################################################
#```{ Set sample kind }
# sample.kind=NULL # if using R 3.5 or earlier
sample.kind='Rounding' # if using R 3.6 or later


#```{ Set to run or not run very heavy computations }
# Many computations were very heavy;
# thus, we provided the RDS files as subtitutes and load only ones
# that can be ran in most computers.
# Set to TRUE if you want to run the very heavy computations.
run_heavy_computation=FALSE






################################################################################
#```{ Install and set specific version of Bioconductor }
# Install devtools to install specific version of BiocManager
if(!require(devtools)) install.packages('devtools')

# Install specific version of BiocManager and Bioconductor
if(!require(BiocManager) | packageVersion('BiocManager')!='1.30.10'){
  devtools::install_version('BiocManager',version='1.30.10')
 }
if(version()!='3.11') install(version='3.11',update=TRUE,ask=FALSE)

# Unload devtools & BiocManager to prevent overlapped functions with others
detach('package:devtools',unload=T)
detach('package:BiocManager',unload=T)






################################################################################
#```{ Install specific packages with specific Bioconductor }
BiocManager::install(
  c('tidyverse','lubridate','dslabs','pbapply','broom','caret','igraph'
    ,'data.table','kableExtra','zeallot','doParallel','WCGNA','MLeval','glmnet'
    ,'Rborist','lda'
  ),
  update=F
)






################################################################################
#```{ Load packages }
library(tidyverse)
options(dplyr.summarise.inform=FALSE)
dslabs::ds_theme_set()
library(parallel)
library(pbapply)
library(lubridate)
library(broom)
library(caret)
library(igraph)
library(data.table)
library(kableExtra)
library(zeallot)
library(doParallel)
library(WGCNA)
library(MLeval)
library(gbm)






################################################################################
#```{ Load preselected data }
selection=readRDS('data/selection.rds')
target_population=readRDS('data/target_population.rds')






################################################################################
#```{ Create medical history by code encounter for each visit }
if(run_heavy_computation){
  
  medical_history=
    
    # Join all codes of diagnosis/procedure
    suppressWarnings(separate(
      readRDS('data/target_population.rds') %>%
        mutate(seq=seq(nrow(.))) %>%
        left_join(readRDS('data/admission_diagnosis.rds'),by='visit_id') %>%
        left_join(readRDS('data/secondary_diagnoses.rds'),by='visit_id') %>%
        left_join(readRDS('data/procedures.rds'),by='visit_id')
      ,icd9_code_desc,c('icd9_code','icd9_desc'),sep='\\s?-\\s?'
    )) %>%
    
    # Differentiate if a provider is primary care
    mutate(
      est_strata_id=
        est_strata_id %>%
        paste0('.',ifelse(str_detect(insurance_model,'Primary care'),1,2))
    ) %>%
    
    # Select diagnosis/procedure-related attributes only
    select(
      subject_id
      ,day_to_event
      ,est_strata_id
      ,icd10_code
      ,admission_icd10_code
      ,secondary_icd10_code
      ,icd9_code
      ,seq
    ) %>%
    
    # Pivot longer the codes
    pivot_longer(
      colnames(.) %>%
        .[!.%in%c('subject_id','day_to_event','est_strata_id','seq')]
      ,names_to='code_type'
      ,values_to='code'
    ) %>%
    
    # Give a sign if the pivoted code column is missing.
    # This is because not all visits have a particular type of code
    mutate(
      code=ifelse(is.na(code)|code=='','NA',code)
      ,value=ifelse(code=='NA',0,1)
    )
  
 }


    
    
    
    
################################################################################
#```{ Build a function to unevenly split a vector given a length }
split_len=function(x,len){
  split_idx=round(seq(1,length(x),len=len))
  lapply(X=1:(length(split_idx)-1),Y=split_idx,Z=x,function(X,Y,Z){
    Z[Y[X]:(Y[X+1]-ifelse(X==(length(split_idx)-1),0,1))]
   })
 }


      
      
      
      
################################################################################
#```{ Compute nationwide cumulated code encounter count of ... }
if(run_heavy_computation){
  cat('Compute nationwide cumulated encounter count of a code for',append=T)
  cat('each subject every visit\n')
  cat('Started:',as.character(now()),'\n')
  cl=makeCluster(detectCores()-1)
  clusterEvalQ(cl,{
    library('tidyverse')
    library('pbapply')
   })
    
    mh_nationwide=
      
      # Compute the number of days a code was encountered for a subject
      # up to each visit from any kind of codes and any healthcare providers
      medical_history %>%
      group_by(subject_id,day_to_event,code,seq) %>%
      summarize(
        d_enc_adm=sum(ifelse(code_type=='admission_icd10_code',value,0))
        ,d_enc_dur=sum(ifelse(code_type=='secondary_icd10_code',value,0))
        ,d_enc_dis=sum(ifelse(code_type=='icd10_code',value,0))
        ,p_enc=sum(ifelse(code_type=='icd9_code',value,0))
      ) %>%
      ungroup() %>%
      mutate(
        any_enc=as.integer(any(c(d_enc_adm>0,d_enc_dur>0,d_enc_dis>0,p_enc>0)))
      ) %>%
      select(-d_enc_adm,-d_enc_dur,-d_enc_dis,-p_enc) %>%
      spread(code,any_enc,fill=0) %>%
      select(-seq,-`NA`) %>%
      arrange(subject_id,desc(day_to_event))
    gc()
    
    mh_nationwide=
      
      # Compute cumulated number of the days for each code
      mh_nationwide %>% 
      .$subject_id %>%
      .[!duplicated(.)] %>%
      split_len(len=500) %>%
      pblapply(X=seq(length(.)),Y=.,Z=mh_nationwide,cl=cl,function(X,Y,Z){
        filter(Z,subject_id%in%Y[[X]]) %>%
          group_by(subject_id) %>%
          mutate_at(
            colnames(.) %>% .[!.%in%c('subject_id','day_to_event')]
            ,cumsum
          ) %>%
          ungroup()
       }) %>%
      do.call(rbind,.)

  stopCluster(cl)
  rm(cl)
  gc()
  cat('End:',as.character(now()))
  saveRDS(mh_nationwide,'data/mh_nationwide.rds')
 }else{
  cat(readRDS('data/log.rds')[['mh_nationwide']])
 }






################################################################################
#```{ Compute cumulated encounter ... by provider database }
if(run_heavy_computation){
  cat('Compute cumulated encounter count of a code for ',append=T)
  cat('each subject every visit isolated by provider database\n')
  cat('Started:',as.character(now()),'\n')
  cl=makeCluster(detectCores()-1)
  clusterEvalQ(cl,{
    library('tidyverse')
    library('pbapply')
   })
    
    mh_provider=
      
      # Compute the number of days a code was encountered for a subject
      # up to each visit from any kind of codes to each healthcare provider
      medical_history %>%
      group_by(subject_id,day_to_event,est_strata_id,code,seq) %>%
      summarize(
        d_enc_adm=sum(ifelse(code_type=='admission_icd10_code',value,0))
        ,d_enc_dur=sum(ifelse(code_type=='secondary_icd10_code',value,0))
        ,d_enc_dis=sum(ifelse(code_type=='icd10_code',value,0))
        ,p_enc=sum(ifelse(code_type=='icd9_code',value,0))
      ) %>%
      ungroup() %>%
      mutate(
        any_enc=as.integer(any(c(d_enc_adm>0,d_enc_dur>0,d_enc_dis>0,p_enc>0)))
      ) %>%
      select(-d_enc_adm,-d_enc_dur,-d_enc_dis,-p_enc) %>%
      spread(code,any_enc,fill=0) %>%
      select(-seq,-`NA`) %>%
      arrange(subject_id,desc(day_to_event),est_strata_id)
    gc()
    
    mh_provider=
      
      # Compute cumulated number of the days for each code
      mh_provider %>% 
      .$subject_id %>%
      .[!duplicated(.)] %>%
      split_len(len=500) %>%
      pblapply(X=seq(length(.)),Y=.,Z=mh_provider,cl=cl,function(X,Y,Z){
        filter(Z,subject_id%in%Y[[X]]) %>%
          group_by(subject_id,est_strata_id) %>%
          mutate_at(
            colnames(.) %>% .[!.%in%c('subject_id','day_to_event')]
            ,cumsum
          ) %>%
          ungroup()
       }) %>%
      do.call(rbind,.)

  stopCluster(cl)
  rm(cl)
  gc()
  cat('End:',as.character(now()))
  saveRDS(mh_provider,'data/mh_provider.rds')
 }else{
  cat(readRDS('data/log.rds')[['mh_provider']])
 }






################################################################################
#```{ Sanity check for the number of visits and subjects }
if(run_heavy_computation){
  rbind(
    
    # Original dataset
    readRDS('data/target_population.rds') %>%
      lapply(X=1:2,Y=.,function(X,Y){
        
        # Visit data
        if(X==1){
          group_by(Y,outcome) %>%
            summarize(n=n()) %>%
            spread(outcome,n) %>%
            mutate(total=event+`non-event`) %>%
            mutate(type='visit')
        
        # Subject data
         }else{
          select(Y,subject_id,outcome) %>%
            .[!duplicated(.),] %>%
            group_by(outcome) %>%
            summarize(n=n()) %>%
            spread(outcome,n) %>%
            mutate(total=event+`non-event`) %>%
            mutate(type='subject')
         }
       }) %>%
      do.call(rbind,.) %>%
      mutate(dataset='original')
    
    # Nationwide medical history
    ,mh_nationwide %>%
      lapply(X=1:2,Y=.,function(X,Y){
        
        # Visit data
        if(X==1){
          Y %>%
            left_join(
              target_population %>%
                select(subject_id,outcome) %>%
                .[!duplicated(.),]
              ,by='subject_id'
            ) %>%
            group_by(outcome) %>%
            summarize(n=n()) %>%
            spread(outcome,n) %>%
            mutate(total=event+`non-event`) %>%
            mutate(type='visit')
        
        # Subject data
         }else{
          Y %>%
            left_join(
              target_population %>%
                select(subject_id,outcome) %>%
                .[!duplicated(.),]
              ,by='subject_id'
            ) %>%
            select(subject_id,outcome) %>%
            .[!duplicated(.),] %>%
            group_by(outcome) %>%
            summarize(n=n()) %>%
            spread(outcome,n) %>%
            mutate(total=event+`non-event`) %>%
            mutate(type='subject')
         }
       }) %>%
      do.call(rbind,.) %>%
      mutate(dataset='nationwide')
    
    # Provider medical history
    ,mh_provider %>%
      lapply(X=1:2,Y=.,function(X,Y){
        
        # Visit data
        if(X==1){
          Y %>%
            left_join(
              target_population %>%
                select(subject_id,outcome) %>%
                .[!duplicated(.),]
              ,by='subject_id'
            ) %>%
            group_by(outcome) %>%
            summarize(n=n()) %>%
            spread(outcome,n) %>%
            mutate(total=event+`non-event`) %>%
            mutate(type='visit')
        
        # Subject data
         }else{
          Y %>%
            left_join(
              target_population %>%
                select(subject_id,outcome) %>%
                .[!duplicated(.),]
              ,by='subject_id'
            ) %>%
            select(subject_id,outcome) %>%
            .[!duplicated(.),] %>%
            group_by(outcome) %>%
            summarize(n=n()) %>%
            spread(outcome,n) %>%
            mutate(total=event+`non-event`) %>%
            mutate(type='subject')
         }
       }) %>%
      do.call(rbind,.) %>%
      mutate(dataset='provider')
  ) %>%
  saveRDS('data/sanity_check.rds')
 }






################################################################################
#```{ Exclude sufficient-sample city for ext. validation }
suppressWarnings(set.seed(33,sample.kind=sample.kind))
extv_city=
  
  # Group visits by city of healthcare provider
  # of which a subject registered to (not a subject visit to)
  target_population %>%
  group_by(reghc_province,reghc_city,subject_id) %>%
  summarize(
    event=as.integer(sum(ifelse(outcome=='event',1,0))>0)
    ,non_event=as.integer(sum(ifelse(outcome=='non-event',1,0))>0)
  ) %>%
  
  # Compute the prevalence of event in each of the cities
  group_by(reghc_province,reghc_city) %>%
  summarize(
    event=sum(event)
    ,non_event=sum(non_event)
  ) %>%
  ungroup() %>%
  mutate(prevalence=event/(event+non_event)) %>%
  
  # Exclude cities and show the prevalence
  lapply(X=1,Y=.,function(X,Y){
    
    # Randomly exclude cities
    ## Filter cities of which event prevalence is between 0 and 1.
    ## Then, sample 22.5% cities for each province without replacement
    Z=filter(Y,prevalence>0 & prevalence<1) %>%
      group_by(reghc_province) %>%
      summarize(extv_city=sample(reghc_city,round(n()*0.225),F))
    
    ## From the sampled cities, sample 40% to overlap with
    ## those of temporal split.
    suppressWarnings(set.seed(33,sample.kind=sample.kind))
    Z %>%
      group_by(reghc_province) %>%
      summarize(extv_city=sample(extv_city,round(n()*0.4),F)) %>%
      ungroup() %>%
      mutate(bgt_city=1) %>%
      right_join(Z,by=c('reghc_province','extv_city')) %>%
      mutate(bgt_city=ifelse(is.na(bgt_city),0,1)) %>%
      arrange(reghc_province,extv_city)
    
   }) %>%
  .[[1]]






################################################################################
#```{ Exclude time-outcome-independent period for ext. validation }
suppressWarnings(set.seed(33,sample.kind=sample.kind))
extv_date=
  
  # Group visits by outcome and month of the event
  target_population %>%
  select(subject_id,event_t,outcome) %>%
  .[!duplicated(.),] %>%
  mutate(month_t=round_date(event_t,unit='month')) %>%
  mutate(x=round((month_t-as_date('2015-01-01'))/dmonths(1))) %>%
  group_by(outcome,month_t,x) %>%
  summarize(subject=n()) %>%
  ungroup() %>%
  
  # Create variable season
  mutate(
    season=case_when(
      between(month_t,as_date('2015-01-01'),as_date('2015-03-20'))~'winter'
      ,between(month_t,as_date('2015-03-21'),as_date('2015-06-21'))~'spring'
      ,between(month_t,as_date('2015-06-22'),as_date('2015-09-23'))~'summer'
      ,between(month_t,as_date('2015-09-24'),as_date('2015-12-22'))~'autumn'
      ,between(month_t,as_date('2015-12-23'),as_date('2016-03-20'))~'winter'
      ,between(month_t,as_date('2016-03-21'),as_date('2016-06-20'))~'spring'
      ,between(month_t,as_date('2016-06-21'),as_date('2016-09-22'))~'summer'
      ,between(month_t,as_date('2016-09-23'),as_date('2016-12-31'))~'autumn'
      ,TRUE~'none'
    )
  ) %>%
  
  # Make sure only including visits on 2015 and 2016
  mutate(year_t=substr(as.character(month_t),1,4)) %>%
  filter(year_t%in%2015:2016) %>%
  
  # For each season, sample a 30-day window to exclude
  lapply(X=seq(sum(!duplicated(paste0(.$year_t,'_',.$season))))
         ,Y=paste0(.$year_t,'_',.$season) %>% .[!duplicated(.)]
         ,Z=.
         ,function(X,Y,Z){
    K=Z %>%
      filter(
        year_t==str_split_fixed(Y[X],'_',2)[,1] &
        season==str_split_fixed(Y[X],'_',2)[,2]
      ) %>%
      mutate(
        min_date=min(month_t)
        ,max_date=ceiling_date(max(month_t),unit='month')-ddays(1)
      )
    L=sample(as_date(K$min_date[1]:K$max_date[1]),1,F)
    
    K %>%
      mutate(min_date=L-days(15),max_date=L+days(15)) %>%
      .[!duplicated(.),]
   }) %>%
  do.call(rbind,.) %>%
  
  # Randomly exclude time-outcome-independent period 
  select(year_t,season,min_date,max_date) %>%
  .[!duplicated(.),] %>%
  filter(!((year_t=='2015' & season=='winter')|
           (year_t=='2016' & season=='autumn')))






################################################################################
#```{ Conduct data partition for int. and ext. validation }
if(run_heavy_computation){
  set=list()
  
  # Select attributes needed for partition
  set$intv=
    target_population %>%
    select(subject_id,outcome,reghc_city,event_t) %>%
    .[!duplicated(.),] %>%
    mutate(outcome=factor(outcome,c('non-event','event')))
  
  # Exclude visits from subjects registered to providers in the excluded cities
  # for geographical split
  set$extv_geo=
    set$intv %>%
    filter(reghc_city %in% extv_city$extv_city)
  
  # Exclude visits from subjects with event date in the excluded period
  # for temporal split
  set$extv_tem=
    set$intv %>%
    lapply(X=1:nrow(extv_date),Y=extv_date,Z=.,function(X,Y,Z){
      filter(Z,between(event_t,Y$min_date[X],Y$max_date[X]))
     }) %>%
    do.call(rbind,.)
  
  # Filter geographically-excluded visits with those of secondly excluded cities
  # and those with event in the excluded period
  # for geotemporal split
  set$extv_bgt=
    set$extv_geo %>%
    lapply(X=1:nrow(extv_date),Y=extv_date,Z=.,function(X,Y,Z){
      Z %>%
        filter(
          reghc_city %in% filter(extv_city,bgt_city==1)$extv_city &
          between(event_t,Y$min_date[X],Y$max_date[X])
        )
     }) %>%
    do.call(rbind,.)
  
  # Hold out the data from any splits
  # for internal validation set
  set$intv=
    set$intv %>%
    filter(!((reghc_city %in% set$extv_geo$reghc_city)|
             (event_t %in% set$extv_tem$event_t)))
  
  # From geographical split, subtract geotemporal split
  set$extv_geo=
    set$extv_geo %>%
    filter(!(reghc_city %in% set$extv_bgt$reghc_city))
  
  # From temporal split, subtract geotemporal split
  set$extv_tem=
    set$extv_tem %>%
    filter(!(event_t %in% set$extv_bgt$event_t))
  
  # From internal validation set, get 5 of 6 parts for training
  # the remaining for validation.
  # This will be cross validation set.
  suppressWarnings(set.seed(33,sample.kind=sample.kind))
  set$cv_idx=
    createDataPartition(set$intv$outcome,times=6,p=5/6,list=F) %>%
    lapply(X=1:ncol(.),Y=.,Z=set$intv,function(X,Y,Z){
      train_set=Z[!seq(nrow(Z))%in%Y[,X],]
      test_set=
        Z[Y[,X],] %>%
        filter(!((reghc_city %in% train_set$reghc_city)|
                 (event_t %in% train_set$event_t)))
      test_set$subject_id
     })
  
  # Summarize prevalence of interval validation set
  set$stats[[1]]=
    data.frame(
      set_subject='Internal validation, cross validation'
      ,event=NA
      ,non_event=NA
      ,prevalence=
        sapply(set$cv_idx,function(x){
          outcome=filter(set$intv,subject_id %in% x)$outcome
          sum(outcome=='event')/length(outcome)
         })
    ) %>%
    group_by(set_subject,event,non_event) %>%
    summarize(
      lb=mean(prevalence)-qnorm(0.975)*sd(prevalence)/sqrt(n())
      ,ub=mean(prevalence)+qnorm(0.975)*sd(prevalence)/sqrt(n())
      ,prevalence=mean(prevalence)
    ) %>%
    ungroup() %>%
    select(set_subject,event,non_event,prevalence,lb,ub)
  
  # Summarize prevalences of external validation sets
  set$stats[[2]]=
    lapply(X=1:4,Y=set,function(X,Y){
      group_by(Y[[X]],outcome) %>%
        summarize(
          set_subject=
            c('Internal validation'
              ,paste('External validation,',c('geographical'
                                              ,'temporal'
                                              ,'geotemporal'),'split')
              )[X]
          ,n=n()
        ) %>%
        spread(outcome,n) %>%
        ungroup() %>%
        rename(non_event=`non-event`) %>%
        mutate(prevalence=event/(event+non_event)) %>%
        mutate(lb=NA,ub=NA)
     }) %>%
    do.call(rbind,.)
  
  # Summarize prevalences by subject
  set$stats[['subject']]=
    set$stats %>%
    do.call(rbind,.) %>%
    .[c(2,1,3:5),]
  
  # Summarize visits by validation set
  set$stats[['visit']]=
    lapply(
      X=1:3
      ,Y=c(paste('External validation,',c('geographical'
                                           ,'temporal'
                                           ,'geotemporal'),'split')
           )
      ,Z=c(paste0('extv_',c('geo','tem','bgt')))
      ,FUN=function(X,Y,Z){
        data.frame(
          set_visit=Y[X]
          ,event=sum(
            mh_nationwide$subject_id %in%
              set[[Z[X]]]$subject_id[set[[Z[X]]]$outcome=='event']
          )
          ,non_event=sum(
            mh_nationwide$subject_id %in%
              set[[Z[X]]]$subject_id[set[[Z[X]]]$outcome=='non-event']
          )
        )
     }) %>%
    do.call(rbind,.) %>%
    
    # Add row for total visit of external validation sets
    add_row(
      set_visit='External validation'
      ,event=
        sum(filter(.,str_detect(set_visit,'External validation'))$event)
      ,non_event=
        sum(filter(.,str_detect(set_visit,'External validation'))$non_event)
    ) %>%
    
    filter(set_visit!='Internal validation') %>%
    
    # Add row for visit of internal validation set
    add_row(
      set_visit='Internal validation'
      ,event=
        (sum(target_population$outcome=='event')-
         filter(.,set_visit=='External validation')$event)
      ,non_event=
        (sum(target_population$outcome=='non-event')-
         filter(.,set_visit=='External validation')$non_event)
    ) %>%
    
    # Add row for visit of internal validation set, training split
    add_row(
      set_visit='Internal validation, training split'
      ,event=
        round(5/6*filter(.,set_visit=='Internal validation')$event)
      ,non_event=
        round(5/6*filter(.,set_visit=='Internal validation')$non_event)
    ) %>%
    
    # Add row for visit of internal validation set, validation split
    add_row(
      set_visit='Internal validation, validation split'
      ,event=
        filter(.,set_visit=='Internal validation')$event-
        filter(.,set_visit=='Internal validation, training split')$event
      ,non_event=
        filter(.,set_visit=='Internal validation')$non_event-
        filter(.,set_visit=='Internal validation, training split')$non_event
    ) %>%
    
    # Compute proportion of each set
    mutate(
      subtotal=event+non_event
      ,total=nrow(target_population)
    ) %>%
    mutate(
      p=subtotal/total
    ) %>%
  
  # Show the number and proportion of data partition
  arrange(factor(set_visit,c(
    'Internal validation'
    ,'Internal validation, training split'
    ,'Internal validation, validation split'
    ,'External validation'
    ,'External validation, geographical split'
    ,'External validation, temporal split'
    ,'External validation, geotemporal split'
  )))
  
  set$stats[1:2]=NULL
  
  saveRDS(set,'data/set.rds')
 }else{
  set=readRDS('data/set.rds')
 }






################################################################################
#```{ Create predictors table for codes and the descriptions }
if(run_heavy_computation){
  predictors=list()
  
  # Get column names form nationwide medical history
  predictors[[1]]=data.frame(predictor=colnames(mh_nationwide)[-1:-2])
  
  # Overlap the names with code description in discharge (primary) diagnosis
  predictors[[2]]=
    predictors[[1]] %>%
    inner_join(
      readRDS('data/target_population.rds') %>%
        select(icd10_code,icd10_3mer_desc) %>%
        setNames(c('predictor','description')) %>%
        .[!duplicated(.),]
      ,by='predictor'
    )
  
  # Overlap the names with code description in admission diagnosis
  predictors[[3]]=
    predictors[[1]] %>%
    inner_join(
      readRDS('data/admission_diagnosis.rds') %>%
        select(admission_icd10_code,admission_icd10_3mer_desc) %>%
        setNames(c('predictor','description')) %>%
        .[!duplicated(.),]
      ,by='predictor'
    )
  
  # Overlap the names with code description in discharge (secondary) diagnosis
  predictors[[4]]=
    predictors[[1]] %>%
    inner_join(
      readRDS('data/secondary_diagnoses.rds') %>%
        select(secondary_icd10_code,secondary_icd10_3mer_desc) %>%
        setNames(c('predictor','description')) %>%
        .[!duplicated(.),]
      ,by='predictor'
    )
  
  # Overlap the names with code description in procedure
  predictors[[5]]=
    predictors[[1]] %>%
    inner_join(
      suppressWarnings(separate(
        readRDS('data/procedures.rds') %>%
          select(icd9_code_desc) %>%
          .[!duplicated(.),,drop=F]
        ,icd9_code_desc,c('predictor','description')
        ,sep=' - '
      ))
      ,by='predictor'
    )
  
  # Get the overlaps only
  predictors=
    predictors[2:5] %>%
    do.call(rbind,.) %>%
    .[!duplicated(.),] %>%
    arrange(predictor)
  
  # Add non-medical history predictors
  predictors=
    readRDS('data/target_population.rds') %>%
    select(
      marital_status
      ,insurance_class
      ,occupation_segment
      ,healthcare_class
    ) %>%
    gather(predictor,value) %>%
    .[!duplicated(.),] %>%
    filter(!is.na(value)) %>%
    mutate(value=str_replace_all(value,'[:punct:]|[:space:]','_')) %>%
    unite(predictor,predictor,value,sep='.') %>%
    mutate(
      description=
        sapply(
          predictor
          ,function(x)paste(str_split(x,'_|\\.')[[1]],collapse=' ')
        )
    ) %>%
    rbind(data.frame(predictor='age',description='age')) %>%
    rbind(predictors)
  
  saveRDS(predictors,'data/predictors.rds')
 }else{
  predictors=readRDS('data/predictors.rds')
 }






################################################################################
#```{ Load DAG list for causal diagram and inference }
dag=readRDS('data/dag.rds')










################################################################################
#```{ Conduct feature representation based on causal diagram }
if(run_heavy_computation){
  
  cat('Conduct feature representation based on causal diagram\n')
  cat('Started:',as.character(now()),'\n')
  
  set$inference=
    
    # Standardize age, then normalize -1.96 to 1.96 into 0 to 1
    readRDS('data/target_population.rds') %>%
    select(subject_id,day_to_event,outcome,age) %>%
    .[!duplicated(.),] %>%
    filter(subject_id %in% set$intv$subject_id) %>%
    mutate(
      age_m=mean(age)
      ,age_s=sd(age)
      ,age=(qnorm(0.975)+(age-age_m)/age_s)/(2*qnorm(0.975))
    ) %>%
    
    # Add non-medical history categorical predictors
    left_join(
      
      # Group by subject and summarize the predictors
      readRDS('data/target_population.rds') %>%
        select(
          subject_id
          ,marital_status
          ,insurance_class
          ,occupation_segment
          ,healthcare_class
        ) %>%
        .[!duplicated(.),] %>%
        group_by(subject_id) %>%
        summarize_all(function(x){
          x[!duplicated(x)] %>%
            .[which.max(
              sapply(x[!duplicated(x)],function(x2)sum(x2==x,na.rm=T))
            )]
         }) %>%
        
        # Spread each categorical predictor as several binary predictors
        lapply(X=seq(ncol(.)),Y=.,function(X,Y){
          if(colnames(Y)[X] %in% c('marital_status'
                                   ,'insurance_class'
                                   ,'occupation_segment'
                                   ,'healthcare_class')){
            Y[,c('subject_id',colnames(Y)[X]),drop=F] %>%
              rename_at(colnames(Y)[X],function(x)'key') %>%
              mutate(
                key=str_replace_all(key,'[:punct:]|[:space:]','_')
                ,value=1
              ) %>%
              mutate_at('key',function(x)paste0(colnames(Y)[X],'.',x)) %>%
              spread(key,value) %>%
              select(-subject_id)
           }else{
            Y[,colnames(Y)[X],drop=F]
           }
         }) %>%
        do.call(cbind,.) %>%
        
        # If a category was applied (non-missing), then 1; otherwise 0
        mutate_at(
          colnames(.) %>% .[.!='subject_id']
          ,function(x)as.integer(!is.na(x))
        )
      ,by='subject_id'
    ) %>%
    
    # Join with nationwide medical history with censoring probability up to 2 y
    right_join(
      mh_nationwide %>%
        filter(subject_id %in% set$intv$subject_id) %>%
        group_by(subject_id) %>%
        mutate(cens_p=(365*2-max(day_to_event))/(365*2)) %>%
        ungroup()
      ,by=c('subject_id','day_to_event')
    ) %>%
    
    # Convert outcome as 0 and 1, then re-arrange columns
    mutate(outcome=as.integer(outcome=='event')) %>%
    select(subject_id,day_to_event,outcome,cens_p,everything()) %>%
    
    # Make  a sequence column and select codes
    # that were used for making causal predictors
    mutate(seq=seq(nrow(.))) %>%
    select_at(unique(c(
      'subject_id','day_to_event','outcome','cens_p','age_m','age_s','seq'
      ,colnames(.) %>%
        .[str_detect(str_to_upper(.)
                     ,paste(dag$measure_nodes$name,collapse='|')
        )]
    ))) %>%
    
    # Stack predictor columns
    pivot_longer(
      colnames(.) %>%
        .[!.%in%c('subject_id'
                  ,'day_to_event'
                  ,'outcome'
                  ,'cens_p'
                  ,'age_m'
                  ,'age_s'
                  ,'age'
                  ,'seq')]
      ,names_to='predictor',values_to='value'
    ) %>%
    as.data.frame() %>%
    
    # Make causal predictors from the codes,
    # but set the values as 0 (no) or 1 (yes), instead of a frequency
    pblapply(X=seq(nrow(dag$measure_nodes))
             ,Y=dag$measure_nodes$name
             ,Z=dag$measure_nodes$label
             ,K=.
             ,function(X,Y,Z,K){
               if(!(Z[X]=='Y01*'|Y[X]=='AGE')){
                 filter(K,str_detect(str_to_upper(predictor),Y[X])) %>%
                   group_by(
                     subject_id
                     ,day_to_event
                     ,outcome
                     ,cens_p
                     ,age_m
                     ,age_s
                     ,age
                     ,seq
                   ) %>%
                   summarize(
                     predictor=str_remove_all(Z[X],'\\*')
                     ,value=as.integer(sum(value)>0)
                   ) %>%
                   ungroup()
                }
              }) %>%
    do.call(rbind,.) %>%
    
    # If yes, turn values as probabilities up to 2 y
    mutate(value=value*(365*2-day_to_event)/(365*2)) %>%
    
    # Spread the predictors as individual columns
    spread(predictor,value) %>%
    select(-seq) %>%
    
    # Finishing touch
    rename(A20=age) %>%
    select_at(colnames(.) %>% .[order(.)]) %>%
    select(subject_id,day_to_event,outcome,cens_p,age_m,age_s,everything())
  
  cat('End:',as.character(now()))
  saveRDS(set$inference,'data/inference.rds')
 }else{
  cat(readRDS('data/log.rds')[['inference']])
  set$inference=readRDS('data/inference.rds')
 }






################################################################################
#```{ Prepare formula + binary data + G-estimation function }
# Formula
dag$formula$A02=
  outcome~
  A02+A20+A25 # backdoor: A21
dag$formula$A03=
  outcome~
  A03+A15+L37+A28+A04+A05+A10+A11+A13+A02+
  L30+A20+A24+L34+A14+L33+L35+A25+A19 # backdoor: A08, A21, A23, A22
dag$formula$A04=
  outcome~
  A04+A02+A05+A10+A11+A13+L30+A20+A25+A24+A14+
  L34+L33+L35+L37+L38+A15+A28+A19 # backdoor: A21, A23, A26, A22, A08
dag$formula$A05=
  outcome~
  A05+A20+A10+A24+L30+A14+L34+A19+A15+A28 # backdoor: A21, A23, A26
dag$formula$A06=
  outcome~
  A06+A19+A09+A20+A25+L30+A05+A10+A24+
  L34+A14+A19+A15+A28 # backdoor: A08, A21, A27, A22, A23, A26
dag$formula$A09=
  outcome~
  A09+A19+A20 # backdoor: A08, A21, A22, A23
dag$formula$A10=
  outcome~
  A10+L34+A14+A19+A15+A28 # backdoor: A21
dag$formula$A11=
  outcome~
  A11+A15+A20+L33+
  L34+L35+L37+A28 # backdoor: A22, A23, A08
dag$formula$A12=
  outcome~
  A12+A02+A25+
  L34+L35+A20 # backdoor: A21
dag$formula$A13=
  outcome~
  A13+A15+A28+
  L35+L37+L38 # backdoor:
dag$formula$A14=
  outcome~
  A14+A15+A28+A19+
  L37 # backdoor: A21
dag$formula$A15=
  outcome~
  A15+A28+
  L37 # backdoor: 
dag$formula$A19=
  outcome~
  A19 # backdoor: 
dag$formula$A20=
  outcome~
  A20 # backdoor: 
dag$formula$A24=
  outcome~
  A24 # backdoor: 
dag$formula$A25=
  outcome~
  A25 # backdoor:
dag$formula$A28=
  outcome~
  A28 # backdoor:

# Binary data
bn_exp_data=
  set$inference %>%
  lapply(X=1:2,Y=.,function(X,Y){
    if(X==1){
      select(Y,outcome,A20)
     }else{
      Y %>%
        select(
          -subject_id
          ,-day_to_event
          ,-outcome
          ,-cens_p
          ,-age_m
          ,-age_s
          ,-A20
        ) %>%
        mutate_all(function(x)as.integer(x>0))
     }
   }) %>%
  do.call(cbind,.) %>%
  select_at(c('outcome',colnames(.) %>% .[.!='outcome'] %>% .[order(.)])) %>%
  mutate(
    A20=as.integer(!between(
      A20
      ,(18-set$inference$age_m[1])/set$inference$age_s[1]
      ,(35-set$inference$age_m[1])/set$inference$age_s[1]
    ))
  )

# G-estimation function
g_estimation=function(formula
                      ,data
                      ,R=5
                      ,par=0
                      ,lower=-10
                      ,upper=10
                      ,parallel.run=T
                      ,ncpus=detectCores()-1
                      ,verbose=T){
  AbsCoefHpsi=function(formula,data){
    Hpsi=function(outcome,exposure){function(psi){outcome-psi*exposure } }
    function(psi){
      outcome=data[[as.character(formula)[2]]]
      exposure=data[[str_split(as.character(formula)[3],' \\+ ')[[1]][1]]]
      cov_names=c(str_split_fixed(as.character(formula)[3],' \\+ ',2))
      cov_names[2]=ifelse(cov_names[2]=='',1,cov_names[2])
      exposure.formula=eval(parse(text=paste0(cov_names,collapse=' ~ ')))
      model=suppressWarnings(
        glm(exposure.formula,data,family=binomial(link='logit'))
      )
      data$..Hpsi..=(Hpsi(outcome,exposure))(psi)
      model2=suppressWarnings(update(model,formula=~.+..Hpsi..))
      tail(coef(model2),1) %>% abs
     }
   }
  bs_func=function(X,formula,data,AbsCoefHpsi){
    data=data %>% .[sample(seq(nrow(.)),nrow(.),T),]
    AbsCoefHpsiFun=AbsCoefHpsi(formula,data)
    results=suppressWarnings(
      optim(par=par,AbsCoefHpsiFun,lower=lower,upper=upper)
    )
    c(results$par
      ,results$value
      ,results$counts[1]
      ,results$counts[2]
      ,results$convergence)
   }
  if(parallel.run){
    cl=makeCluster(ncpus)
    clusterEvalQ(cl,{
      library('parallel')
      library('tidyverse')
     })
    b=pbsapply(
      X=seq(R)
      ,formula=formula
      ,data=data
      ,AbsCoefHpsi=AbsCoefHpsi
      ,cl=cl
      ,bs_func
    )
    stopCluster(cl)
   }else{
    if(verbose){
      b=pbsapply(
        X=seq(R)
        ,formula=formula
        ,data=data
        ,AbsCoefHpsi=AbsCoefHpsi
        ,bs_func
      )
     }else{
      b=sapply(
        X=seq(R)
        ,formula=formula
        ,data=data
        ,AbsCoefHpsi=AbsCoefHpsi
        ,bs_func
      )
     }
   }
  t(b) %>%
    `dimnames<-`(
      list(NULL,c('par','psi','counts.function','counts.gradient','converge'))
    ) %>%
    as.data.frame() %>%
    lapply(X=1,Y=.,function(X,Y){
      rbind(
        Y %>%
          summarize_all(mean) %>%
          mutate(metric='G_estimate')
        ,Y %>%
          summarize_all(function(x){
            mean(x)-qnorm(0.975)*sd(x)/sqrt(length(x))
           }) %>%
          mutate(metric='G_lb')
        ,Y %>%
          summarize_all(function(x){
            mean(x)+qnorm(0.975)*sd(x)/sqrt(length(x))
           }) %>%
          mutate(metric='G_ub')
      ) %>%
        column_to_rownames(var='metric') %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column(var='metric')
     }) %>%
    .[[1]] %>%
    mutate(OR=exp(G_estimate),OR_lb=exp(G_lb),OR_ub=exp(G_ub)) %>%
    mutate(Pr=OR/(1+OR),Pr_lb=OR_lb/(1+OR_lb),Pr_ub=OR_ub/(1+OR_ub))
 }






################################################################################
#```{ Conduct logistic regression for causal inference }
dag$glm=
  dag$formula %>%
  lapply(glm,binomial(link='logit'),bn_exp_data)

dag$coef=
  dag$glm %>%
  lapply(tidy) %>%
  lapply(
    mutate
    ,OR=exp(estimate)
    ,OR_lb=exp(estimate-qnorm(0.975)*std.error)
    ,OR_ub=exp(estimate+qnorm(0.975)*std.error)
    ,Pr=OR/(1+OR)
    ,Pr_lb=OR_lb/(1+OR_lb)
    ,Pr_ub=OR_ub/(1+OR_ub)
  )


#```{ Conduct G-estimation for causal inference ... }
if(run_heavy_computation){
  cat('Conduct G-estimation for causal inference',append=T)
  cat('based on the theoretical causal diagram\n')
  cat('Started:',as.character(now()),'\n')
  
  suppressWarnings(set.seed(33,sample.kind=sample.kind))
  dag$gest=
    dag$formula %>%
    pblapply(g_estimation,data=bn_exp_data,R=30,ncpus=detectCores()-1)
  
  cat('End:',as.character(now()))
  saveRDS(dag$gest,'data/g_estimation.rds')
 }else{
  cat(readRDS('data/log.rds')[['g_estimation']])
  dag$gest=readRDS('data/g_estimation.rds')
 }

dag$sig=
  lapply(X=1:2,Y=dag$coef,Z=dag$gest,function(X,Y,Z){
    if(X==1){
      sapply(Y,function(x){
        ifelse(
          between(1,round(x$OR_lb[2],4),round(x$OR_ub[2],4))|is.na(x$OR[2])
          ,0,1
        )
       })
     }else{
      sapply(Z,function(x){
        ifelse(
          between(1,round(x$OR_lb[2],4),round(x$OR_ub[2],4))|is.na(x$OR[2])
          ,0,1
        )
       })
     }
   }) %>%
  do.call(rbind,.) %>%
  `rownames<-`(c('glm','gest'))





################################################################################
#```{ Construct a class-balanced training set by provider ... }
if(run_heavy_computation){
  cat('Construct a class-balanced training set by provider',append=T)
  cat('with all and represented features\n')
  cat('Started:',as.character(now()),'\n')
  
  pb=startpb(0,5)
  on.exit(closepb(pb))
  setpb(pb,0)
  
  # Use provider medical history
  setpb(pb,1)
  temp=
    mh_provider %>%
    select(-subject_id,-day_to_event,-est_strata_id) %>%
    t() %>%
    t()
  
  temp=
    (temp>0) %>%
    `storage.mode<-`('integer')
  
  # Make the confirmed causal predictors defined by provider medical histories
  setpb(pb,2)
  temp2=
    lapply(
      X=seq(sum(
        str_remove_all(dag$measure_nodes$label,'\\*') %in%
          colnames(dag$sig)[dag$sig[2,]==1]
      ))
      ,Y=
        dag$measure_nodes %>%
        filter(
          str_remove_all(label,'\\*') %in% colnames(dag$sig)[dag$sig[2,]==1]
        ) %>%
        pull(name)
      ,Z=
        dag$measure_nodes %>%
        filter(
          str_remove_all(label,'\\*') %in% colnames(dag$sig)[dag$sig[2,]==1]
        ) %>%
        pull(label)
      ,K=temp
      ,function(X,Y,Z,K){
        K=K %>%
          .[,str_detect(str_to_upper(colnames(.)),Y[X]),drop=F] %>%
          rowSums() %>%
          matrix() %>%
          `dimnames<-`(
            list(rownames(K),paste0('causal_',str_remove_all(Z[X],'\\*')))
          )
        
        K=(K>0) %>%
          `storage.mode<-`('integer')
        
       }) %>%
    do.call(cbind,.)
  
  # Combine medical history and the causal predictors,
  # then convert frequency into probability or proportion of day being
  # encountered up to 2 years before the event
  setpb(pb,3)
  temp=
    cbind(temp2,temp) %>%
    sweep(1,((365*2-mh_provider$day_to_event)/(365*2)),'*')
  rm(temp2)
  
  
  # Join the results with other attributes for data partition
  temp=
    mh_provider %>%
    select(subject_id,day_to_event,est_strata_id) %>%
    cbind(as.data.frame(temp))
  
  # Get non-medical history predictors
  setpb(pb,4)
  set$training=
    
    ## Numerical predictors (age), standardize, and normalize intor 0 to 1
    readRDS('data/target_population.rds') %>%
    select(subject_id,day_to_event,outcome,age) %>%
    .[!duplicated(.),] %>%
    filter(subject_id %in% set$intv$subject_id) %>%
    mutate(
      age_m=mean(age)
      ,age_s=sd(age)
      ,age=(qnorm(0.975)+(age-age_m)/age_s)/(2*qnorm(0.975))
    ) %>%
    
    ## Categorical predictors and getting these cleaned
    left_join(
      readRDS('data/target_population.rds') %>%
        select(
          subject_id
          ,marital_status
          ,insurance_class
          ,occupation_segment
          ,healthcare_class
        ) %>%
        .[!duplicated(.),] %>%
        group_by(subject_id) %>%
        summarize_all(function(x){
          x[!duplicated(x)] %>%
            .[which.max(
              sapply(x[!duplicated(x)],function(x2)sum(x2==x,na.rm=T))
            )]
         }) %>%
        lapply(X=seq(ncol(.)),Y=.,function(X,Y){
          if(colnames(Y)[X] %in% c('marital_status'
                                   ,'insurance_class'
                                   ,'occupation_segment'
                                   ,'healthcare_class')){
            Y[,c('subject_id',colnames(Y)[X]),drop=F] %>%
              rename_at(colnames(Y)[X],function(x)'key') %>%
              mutate(
                key=str_replace_all(key,'[:punct:]|[:space:]','_'),value=1
              ) %>%
              mutate_at('key',function(x)paste0(colnames(Y)[X],'.',x)) %>%
              spread(key,value) %>%
              select(-subject_id)
           }else{
            Y[,colnames(Y)[X],drop=F]
           }
         }) %>%
        do.call(cbind,.) %>%
        mutate_at(
          colnames(.) %>% .[.!='subject_id']
          ,function(x)as.integer(!is.na(x))
        )
      ,by='subject_id'
    ) %>%
    
    ## Subset internal validation set and add censoring probability
    right_join(
      temp %>%
        filter(subject_id %in% set$intv$subject_id) %>%
        group_by(subject_id) %>%
        mutate(cens_p=(365*2-max(day_to_event))/(365*2)) %>%
        ungroup()
      ,by=c('subject_id','day_to_event')
    ) %>%
    
    # Convert outcome as 0 and 1
    mutate(outcome=as.integer(outcome=='event')) %>%
    
    # Pre-finishing touch
    mutate(causal_A20=age) %>%
    select(subject_id,day_to_event,est_strata_id,outcome,cens_p,everything())
  rm(temp)
  
  
  # Upsample minority to balance the outcome
  setpb(pb,5)
  suppressWarnings(set.seed(33,sample.kind=sample.kind))
  set$training=
    set$training %>%
    upSample(factor(.$outcome)) %>%
    mutate(
      Class=
        ifelse(Class==0,'non_event','event') %>%
        factor(c('event','non_event'))
    ) %>%
    select(
      subject_id
      ,day_to_event
      ,outcome
      ,Class
      ,age_m
      ,age_s
      ,cens_p
      ,everything()
    )
  
  cat('\nEnd:',as.character(now()))
  rm(pb)
  saveRDS(set$training,'data/training.rds')
 }else{
  cat(readRDS('data/log.rds')[['training']])
 }






################################################################################
#```{ Create an empty list to save model }
if(run_heavy_computation){
  model=list()
 }






################################################################################
#```{ Conduct ridge regression by parallel computing }
if(run_heavy_computation){
  cat('Conduct ridge regression by parallel computing\n')
  cat('Started:',as.character(now()),'\n')
  pb=startpb(0,3)
  on.exit(closepb(pb))
  setpb(pb,0)
  cl=makePSOCKcluster(detectCores()-1)
  registerDoParallel(cl)
  clusterEvalQ(cl,{
    library('parallel')
    library('doParallel')
    library('tidyverse')
    library('pbapply')
    library('caret')
    library('glmnet')
   })
  
  setpb(pb,1)
  suppressWarnings(set.seed(33,sample.kind=sample.kind))
  model$ridge=
    train(
      Class~.
      ,data=
        set$training %>%
        .[,c('Class','cens_p',colnames(.) %>% .[str_detect(.,'causal')])]
      ,method='glmnet'
      ,metric='ROC'
      ,trControl=
        trainControl(
          'cv'
          ,number=6
          ,summaryFunction=twoClassSummary
          ,classProbs=T
          ,savePredictions=T
          ,allowParallel=T
        )
      ,tuneGrid=data.frame(alpha=0,lambda=10^seq(-9,0,len=10))
    )
  
  setpb(pb,2)
  model$ridge=
    train(
      Class~.
      ,data=
        set$training %>%
        .[,c('Class','cens_p',colnames(.) %>% .[str_detect(.,'causal')])]
      ,method='glmnet'
      ,metric='ROC'
      ,trControl=
        trainControl(
          method='boot'
          ,number=30
          ,summaryFunction=twoClassSummary
          ,classProbs=T
          ,savePredictions=T
          ,allowParallel=T
        )
      ,tuneGrid=data.frame(alpha=0,lambda=model$ridge$bestTune$lambda)
    )
  
  stopCluster(cl)
  registerDoSEQ()
  rm(cl)
  gc()
  setpb(pb,3)
  cat('\nEnd:',as.character(now()))
  rm(pb)
  saveRDS(model$ridge,paste0('data/ridge.rds'))
 }else{
  cat(readRDS('data/log.rds')[['ridge']])
 }






################################################################################
#```{ Only use predictor with no perfect separation }
if(run_heavy_computation){
  set$no_per_sep=
    set$training %>%
    select(-subject_id,-day_to_event,-outcome,-age_m,-age_s,-est_strata_id) %>%
    group_by(Class) %>%
    summarize_all(sd) %>%
    gather(predictor,value,-Class) %>%
    spread(Class,value) %>%
    filter(!(event==0 | non_event==0))
 }






################################################################################
#```{ Conduct NPS cross-validated PCA }
if(run_heavy_computation){
  cat('Conduct NPS cross-validated PCA\n')
  cat('Started:',as.character(now()),'\n')
  
  model$pca_nps=
    set$training %>%
    .[,c('subject_id',set$no_per_sep$predictor)] %>%
    pblapply(X=seq(length(set$cv_idx)),Y=set$cv_idx,Z=.,function(X,Y,Z){
      K=Z %>%
        filter(!subject_id %in% Y[[X]]) %>%
        select(-subject_id)
      
      K=K[,K %>%
            summarise_all(sd) %>%
            gather() %>%
            filter(value!=0) %>%
            pull(key)
      ]
      
      L=summarize_all(K,mean) %>%
        gather()
      L=setNames(L$value,L$key)
      M=summarize_all(K,sd) %>%
        gather()
      M=setNames(M$value,M$key)
      N=as.matrix(K) %>%
        sweep(2,L,'-') %>%
        sweep(2,M,'/') %>%
        prcomp()
      
      list(mean=L,sd=M,prcomp=N)
     })
  
  gc()
  cat('End:',as.character(now()))
  saveRDS(model$pca_nps,'data/pca_nps.rds')
 }else{
  cat(readRDS('data/log.rds')[['pca_nps']])
 }






################################################################################
#```{ Save min. NPS PCs to get 50% cum. pct. of var. explained }
if(run_heavy_computation){
  pca_nps_pve50=
    model$pca_nps %>%
    sapply(X=seq(length(.)),Y=.,function(X,Y){
      cum_pve=cumsum(Y[[X]]$prcomp$sdev^2/sum(Y[[X]]$prcomp$sdev^2))
      print(plot(cum_pve))
      which(cum_pve>=0.5)[1]
     }) %>%
    max()
 }






################################################################################
#```{ A function to represent predictors as the CV PCA }
pc_converter=function(data,pca,npc=NULL){
  
  pb=startpb(0,7)
  on.exit(closepb(pb))
  setpb(pb,0)
  
  # Get rotation of the PC matrix (PC weights as columns)
  if(is.null(npc)) npc=seq(ncol(pca[[1]]$prcomp$rotation))
  
  # Subset the top PCs and stack all versions of PCs
  setpb(pb,1)
  rotated_pc=
    pca %>%
    lapply(X=seq(length(.)),Y=.,function(X,Y){
      as.data.frame(Y[[X]]$prcomp$rotation[,npc]) %>%
        rownames_to_column(var='predictor')
     }) %>%
    do.call(rbind,.)
  
  # Group by predictor and average the PC weights over all versions
  setpb(pb,2)
  rotated_pc=
    rotated_pc %>%
    group_by(predictor) %>%
    summarize_all(function(x)mean(x,na.rm=T)) %>%
    ungroup() %>%
    column_to_rownames(var='predictor')
  
  # Get mean and standard deviation for each version
  setpb(pb,3)
  scaler=
    pca %>%
    lapply(X=seq(length(.)),Y=.,function(X,Y){
      data.frame(predictor=names(Y[[X]]$mean),mean=Y[[X]]$mean,sd=Y[[X]]$sd) %>%
        `rownames<-`(NULL)
     }) %>%
    do.call(rbind,.)
  
  # Average the mean and standard deviation of all versions
  setpb(pb,4)
  scaler=
    scaler %>%
    group_by(predictor) %>%
    summarize_all(function(x)mean(x,na.rm=T)) %>%
    ungroup() %>%
    column_to_rownames(var='predictor')
  
  # Standardized predictors using the estimate of mean and standard deviation
  setpb(pb,5)
  scaled_data=
    data %>%
    .[,rownames(rotated_pc)] %>%
    sweep(2,scaler$mean,'-') %>%
    sweep(2,scaler$sd,'/')
  
  # Compute dot product of standardized predictors with the averaged rotated PCs
  setpb(pb,6)
  pc_data=as.matrix(scaled_data) %*% as.matrix(rotated_pc)
  
  setpb(pb,7)
  as.data.frame(pc_data)
  
 }






################################################################################
#```{ Represented NPS predictors as the cross-validated PCA }
if(run_heavy_computation){
  cat('Represent NPS predictors as the cross-validated PCA\n')
  cat('Started:',as.character(now()),'\n')
  
  set$training_pc_nps=pc_converter(set$training,model$pca_nps,npc=1:pca_nps_pve50)
  
  gc()
  cat('End:',as.character(now()))
  saveRDS(set$training_pc_nps,'data/training_pc_nps.rds')
 }else{
  cat(readRDS('data/log.rds')[['training_pc_nps']])
 }






################################################################################
#```{ Conduct NPS-PC elastic net regression by parallel computing }
if(run_heavy_computation){
  cat('Conduct NPS-PC elastic net regression by parallel computing\n')
  cat('Started:',as.character(now()),'\n')
  pb=startpb(0,3)
  on.exit(closepb(pb))
  setpb(pb,0)
  cl=makePSOCKcluster(floor(detectCores()/2))
  registerDoParallel(cl)
  clusterEvalQ(cl,{
    library('parallel')
    library('doParallel')
    library('tidyverse')
    library('pbapply')
    library('caret')
    library('glmnet')
   })
  
  setpb(pb,1)
  suppressWarnings(set.seed(33,sample.kind=sample.kind))
  model$pc_nps_elnet=
    train(
      Class~.
      ,data=
        set$training %>%
        select(Class,cens_p) %>%
        cbind(set$training_pc_nps)
      ,method='glmnet'
      ,metric='ROC'
      ,trControl=
        trainControl(
          'cv'
          ,number=6
          ,summaryFunction=twoClassSummary
          ,classProbs=T
          ,savePredictions=T
          ,allowParallel=T
        )
      ,tuneLength=10
    )
  
  setpb(pb,2)
  model$pc_nps_elnet=
    train(
      Class~.
      ,data=
        set$training %>%
        select(Class,cens_p) %>%
        cbind(set$training_pc_nps)
      ,method='glmnet'
      ,metric='ROC'
      ,trControl=
        trainControl(
          method='boot'
          ,number=30
          ,summaryFunction=twoClassSummary
          ,classProbs=T
          ,savePredictions=T
          ,allowParallel=T
        )
      ,tuneGrid=
        data.frame(
          alpha=model$pc_nps_elnet$bestTune$alpha
          ,lambda=model$pc_nps_elnet$bestTune$lambda
        )
    )
  
  stopCluster(cl)
  registerDoSEQ()
  rm(cl)
  gc()
  setpb(pb,3)
  cat('\nEnd:',as.character(now()))
  rm(pb)
  saveRDS(model$pc_nps_elnet,'data/pc_nps_elnet.rds')
 }else{
  cat(readRDS('data/log.rds')[['pc_nps_elnet']])
 }






################################################################################
#```{ The selected NPS-PC by the elastic net regression }
if(run_heavy_computation){
  set$selected_pc=
    coef(model$pc_nps_elnet$finalModel,model$pc_nps_elnet$bestTune$lambda) %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column(var='predictor') %>%
    mutate(predictor=str_remove_all(predictor,'\\`|\\(|\\)')) %>%
    setNames(c('predictor','estimate')) %>%
    filter(estimate!=0) %>%
    arrange(desc(abs(estimate))) %>%
    filter(!predictor %in% c('Intercept','cens_p')) %>%
    slice(1:20) %>%
    mutate(idx=str_remove_all(predictor,'PC') %>% as.integer())
  
  saveRDS(set$selected_pc,'data/selected_pc.rds')
 }else{
  set$selected_pc=readRDS('data/selected_pc.rds')
 }






################################################################################
#```{ Represented NPS predictors as the selected and CV PCA }
if(run_heavy_computation){
  cat('Represent NPS predictors as the selected and cross-validated PCA\n')
  cat('Started:',as.character(now()),'\n')
  
  set$training_spc_nps=pc_converter(set$training,model$pca_nps,npc=sort(set$selected_pc$idx))
  
  gc()
  cat('End:',as.character(now()))
  saveRDS(set$training_spc_nps,'data/training_spc_nps.rds')
 }else{
  cat(readRDS('data/log.rds')[['training_spc_nps']])
 }






################################################################################
#```{ Conduct selected NPS-PC random forest by parallel computing }
if(run_heavy_computation){
  cat('Conduct selected NPS-PC random forest by parallel computing\n')
  cat('Started:',as.character(now()),'\n')
  pb=startpb(0,3)
  on.exit(closepb(pb))
  setpb(pb,0)
  cl=makePSOCKcluster(floor(detectCores()/2))
  registerDoParallel(cl)
  clusterEvalQ(cl,{
    library('parallel')
    library('doParallel')
    library('tidyverse')
    library('pbapply')
    library('caret')
    library('Rborist')
   })
  
  setpb(pb,1)
  suppressWarnings(set.seed(33,sample.kind=sample.kind))
  model$spc_nps_rf=
    train(
      Class~.
      ,data=
        set$training %>%
        select(Class,cens_p) %>%
        cbind(set$training_spc_nps)
      ,method='Rborist'
      ,metric='ROC'
      ,trControl=
        trainControl(
          'cv'
          ,number=6
          ,summaryFunction=twoClassSummary
          ,classProbs=T
          ,savePredictions=T
          ,allowParallel=T
        )
      ,tuneLength=10
    )
  
  setpb(pb,2)
  model$spc_nps_rf=
    train(
      Class~.
      ,data=
        set$training %>%
        select(Class,cens_p) %>%
        cbind(set$training_spc_nps)
      ,method='Rborist'
      ,metric='ROC'
      ,trControl=
        trainControl(
          method='boot'
          ,number=30
          ,summaryFunction=twoClassSummary
          ,classProbs=T
          ,savePredictions=T
          ,allowParallel=T
        )
      ,tuneGrid=
        data.frame(
          predFixed=model$spc_nps_rf$bestTune$predFixed
          ,minNode=model$spc_nps_rf$bestTune$minNode
        )
    )
  
  stopCluster(cl)
  registerDoSEQ()
  rm(cl)
  gc()
  setpb(pb,3)
  cat('\nEnd:',as.character(now()))
  rm(pb)
  saveRDS(model$spc_nps_rf,'data/spc_nps_rf.rds')
 }else{
  cat(readRDS('data/log.rds')[['spc_nps_rf']])
 }






################################################################################
#```{ Conduct selected NPS-PC LDA by parallel computing }
if(run_heavy_computation){
  cat('Conduct selected NPS-PC LDA by parallel computing\n')
  cat('Started:',as.character(now()),'\n')
  pb=startpb(0,3)
  on.exit(closepb(pb))
  setpb(pb,0)
  cl=makePSOCKcluster(floor(detectCores()/2))
  registerDoParallel(cl)
  clusterEvalQ(cl,{
    library('parallel')
    library('doParallel')
    library('tidyverse')
    library('pbapply')
    library('caret')
    library('lda')
   })
  
  setpb(pb,1)
  suppressWarnings(set.seed(33,sample.kind=sample.kind))
  model$spc_nps_lda=
    train(
      Class~.
      ,data=
        set$training %>%
        select(Class,cens_p) %>%
        cbind(set$training_spc_nps)
      ,method='lda'
      ,metric='ROC'
      ,trControl=trainControl(
        'cv'
        ,number=6
        ,summaryFunction=twoClassSummary
        ,classProbs=T
        ,savePredictions=T
        ,allowParallel=T
      )
      ,tuneLength=10
    )
  
  setpb(pb,2)
  model$spc_nps_lda=
    train(
      Class~.
      ,data=
        set$training %>%
        select(Class,cens_p) %>%
        cbind(set$training_spc_nps)
      ,method='lda'
      ,metric='ROC'
      ,trControl=
        trainControl(
          method='boot'
          ,number=30
          ,summaryFunction=twoClassSummary
          ,classProbs=T
          ,savePredictions=T
          ,allowParallel=T
        )
      ,tuneGrid=data.frame(parameter=model$spc_nps_lda$bestTune$parameter)
    )
  
  stopCluster(cl)
  registerDoSEQ()
  rm(cl)
  gc()
  setpb(pb,3)
  cat('\nEnd:',as.character(now()))
  rm(pb)
  saveRDS(model$spc_nps_lda,'data/spc_nps_lda.rds')
 }else{
  cat(readRDS('data/log.rds')[['spc_nps_lda']])
 }






################################################################################
#```{ Conduct selected NPS-PC GBM by parallel computing }
if(run_heavy_computation){
  cat('Conduct selected NPS-PC GBM by parallel computing\n')
  cat('Started:',as.character(now()),'\n')
  pb=startpb(0,3)
  on.exit(closepb(pb))
  setpb(pb,0)
  cl=makePSOCKcluster(floor(detectCores()/2))
  registerDoParallel(cl)
  clusterEvalQ(cl,{
    library('parallel')
    library('doParallel')
    library('tidyverse')
    library('pbapply')
    library('caret')
    library('gbm')
   })
  
  setpb(pb,1)
  suppressWarnings(set.seed(33,sample.kind=sample.kind))
  model$spc_nps_gbm=
    train(
      Class~.
      ,data=
        set$training %>%
        select(Class,cens_p) %>%
        cbind(set$training_spc_nps)
      ,method='gbm'
      ,metric='ROC'
      ,trControl=
        trainControl(
          'cv'
          ,number=6
          ,summaryFunction=twoClassSummary
          ,classProbs=T
          ,savePredictions=T
          ,allowParallel=T
        )
      ,tuneLength=10
    )
  
  setpb(pb,2)
  model$spc_nps_gbm=
    train(
      Class~.
      ,data=
        set$training %>%
        select(Class,cens_p) %>%
        cbind(set$training_spc_nps)
      ,method='gbm'
      ,metric='ROC'
      ,trControl=
        trainControl(
          method='boot'
          ,number=30
          ,summaryFunction=twoClassSummary
          ,classProbs=T
          ,savePredictions=T
          ,allowParallel=T
        )
      ,tuneGrid=data.frame(
        n.trees=model$spc_nps_gbm$bestTune$n.trees
        ,interaction.depth=model$spc_nps_gbm$bestTune$interaction.depth
        ,shrinkage=model$spc_nps_gbm$bestTune$shrinkage
        ,n.minobsinnode=model$spc_nps_gbm$bestTune$n.minobsinnode
      )
    )
  
  stopCluster(cl)
  registerDoSEQ()
  rm(cl)
  gc()
  setpb(pb,3)
  cat('\nEnd:',as.character(now()))
  rm(pb)
  saveRDS(model$spc_nps_gbm,'data/spc_nps_gbm.rds')
 }else{
  cat(readRDS('data/log.rds')[['spc_nps_gbm']])
 }





################################################################################
#```{ Create empty list to save reported data }
reported_data=list()






################################################################################
#```{ Create calibration plot before calibration }
if(run_heavy_computation){
  model %>%
    lapply(X=names(.) %>% .[!.%in%c('pca','pca_nps')],Y=.,function(X,Y){
      mutate(
        Y[[X]]$pred
        ,event=round(event,1)
        ,obs=as.integer(obs=='event')
      ) %>%
        group_by(event) %>%
        summarize(
          model=factor(X,names(Y) %>% .[!.%in%c('pca','pca_nps')])
          ,obs=mean(obs)
          ,lb=mean(obs)-qnorm(0.975)*sqrt(mean(obs)*(1-mean(obs))/n())
          ,ub=mean(obs)+qnorm(0.975)*sqrt(mean(obs)*(1-mean(obs))/n())
        )
     }) %>%
    do.call(rbind,.) %>%
    saveRDS('data/report_pre_calib.rds')
 }else{
  reported_data$pre_calib=readRDS('data/report_pre_calib.rds')
 }






################################################################################
#```{ Compute probability distribution before calibration }
if(run_heavy_computation){
  model %>%
    lapply(X=names(.) %>% .[!.%in%c('pca','pca_nps')],Y=.,function(X,Y){
      mutate(
        Y[[X]]$pred
        ,event=round(event,1)
        ,obs=as.integer(obs=='event')
      ) %>%
        group_by(event) %>%
        summarize(
          model=factor(X,names(Y) %>% .[!.%in%c('pca','pca_nps')])
          ,n=n()
        )
     }) %>%
    do.call(rbind,.) %>%
    saveRDS('data/report_pre_caldist.rds')
 }else{
  reported_data$pre_caldist=readRDS('data/report_pre_caldist.rds')
 }






################################################################################
#```{ Compute calibration intercept and slope before calibration }
if(run_heavy_computation){
  model %>%
    lapply(X=names(.) %>% .[!.%in%c('pca','pca_nps')],Y=.,function(X,Y){
      mutate(
        Y[[X]]$pred
        ,event=round(event,2)
        ,obs=as.integer(obs=='event')
      ) %>%
        group_by(event) %>%
        summarize(
          model=factor(X,names(Y) %>% .[!.%in%c('pca','pca_nps')])
          ,obs=mean(obs)
        )
     }) %>%
    do.call(rbind,.) %>%
    group_by(model) %>%
    do(tidy(lm(obs~event,data=.))) %>%
    mutate(term=str_remove_all(term,'\\(|\\)') %>% str_to_lower()) %>%
    pivot_wider(
      model
      ,names_from=c('term','term')
      ,values_from=c('estimate','std.error')
    ) %>%
    select(model,estimate_intercept,std.error_intercept,everything()) %>%
    mutate_at(colnames(.) %>% .[.!='model'],function(x)round(x,2)) %>%
    saveRDS('data/report_pre_calislope.rds')
 }else{
  reported_data$pre_calislope=readRDS('data/report_pre_calislope.rds')
 }






################################################################################
#```{ Create ROC curve among models before calibration }
if(run_heavy_computation){
  model %>%
    lapply(X=names(.) %>% .[!.%in%c('pca','pca_nps')],Y=.,function(X,Y){
      Y[[X]]$results %>%
        mutate(ROC_lb=ROC-qnorm(0.975)*ROCSD,ROC_ub=ROC+qnorm(0.975)*ROCSD) %>%
        mutate(model=X) %>%
        select(model,ROC,ROC_lb,ROC_ub)
     }) %>%
    do.call(rbind,.) %>%
    mutate(model=reorder(model,ROC)) %>%
    saveRDS('data/report_pre_roc.rds')
 }else{
  reported_data$pre_roc=readRDS('data/report_pre_roc.rds')
 }






################################################################################
#```{ Compute AUC of ROC among models before calibration }
if(run_heavy_computation){
  model %>%
    lapply(X=names(.) %>% .[!.%in%c('pca','pca_nps')],Y=.,function(X,Y){
      Y[[X]]$results %>%
        mutate(ROC_lb=ROC-qnorm(0.975)*ROCSD,ROC_ub=ROC+qnorm(0.975)*ROCSD) %>%
        mutate(model=X) %>%
        select(model,ROC,ROC_lb,ROC_ub)
     }) %>%
    do.call(rbind,.) %>%
    saveRDS('data/report_pre_auroc.rds')
 }else{
  reported_data$pre_auroc=readRDS('data/report_pre_auroc.rds')
 }






################################################################################
#```{ Get predicted probability of training data }
if(run_heavy_computation){
  cat('Get predicted probability of training data\n')
  cat('Started:',as.character(now()),'\n')
  
  set$training_calib=
    model %>%
    pblapply(X=names(.) %>% .[!.%in%c('pca','pca_nps')],Y=.,function(X,Y){
      Y[[X]] %>%
        predict(Y[[X]]$trainingData,type='prob') %>%
        cbind(select(Y[[X]]$trainingData,.outcome)) %>%
        select(event,.outcome) %>%
        setNames(c('event','obs')) %>%
        mutate(obs=as.integer(obs=='event'))
     }) %>%
    setNames(names(model) %>% .[!.%in%c('pca','pca_nps')])
  
  cat('End:',as.character(now()))
  saveRDS(set$training_calib,'data/training_calib.rds')
 }else{
  cat(readRDS('data/log.rds')[['training_calib']])
 }






################################################################################
#```{ Conduct calibration using predicted probability ... }
if(run_heavy_computation){
  cat('Conduct calibration using predicted probability',append=T)
  cat('and outcome by logistic regression\n')
  cat('Started:',as.character(now()),'\n')
  
  calib_model=
    model %>%
    pblapply(X=names(.) %>% .[!.%in%c('pca','pca_nps')]
             ,Y=.
             ,Z=set$training_calib
             ,function(X,Y,Z){
               suppressWarnings(set.seed(33,sample.kind=sample.kind))
               suppressWarnings(train(
                 obs~event
                 ,data=
                   Z[[X]] %>%
                   mutate(
                     obs=factor(
                       ifelse(obs==1,'event','non_event')
                       ,c('event','non_event')
                     )
                   )
                 ,method='glm'
                 ,metric='ROC'
                 ,family=binomial(link='logit')
                 ,trControl=
                   trainControl(
                     method='boot'
                     ,number=30
                     ,summaryFunction=twoClassSummary
                     ,classProbs=T
                     ,savePredictions=T
                   )
               ))
              }) %>%
    setNames(names(model) %>% .[!.%in%c('pca','pca_nps')])
  
  cat('End:',as.character(now()))
  saveRDS(calib_model,'data/calib_model.rds')
 }else{
  cat(readRDS('data/log.rds')[['calib_model']])
 }






################################################################################
#```{ Create calibration plot after calibration }
if(run_heavy_computation){
  calib_model %>%
    lapply(X=names(.),Y=.,function(X,Y){
      mutate(
        Y[[X]]$pred
        ,event=round(event,1)
        ,obs=as.integer(obs=='event')
      ) %>%
        group_by(event) %>%
        summarize(
          model=factor(X,names(Y))
          ,obs=mean(obs)
          ,lb=mean(obs)-qnorm(0.975)*sqrt(mean(obs)*(1-mean(obs))/n())
          ,ub=mean(obs)+qnorm(0.975)*sqrt(mean(obs)*(1-mean(obs))/n())
        )
     }) %>%
    do.call(rbind,.) %>%
    saveRDS('data/report_post_calib.rds')
 }else{
  reported_data$post_calib=readRDS('data/report_post_calib.rds')
 }






################################################################################
#```{ Compute probability distribution after calibration }
if(run_heavy_computation){
  calib_model %>%
    lapply(X=names(.),Y=.,function(X,Y){
      mutate(
        Y[[X]]$pred
        ,event=round(event,1)
        ,obs=as.integer(obs=='event')
      ) %>%
        group_by(event) %>%
        summarize(
          model=factor(X,names(Y))
          ,n=n()
        )
     }) %>%
    do.call(rbind,.) %>%
    saveRDS('data/report_post_caldist.rds')
 }else{
  reported_data$post_caldist=readRDS('data/report_post_caldist.rds')
 }






################################################################################
#```{ Compute calibration intercept and slope after calibration }
if(run_heavy_computation){
  calib_model %>%
    lapply(X=names(.),Y=.,function(X,Y){
      mutate(
        Y[[X]]$pred
        ,event=round(event,1)
        ,obs=as.integer(obs=='event')
      ) %>%
        group_by(event) %>%
        summarize(
          model=factor(X,names(Y))
          ,obs=mean(obs)
        )
     }) %>%
    do.call(rbind,.) %>%
    group_by(model) %>%
    do(tidy(lm(obs~event,data=.))) %>%
    mutate(term=str_remove_all(term,'\\(|\\)') %>% str_to_lower()) %>%
    pivot_wider(
      model
      ,names_from=c('term','term')
      ,values_from=c('estimate','std.error')
    ) %>%
    select(model,estimate_intercept,std.error_intercept,everything()) %>%
    mutate_at(colnames(.) %>% .[.!='model'],function(x)round(x,2)) %>%
    saveRDS('data/report_post_calislope.rds')
 }else{
  reported_data$post_calislope=readRDS('data/report_post_calislope.rds')
 }






################################################################################
#```{ Create ROC curve among models after calibration }
if(run_heavy_computation){
  calib_model %>%
    lapply(X=names(.) %>% .[!.%in%c('pca','pca_nps')],Y=.,function(X,Y){
      Y[[X]]$results %>%
        mutate(ROC_lb=ROC-qnorm(0.975)*ROCSD,ROC_ub=ROC+qnorm(0.975)*ROCSD) %>%
        mutate(model=X) %>%
        select(model,ROC,ROC_lb,ROC_ub)
     }) %>%
    do.call(rbind,.) %>%
    mutate(model=reorder(model,ROC)) %>%
    saveRDS('data/report_post_roc.rds')
 }else{
  reported_data$post_roc=readRDS('data/report_post_roc.rds')
 }






################################################################################
#```{ Compute AUC of ROC among models after calibration }
if(run_heavy_computation){
  calib_model %>%
    lapply(X=names(.) %>% .[!.%in%c('pca','pca_nps')],Y=.,function(X,Y){
      Y[[X]]$results %>%
        mutate(ROC_lb=ROC-qnorm(0.975)*ROCSD,ROC_ub=ROC+qnorm(0.975)*ROCSD) %>%
        mutate(model=X) %>%
        select(model,ROC,ROC_lb,ROC_ub)
     }) %>%
    do.call(rbind,.) %>%
    saveRDS('data/report_post_auroc.rds')
 }else{
  reported_data$post_auroc=readRDS('data/report_post_auroc.rds')
 }






################################################################################
#```{ Construct a class-balanced testing set by provider ... }
if(run_heavy_computation){
  mh_provider=readRDS('data/mh_provider.rds')
  construct_testing_set=function(testing_idx_db){
    
    # Use provider medical history
    temp=
      mh_provider %>%
      select(-subject_id,-day_to_event,-est_strata_id) %>%
      t() %>%
      t()
    
    temp=
      (temp>0) %>%
      `storage.mode<-`('integer')
    
    # Make the confirmed causal predictors defined by provider medical histories
    temp2=
      lapply(
        X=seq(sum(
          str_remove_all(dag$measure_nodes$label,'\\*') %in%
            colnames(dag$sig)[dag$sig[2,]==1]
        ))
        ,Y=
          dag$measure_nodes %>%
          filter(
            str_remove_all(label,'\\*') %in%
              colnames(dag$sig)[dag$sig[2,]==1]
          ) %>%
          pull(name)
        ,Z=dag$measure_nodes %>%
          filter(
            str_remove_all(label,'\\*') %in%
              colnames(dag$sig)[dag$sig[2,]==1]
          ) %>%
          pull(label)
        ,K=temp
        ,function(X,Y,Z,K){
          K=K %>%
            .[,str_detect(str_to_upper(colnames(.)),Y[X]),drop=F] %>%
            rowSums() %>%
            matrix() %>%
            `dimnames<-`(list(
              rownames(K)
              ,paste0('causal_',str_remove_all(Z[X],'\\*'))
            ))
          
          K=(K>0) %>%
            `storage.mode<-`('integer')
          
         }) %>%
      do.call(cbind,.)
    
    # Combine medical history and the causal predictors,
    # then convert frequency into probability or proportion of day being
    # encountered up to 2 years before the event
    temp=
      cbind(temp2,temp) %>%
      sweep(1,((365*2-mh_provider$day_to_event)/(365*2)),'*')
    rm(temp2)
    
    # Join the results with other attributes for data partition
    temp=
      mh_provider %>%
      select(subject_id,day_to_event,est_strata_id) %>%
      cbind(as.data.frame(temp))
    
    # Compute mean and standard deviation of age based on training set
    training_age_sum=
      readRDS('data/target_population.rds') %>%
      select(subject_id,day_to_event,outcome,age) %>%
      .[!duplicated(.),] %>%
      filter(subject_id %in% set$intv$subject_id) %>%
      summarize(age_m=mean(age),age_s=sd(age))
    
    # Get non-medical history predictors
    testing_set=
      
      ## Numerical predictors (age), standardize, and normalize into 0 to 1
      ## using the training set
      readRDS('data/target_population.rds') %>%
      select(subject_id,day_to_event,outcome,age) %>%
      .[!duplicated(.),] %>%
      filter(subject_id %in% testing_idx_db$subject_id) %>%
      cbind(training_age_sum) %>%
      mutate(age=(qnorm(0.975)+(age-age_m)/age_s)/(2*qnorm(0.975))) %>%
      
      ## Categorical predictors and getting these cleaned
      left_join(
        readRDS('data/target_population.rds') %>%
          select(
            subject_id
            ,marital_status
            ,insurance_class
            ,occupation_segment
            ,healthcare_class
          ) %>%
          .[!duplicated(.),] %>%
          group_by(subject_id) %>%
          summarize_all(function(x){
            x[!duplicated(x)] %>%
              .[which.max(
                sapply(x[!duplicated(x)],function(x2)sum(x2==x,na.rm=T))
              )]
           }) %>%
          lapply(X=seq(ncol(.)),Y=.,function(X,Y){
            if(colnames(Y)[X] %in% c('marital_status'
                                     ,'insurance_class'
                                     ,'occupation_segment'
                                     ,'healthcare_class')){
              Y[,c('subject_id',colnames(Y)[X]),drop=F] %>%
                rename_at(colnames(Y)[X],function(x)'key') %>%
                mutate(
                  key=str_replace_all(key,'[:punct:]|[:space:]','_')
                  ,value=1
                ) %>%
                mutate_at('key',function(x)paste0(colnames(Y)[X],'.',x)) %>%
                spread(key,value) %>%
                select(-subject_id)
             }else{
              Y[,colnames(Y)[X],drop=F]
             }
           }) %>%
          do.call(cbind,.) %>%
          mutate_at(
            colnames(.) %>% .[.!='subject_id']
            ,function(x)as.integer(!is.na(x))
          )
        ,by='subject_id'
      ) %>%
      
      ## Subset internal validation set and add censoring probability
      right_join(
        temp %>%
          filter(subject_id %in% testing_idx_db$subject_id) %>%
          group_by(subject_id) %>%
          mutate(cens_p=(365*2-max(day_to_event))/(365*2)) %>%
          ungroup()
        ,by=c('subject_id','day_to_event')
      ) %>%
      
      # Convert outcome as 0 and 1
      mutate(outcome=as.integer(outcome=='event')) %>%
      
      # Pre-finishing touch
      mutate(causal_A20=age) %>%
      select(subject_id,day_to_event,est_strata_id,outcome,cens_p,everything())
    rm(temp)
    
    # Upsample minority to balance the outcome
    suppressWarnings(set.seed(33,sample.kind=sample.kind))
    testing_set %>%
      upSample(factor(.$outcome)) %>%
      mutate(
        Class=
          ifelse(Class==0,'non_event','event') %>%
          factor(c('event','non_event'))
      ) %>%
      select(
        subject_id
        ,day_to_event
        ,outcome
        ,Class
        ,age_m
        ,age_s
        ,cens_p
        ,everything()
      )
   }
  
  set$testing_geo=construct_testing_set(set$extv_geo)
  set$testing_tem=construct_testing_set(set$extv_tem)
  set$testing_bgt=construct_testing_set(set$extv_bgt)
  set$testing=
    set[c('testing_geo','testing_tem','testing_bgt')] %>%
    do.call(rbind,.)
  pred=list()
  mod_eval=list()
  
  rm(mh_provider)
 }






################################################################################
#```{ ... testing ... ridge regression with causal predictors }
if(run_heavy_computation){
  cat('Conduct prediction and compute ROC on testing sets',append=T)
  cat('for the ridge regression with causal predictors\n')
  cat('Started:',as.character(now()),'\n')
  
  pred$ridge_calib=
    lapply(X=1:5
           ,Y=set[c(
             'training'
             ,'testing_geo'
             ,'testing_tem'
             ,'testing_bgt'
             ,'testing'
           )]
           ,Z=model
           ,K=calib_model
           ,L='ridge'
           ,function(X,Y,Z,K,L){
             M=Y[[X]] %>%
               .[,c('Class','cens_p',colnames(.) %>% .[str_detect(.,'causal')])]
             
             suppressWarnings(set.seed(33,sample.kind=sample.kind))
             data.frame(
               event=predict(Z[[L]],newdata=M,type='prob')$event
               ,obs=M$Class
             ) %>%
               mutate(calib=predict(K[[L]],newdata=.,type='prob'))
            }) %>%
    setNames(c(
      'training'
      ,'testing_geo'
      ,'testing_tem'
      ,'testing_bgt'
      ,'testing'
    ))
  
  mod_eval$ridge_calib=
    pred$ridge_calib %>%
    pblapply(X=names(.),Y=.,function(X,Y){
      data.frame(Y[[X]]$calib,Y[[X]]$obs)
     }) %>%
    pblapply(evalm) %>%
    setNames(c(
      'training'
      ,'testing_geo'
      ,'testing_tem'
      ,'testing_bgt'
      ,'testing'
    ))
  
  cat('End:',as.character(now()))
  saveRDS(pred$ridge_calib,'data/pred_ridge_calib.rds')
  saveRDS(mod_eval$ridge_calib,'data/mod_eval_ridge_calib.rds')
 }else{
  cat(readRDS('data/log.rds')[['pred_eval_ridge_calib']])
 }






################################################################################
#```{ ... testing ... NPS-PC elastic net regression model }
if(run_heavy_computation){
  cat('Conduct prediction and compute ROC on testing sets',append=T)
  cat('for the NPS-PC elastic net regression model\n')
  cat('Started:',as.character(now()),'\n')
  
  pred$pc_nps_elnet_calib=
    lapply(X=1:5
           ,Y=set[c(
             'training'
             ,'testing_geo'
             ,'testing_tem'
             ,'testing_bgt'
             ,'testing'
           )]
           ,Z=model
           ,K=calib_model
           ,L='pc_nps_elnet'
           ,function(X,Y,Z,K,L){
             M=select(Y[[X]],Class,cens_p) %>%
               cbind(
                 pc_converter(
                   Y[[X]]
                   ,model$pca_nps
                   ,npc=as.numeric(str_remove_all(
                     colnames(model$pc_nps_elnet$trainingData)[-1:-2]
                     ,'PC'
                   ))
                 )
               )
             
             suppressWarnings(set.seed(33,sample.kind=sample.kind))
             data.frame(
               event=predict(Z[[L]],newdata=M,type='prob')$event
               ,obs=M$Class
             ) %>%
               mutate(calib=predict(K[[L]],newdata=.,type='prob'))
            }) %>%
    setNames(c(
      'training'
      ,'testing_geo'
      ,'testing_tem'
      ,'testing_bgt'
      ,'testing'
    ))
  
  mod_eval$pc_nps_elnet_calib=
    pred$pc_nps_elnet_calib %>%
    pblapply(X=names(.),Y=.,function(X,Y){
      data.frame(Y[[X]]$calib,Y[[X]]$obs)
     }) %>%
    pblapply(evalm) %>%
    setNames(c(
      'training'
      ,'testing_geo'
      ,'testing_tem'
      ,'testing_bgt'
      ,'testing'
    ))
  
  cat('End:',as.character(now()))
  saveRDS(pred$pc_nps_elnet_calib,'data/pred_pc_nps_elnet_calib.rds')
  saveRDS(mod_eval$pc_nps_elnet_calib,'data/mod_eval_pc_nps_elnet_calib.rds')
 }else{
  cat(readRDS('data/log.rds')[['pred_eval_pc_nps_elnet_calib']])
 }






################################################################################
#```{ ... testing ... random forest model }
if(run_heavy_computation){
  cat('Conduct prediction and compute ROC on testing sets',append=T)
  cat('for the random forest model\n')
  cat('Started:',as.character(now()),'\n')
  
  pred$rf_calib=
    lapply(X=1:5
           ,Y=set[c(
             'training'
             ,'testing_geo'
             ,'testing_tem'
             ,'testing_bgt'
             ,'testing'
           )]
           ,Z=model
           ,K=calib_model
           ,L='spc_nps_rf'
           ,function(X,Y,Z,K,L){
             M=select(Y[[X]],Class,cens_p) %>%
               cbind(
                 pc_converter(
                   Y[[X]]
                   ,model$pca_nps
                   ,npc=sort(set$selected_pc$idx)
                 )
               )
             
             suppressWarnings(set.seed(33,sample.kind=sample.kind))
             data.frame(
               event=predict(Z[[L]],newdata=M,type='prob')$event
               ,obs=M$Class
             ) %>%
               mutate(calib=predict(K[[L]],newdata=.,type='prob'))
            }) %>%
    setNames(c(
      'training'
      ,'testing_geo'
      ,'testing_tem'
      ,'testing_bgt'
      ,'testing'
    ))
  
  mod_eval$rf_calib=
    pred$rf_calib %>%
    pblapply(X=names(.),Y=.,function(X,Y)data.frame(Y[[X]]$calib,Y[[X]]$obs)) %>%
    pblapply(evalm) %>%
    setNames(c(
      'training'
      ,'testing_geo'
      ,'testing_tem'
      ,'testing_bgt'
      ,'testing'
    ))
  
  cat('End:',as.character(now()))
  saveRDS(pred$rf_calib,'data/pred_rf_calib.rds')
  saveRDS(mod_eval$rf_calib,'data/mod_eval_rf_calib.rds')
 }else{
  cat(readRDS('data/log.rds')[['pred_eval_rf_calib']])
 }






################################################################################
#```{ ... testing ... LDA model }
if(run_heavy_computation){
  cat('Conduct prediction and compute ROC on testing sets for the LDA model\n')
  cat('Started:',as.character(now()),'\n')
  
  pred$lda_calib=
    lapply(X=1:5
           ,Y=set[c(
             'training'
             ,'testing_geo'
             ,'testing_tem'
             ,'testing_bgt'
             ,'testing'
           )]
           ,Z=model
           ,K=calib_model
           ,L='spc_nps_lda'
           ,function(X,Y,Z,K,L){
             M=select(Y[[X]],Class,cens_p) %>%
               cbind(
                 pc_converter(
                   Y[[X]]
                   ,model$pca_nps
                   ,npc=sort(set$selected_pc$idx)
                 )
               )
             
             suppressWarnings(set.seed(33,sample.kind=sample.kind))
             data.frame(
               event=predict(Z[[L]],newdata=M,type='prob')$event
               ,obs=M$Class
             ) %>%
               mutate(calib=predict(K[[L]],newdata=.,type='prob'))
            }) %>%
    setNames(c(
      'training'
      ,'testing_geo'
      ,'testing_tem'
      ,'testing_bgt'
      ,'testing'
    ))
  
  mod_eval$lda_calib=
    pred$lda_calib %>%
    pblapply(X=names(.),Y=.,function(X,Y){
      data.frame(Y[[X]]$calib,Y[[X]]$obs)
     }) %>%
    pblapply(evalm) %>%
    setNames(c(
      'training'
      ,'testing_geo'
      ,'testing_tem'
      ,'testing_bgt'
      ,'testing'
    ))
  
  cat('End:',as.character(now()))
  saveRDS(pred$lda_calib,'data/pred_lda_calib.rds')
  saveRDS(mod_eval$lda_calib,'data/mod_eval_lda_calib.rds')
 }else{
  cat(readRDS('data/log.rds')[['pred_eval_lda_calib']])
 }






################################################################################
#```{ ... testing ... PC-GBM model }
if(run_heavy_computation){
  cat('Conduct prediction and compute ROC on testing sets',append=T)
  cat('for the best model or the calibrated GBM model\n')
  cat('Started:',as.character(now()),'\n')
  
  pred$gbm_calib=
    lapply(X=1:5
           ,Y=set[c(
             'training'
             ,'testing_geo'
             ,'testing_tem'
             ,'testing_bgt'
             ,'testing'
           )]
           ,Z=model
           ,K=calib_model
           ,L='spc_nps_gbm'
           ,function(X,Y,Z,K,L){
             M=select(Y[[X]],Class,cens_p) %>%
               cbind(
                 pc_converter(
                   Y[[X]]
                   ,model$pca_nps
                   ,npc=sort(set$selected_pc$idx)
                 )
               )
             
             suppressWarnings(set.seed(33,sample.kind=sample.kind))
             data.frame(
               event=predict(Z[[L]],newdata=M,type='prob')$event
               ,obs=M$Class
             ) %>%
               mutate(calib=predict(K[[L]],newdata=.,type='prob'))
            }) %>%
    setNames(c(
      'training'
      ,'testing_geo'
      ,'testing_tem'
      ,'testing_bgt'
      ,'testing'
    ))
  
  mod_eval$gbm_calib=
    pred$gbm_calib %>%
    pblapply(X=names(.),Y=.,function(X,Y){
      data.frame(Y[[X]]$calib,Y[[X]]$obs)
     }) %>%
    pblapply(evalm) %>%
    setNames(c(
      'training'
      ,'testing_geo'
      ,'testing_tem'
      ,'testing_bgt'
      ,'testing'
    ))
  
  cat('End:',as.character(now()))
  saveRDS(pred$gbm_calib,'data/pred_gbm_calib.rds')
  saveRDS(mod_eval$gbm_calib,'data/mod_eval_gbm_calib.rds')
 }else{
  cat(readRDS('data/log.rds')[['pred_eval_gbm_calib']])
 }






################################################################################
#```{ Compare testing AUROC among models }
if(run_heavy_computation){
  mod_eval %>%
    .[!str_detect(names(.),'pte')] %>%
    lapply(X=names(.),Y=.,function(X,Y){
      Z=Y[[X]] %>%
        lapply(function(x)x$optres$Group1['AUC-ROC',]) %>%
        do.call(rbind,.) %>%
        rownames_to_column(var='set') %>%
        mutate(model=X)
      suppressWarnings(separate(Z,CI,c('lb','ub'),sep='-'))
     }) %>%
    do.call(rbind,.) %>%
    mutate(
      lb=ifelse(lb=='NA',Score,lb)
      ,ub=ifelse(ub=='NA',Score,ub)
    ) %>%
    mutate_at(colnames(.) %>% .[!.%in%c('model','set')],function(x)as.numeric(x)) %>%
    saveRDS('data/report_calib_test.rds')
 }else{
  reported_data$calib_test=readRDS('data/report_calib_test.rds')
 }







################################################################################
#```{ Get rotated PCs as the PC weights }
if(run_heavy_computation){
  rotated_pc=
    
    # Subset the top PCs and stack all versions of PCs
    model$pca_nps %>%
    lapply(X=seq(length(.)),Y=.,function(X,Y){
      as.data.frame(Y[[X]]$prcomp$rotation[,sort(set$selected_pc$idx)]) %>%
        rownames_to_column(var='predictor')
     }) %>%
    do.call(rbind,.) %>%
    
    # Group by predictor and average the PC weights over all versions
    group_by(predictor) %>%
    summarize_all(function(x)mean(x,na.rm=T)) %>%
    ungroup() %>%
    column_to_rownames(var='predictor')
  
  saveRDS(rotated_pc,'data/rotated_pc.rds')
 }else{
  rotated_pc=readRDS('data/rotated_pc.rds')
 }






################################################################################
#```{ Get variable importance of PC-GBM }
if(run_heavy_computation){
  varimp_gbm=
    model$spc_nps_gbm %>%
    varImp() %>%
    .$importance %>%
    rownames_to_column(var='feature') %>%
    rename(importance=Overall)
  
  saveRDS(varimp_gbm,'data/varimp_gbm.rds')
 }else{
  varimp_gbm=readRDS('data/varimp_gbm.rds')
 }
