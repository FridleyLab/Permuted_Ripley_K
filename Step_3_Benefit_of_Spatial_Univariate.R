library(tidyverse)
library(survival)
library(survminer)
library(GGally)

#apply the final_cuts to the data
models_roi_spatabund = lapply(final_cut$Marker, function(i){
  tmp = final_roi_data %>% 
    filter(Marker == i) %>%
    rename(Cell_Percent = paste0('% ', i)) %>%
    mutate(five_group_percent = case_when(is.na(`Degree of Clustering Permutation`) == TRUE~ 'None',
                                          Cell_Percent > final_cut$abundance[final_cut$Marker == i] &
                                            `Degree of Clustering Permutation` > final_cut$clustering[final_cut$Marker == i] ~ 'HH',
                                          Cell_Percent > final_cut$abundance[final_cut$Marker == i] &
                                            `Degree of Clustering Permutation` <=final_cut$clustering[final_cut$Marker == i] ~ 'HL',
                                          Cell_Percent <= final_cut$abundance[final_cut$Marker == i]&
                                            `Degree of Clustering Permutation` >final_cut$clustering[final_cut$Marker == i] ~ 'LH',
                                          Cell_Percent <=final_cut$abundance[final_cut$Marker == i] &
                                            `Degree of Clustering Permutation` <= final_cut$clustering[final_cut$Marker == i] ~ 'LL'),
           #create categorical including spatial and abundance
           five_group_percent = factor(five_group_percent, levels = c('None','HH',"HL", "LH", 'LL')),
           #create categorical including only abundance information
           three_group_percent = case_when(Cell_Percent == 0 ~ 'None',
                                           Cell_Percent > final_cut$abundance[final_cut$Marker == i] ~ 'HIGH',
                                           Cell_Percent <= final_cut$abundance[final_cut$Marker == i]~ 'LOW'),
           three_group_percent = factor(three_group_percent, levels = c('None', 'HIGH', 'LOW')),
           stage = factor(stage)) %>% 
    select(suid, vitalstatus_new, timelastfu_new, stage, refage, stage,
           five_group_percent, three_group_percent) %>% data.frame(check.names = FALSE)
  
  #spatial and abundance in full model
  model = coxph(Surv(timelastfu_new, as.numeric(vitalstatus_new)) ~ refage + stage + five_group_percent + cluster(suid), 
                data = tmp)
  #only abundance in 'reduced'
  model_2 = coxph(Surv(timelastfu_new, as.numeric(vitalstatus_new)) ~ refage + stage + three_group_percent + cluster(suid), 
                  data = tmp)
  lrt_log_lik_full = summary(model)$loglik[2]
  lrt_log_lik_clinical = summary(model_2)$loglik[2]
  lrt = 2*(lrt_log_lik_full - lrt_log_lik_clinical)
  #calculate significance of spatial to model
  print(c(i,lrt, 1-pchisq(lrt,abs(summary(model)$logtest[2]-summary(model_2)$logtest[2]))))
  return(list(model = model,
              final_cut = final_cut, 
              full_model_n = table(tmp$five_group_percent),
              reduced_model_n = table(tmp$three_group_percent),
              lrt = lrt,
              p_value = c(Marker = i,
                          LRT = lrt, 
                          LRT_p_value = 1-pchisq(lrt,summary(model)$logtest[2]-summary(model_2)$logtest[2]))))})
