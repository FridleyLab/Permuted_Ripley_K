library(tidyverse)
library(survival)
library(survminer)
library(GGally)

#using bivariate final_cut
models_roi_spatabund = lapply(1:nrow(final_cut), function(i){
  tmp = final_data %>% 
    filter(anchor == final_cut$anchor[i],
           counted == final_cut$counted[i]) %>%
    mutate(Spatial_HL = case_when(Presence_anchor != 'P' | Presence_counted != 'P' | 
                                    is.nan(`Degree of Clustering Permutation`) ~ 'N',
                                  Presence_anchor == 'P' & Presence_counted == 'P' & 
                                    `Degree of Clustering Permutation` >= final_cut$cut[i] ~ 'H',
                                  Presence_anchor == 'P' & Presence_counted == 'P' & 
                                    `Degree of Clustering Permutation` < final_cut$cut[i] ~ 'L'),
           #Create 5 level variable
           final_group = paste0(Presence_anchor, Presence_counted, Spatial_HL),
           #create reduced 3 level variable
           final_group2 = paste0(Presence_anchor, Presence_counted))
  
  #model with spatial and abundance information
  model = coxph(Surv(timelastfu_new, as.numeric(vitalstatus_new)) ~ refage + stage + final_group + cluster(suid), 
                data = tmp)
  #model with only abundance information
  model_2 = coxph(Surv(timelastfu_new, as.numeric(vitalstatus_new)) ~ refage + stage + final_group2 + cluster(suid), 
                  data = tmp)
  lrt_log_lik_full = summary(model)$loglik[2]
  lrt_log_lik_clinical = summary(model_2)$loglik[2]
  lrt = 2*(lrt_log_lik_full - lrt_log_lik_clinical)
  #calculate contribution of spatial information to the model
  print(c(final_cut$anchor[i], final_cut$counted[i],lrt, 1-pchisq(lrt,summary(model)$logtest[2]-summary(model_2)$logtest[2])))
  return(list(model = model,
              final_cut = final_cut, 
              full_model_n = table(tmp$final_group),
              reduced_model_n = table(tmp$final_group2),
              lrt = lrt,
              p_value = c(Anchor = final_cut$anchor[i], 
                          Counted = final_cut$counted[i],
                          LRT = lrt, 
                          LRT_p_value = 1-pchisq(lrt,summary(model)$logtest[2]-summary(model_2)$logtest[2]))))})
