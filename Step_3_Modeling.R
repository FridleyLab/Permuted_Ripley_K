#This script conducts the cross validation and constructs the final survival models for each marker. This script was used for the intrarumoral ROIs, but
#it can easily be revised to analyze the tumor compartments of the TMAs.

library(spatialTIME)
library(tidyverse)
library(survival)
library(survminer)

#final_roi_data contains all clinical, summary, and spatial data for the intratumoral ROIs.

#Rename marker
final_roi_data = final_roi_data %>% 
  mutate(Marker = case_when(Marker == "CD11b..CD15." ~ "CD11b+ CD15+ Positive Cells",
                            Marker == "CD11b..Opal.620..Positive" ~ "CD11b (Opal 620) Positive Cells",
                            Marker == "CD15..Opal.520..Positive" ~ "CD15 (Opal 520) Positive Cells",
                            Marker == "CD3..CD8." ~ "CD3+ CD8+ Positive Cells",
                            Marker == "CD3..FOXP3." ~ "CD3+ FOXP3+ Positive Cells",
                            Marker == "CD3..Opal.650..Positive" ~ "CD3 (Opal 650) Positive Cells",
                            Marker == "CD8..Opal.570..Positive" ~ "CD8 (Opal 570) Positive Cells",
                            Marker == "FOXP3..Opal.540..Positive" ~ "FOXP3 (Opal 540) Positive Cells")) %>%
  rename(
    "CD11b+ CD15+ Positive Cells" = cd11bplus_cd15plus_cells,
    "CD11b (Opal 620) Positive Cells" = cd11b_opal_620_positive_cells,
    "CD15 (Opal 520) Positive Cells" = cd15_opal_520_positive_cells,
    "CD3+ CD8+ Positive Cells" = cd3plus_cd8plus_cells,
    "CD3+ FOXP3+ Positive Cells" = cd3plus_foxp3plus_cells,
    "CD3 (Opal 650) Positive Cells" = cd3_opal_650_positive_cells,
    "CD8 (Opal 570) Positive Cells" = cd8_opal_570_positive_cells,
    "FOXP3 (Opal 540) Positive Cells" = foxp3_opal_540_positive_cells,
    "% CD11b+ CD15+ Positive Cells" = percent_cd11bplus_cd15plus_positive_cells,
    "% CD11b (Opal 620) Positive Cells" = percent_cd11b_opal_620_positive_cells,
    "% CD15 (Opal 520) Positive Cells" = percent_cd15_opal_520_positive_cells,
    "% CD3+ CD8+ Positive Cells" = percent_cd3plus_cd8plus_positive_cells,
    "% CD3+ FOXP3+ Positive Cells" = percent_cd3plus_foxp3plus_positive_cells,
    "% CD3 (Opal 650) Positive Cells" = percent_cd3_opal_650_positive_cells,
    "% CD8 (Opal 570) Positive Cells" = percent_cd8_opal_570_positive_cells,
    "% FOXP3 (Opal 540) Positive Cells" = percent_foxp3_opal_540_positive_cells
    )


#Extract marker names to loop through
markers = unique(final_roi_data$Marker)
best = data.frame()

for(marker in markers){
  print(marker)
#subset to a particular marker and standardize the name for the sake of looping
  data_by_marker = final_roi_data %>% 
    filter(Marker == marker) %>%
    rename(Cell_Percent = paste('%', marker))
  
#Remove samples with 0 cells since they should not contribute to the spatial component
  no_zeros = data_by_marker %>% filter(!is.na(`Degree of Clustering Permutation`)) 

#Enumerate the combinations of potential cuts for both abundance and degree of clustering (169 combinations)
  grid = expand_grid(abundance = quantile(no_zeros$Cell_Percent, seq(0.2,0.8,0.05)),
                     clustering = quantile(no_zeros$`Degree of Clustering Permutation`, 
                                           seq(0.2,0.8,0.05))) %>%
    mutate(lrt_log_lik_full = NA,
           lrt_log_lik_clinical = NA,
           lrt = NA, big_group = NA)
  
  for(i in 1:10){
    for(b in 1:nrow(grid)){
	#leave out one fold and classify each sample based on abundance and clustering
      tmp = data_by_marker %>% 
        filter(fold != i) %>%
        mutate(five_group_percent = case_when(is.na(`Degree of Clustering Permutation`) == TRUE~ 'None',
                                              Cell_Percent > grid$abundance[b] &
                                                `Degree of Clustering Permutation` > grid$clustering[b] ~ 'HH',
                                              Cell_Percent > grid$abundance[b] &
                                                `Degree of Clustering Permutation` <= grid$clustering[b] ~ 'HL',
                                              Cell_Percent <= grid$abundance[b] &
                                                `Degree of Clustering Permutation` > grid$clustering[b] ~ 'LH',
                                              Cell_Percent <= grid$abundance[b] &
                                                `Degree of Clustering Permutation` <= grid$clustering[b] ~ 'LL'),
               five_group_percent = factor(five_group_percent, levels = c('None','HH',"HL", "LH", 'LL')),
               stage = factor(stage)) %>% 
        select(suid, vitalstatus_new, timelastfu_new, stage, refage, stage,
               five_group_percent) %>% data.frame(check.names = FALSE)
      
	    #full model
      model_lrt_1 = coxph(Surv(timelastfu_new, as.numeric(vitalstatus_new)) ~ refage + stage + five_group_percent + cluster(suid), 
                          data = tmp)
      #reduced model (spatial component is left out)
      model_lrt_2 = coxph(Surv(timelastfu_new, as.numeric(vitalstatus_new)) ~ refage + stage + cluster(suid), 
                          data = tmp)
      
      grid$lrt_log_lik_full[b] = summary(model_lrt_1)$loglik[2] #extract logliklihood value
      grid$lrt_log_lik_clinical[b] = summary(model_lrt_2)$loglik[2] #extract logliklihood value
      grid$lrt[b] = -2*log(grid$lrt_log_lik_full[b] - grid$lrt_log_lik_clinical[b]) #determine the benefit of adding the spatial component to the full model 
      grid$big_group[b] = all(table(tmp$five_group_percent)[1:4] >= 10)#we are going to constrain stratifications where each group has more than 10 samples
      grid$n_HH[b] = table(tmp$five_group_percent)['HH']
      grid$n_HL[b] = table(tmp$five_group_percent)['HL']
      grid$n_LH[b] = table(tmp$five_group_percent)['LH']
      grid$n_LL[b] = table(tmp$five_group_percent)['LL']
      grid$n_None[b] = table(tmp$five_group_percent)['None']
    }
    #record the results for each marker, fold that is left out and cutpoint
    best = plyr::rbind.fill(best,
                            grid %>% 
                              filter(!is.na(lrt)) %>%
                              arrange(lrt) %>%
                              head(1) %>%
                              mutate(left_out = i,
                                     Marker = marker))
  }
}

#Filter to the stratifications that had more than 10 samples and taking the median cutpoints for each marker
final_cut = best %>% group_by(Marker) %>%
  filter(big_group == TRUE) %>%
  select(abundance, clustering) %>%
  summarize_all(~median(.)) 

#Use the values from final_cut to construct the final models and make predicted curves
models_roi = lapply(setNames(final_cut$Marker,final_cut$Marker), function(i){
  #Use the cutpoints that lead to the model where the spatial component explained the most variation to build the final model.
  tmp = final_roi_data %>% 
    filter(Marker == i) %>%
    rename('Cell_Percent' = paste0('% ', i)) %>%
    mutate(five_group_percent = case_when(is.na(`Degree of Clustering Permutation`) == TRUE~ 'None',
                                          Cell_Percent > final_cut$abundance[final_cut$Marker == i] &
                                            `Degree of Clustering Permutation` > final_cut$clustering[final_cut$Marker == i] ~ 'HH',
                                          Cell_Percent > final_cut$abundance[final_cut$Marker == i] &
                                            `Degree of Clustering Permutation` <=final_cut$clustering[final_cut$Marker == i] ~ 'HL',
                                          Cell_Percent <= final_cut$abundance[final_cut$Marker == i]&
                                            `Degree of Clustering Permutation` >final_cut$clustering[final_cut$Marker == i] ~ 'LH',
                                          Cell_Percent <=final_cut$abundance[final_cut$Marker == i] &
                                            `Degree of Clustering Permutation` <= final_cut$clustering[final_cut$Marker == i] ~ 'LL'),
           five_group_percent = factor(five_group_percent, levels = c('None','HH',"HL", "LH", 'LL')),
           stage = factor(stage)) %>% 
    select(suid, vitalstatus_new, timelastfu_new, stage, refage, stage,
           five_group_percent) %>% data.frame(check.names = FALSE)
  
  model = coxph(Surv(timelastfu_new, as.numeric(vitalstatus_new)) ~ refage + stage + five_group_percent + cluster(suid), 
                data = tmp)
  model_2 = coxph(Surv(timelastfu_new, as.numeric(vitalstatus_new)) ~ refage + stage + cluster(suid), 
                data = tmp)
  
  plot = ggadjustedcurves(model, data = tmp, variable = 'five_group_percent', size = 0.75) +   
    ggtitle(i) + labs(col="Group") + xlab('Time (Days)') + 
    theme(legend.title = element_text(size=20),  legend.text = element_text(size = 18))

  lrt_log_lik_full = summary(model)$loglik[2]
  lrt_log_lik_clinical = summary(model_2)$loglik[2]
  lrt = 2*(lrt_log_lik_full - lrt_log_lik_clinical)
  print(c(i,lrt, 1-pchisq(lrt,4)))
  return(list(model = model, plot = plot, final_cut = final_cut, n = table(tmp$five_group_percent)))})

#Getting the model summary for each marker
final_models_output_roi = lapply(final_cut$Marker, function(a){
  model_fits = c(marker = a,
                 coefficients(summary(models_roi[[a]]$model))[4:7, 'exp(coef)'],
                 coefficients(summary(models_roi[[a]]$model))[4:7, 'Pr(>|z|)'],
                 models_roi[[a]]$n)
  return(model_fits)
}) %>% plyr::ldply()

colnames(final_models_output_roi)[6:9] = paste0(colnames(final_models_output_roi)[6:9], '_p.val')
colnames(final_models_output_roi)[10:14] = paste0(colnames(final_models_output_roi)[10:14], '_n')
