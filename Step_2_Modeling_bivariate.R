library(tidyverse)
library(survival)
library(survminer)
library(GGally)

#Assign fold values for each patient, then merge with the sample data
clin_summ = mif_tumor$clinical %>% 
  mutate(fold = sample(1:10, nrow(.), replace = TRUE)) %>%
  inner_join(mif_tumor$sample)
  

#Get all of the data into a single dataframe
#Using a radius of 30, subsetting to the AAS study, 
#and to the HGS histotype, rename anchors and counted 
#Note that the abundance cutoff with not vary in this analysis
#thus we classify images at this stage based on presence
#or absence of each cell type.
final_data = clin_summ %>%
  inner_join(mif_tumor$derived$bivariate_Count %>%
               rename(image_tag = image.tag)) %>%
  filter(r == 30,
         site == 'AAS',
         histotype == "high-grade serous") %>%
  mutate(anchor = case_when(anchor == "CD11b..CD15." ~ "CD11b+ CD15+ Positive Cells",
                            anchor == "CD11b..Opal.620..Positive" ~ "CD11b (Opal 620) Positive Cells",
                            anchor == "CD15..Opal.520..Positive" ~ "CD15 (Opal 520) Positive Cells",
                            anchor == "CD3..CD8." ~ "CD3+ CD8+ Positive Cells",
                            anchor == "CD3..FOXP3." ~ "CD3+ FOXP3+ Positive Cells",
                            anchor == "CD3..Opal.650..Positive" ~ "CD3 (Opal 650) Positive Cells",
                            anchor == "CD8..Opal.570..Positive" ~ "CD8 (Opal 570) Positive Cells",
                            anchor == "FOXP3..Opal.540..Positive" ~ "FOXP3 (Opal 540) Positive Cells"),
         counted = case_when(counted == "CD11b..CD15." ~ "CD11b+ CD15+ Positive Cells",
                            counted == "CD11b..Opal.620..Positive" ~ "CD11b (Opal 620) Positive Cells",
                            counted == "CD15..Opal.520..Positive" ~ "CD15 (Opal 520) Positive Cells",
                            counted == "CD3..CD8." ~ "CD3+ CD8+ Positive Cells",
                            counted == "CD3..FOXP3." ~ "CD3+ FOXP3+ Positive Cells",
                            counted == "CD3..Opal.650..Positive" ~ "CD3 (Opal 650) Positive Cells",
                            counted == "CD8..Opal.570..Positive" ~ "CD8 (Opal 570) Positive Cells",
                            counted == "FOXP3..Opal.540..Positive" ~ "FOXP3 (Opal 540) Positive Cells")) %>%
  rename("CD11b+ CD15+ Positive Cells" = cd11bplus_cd15plus_cells,
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
         "% FOXP3 (Opal 540) Positive Cells" = percent_foxp3_opal_540_positive_cells) %>%
  filter(anchor %in% c("CD3+ CD8+ Positive Cells", "CD3+ FOXP3+ Positive Cells") & 
           counted %in% c("CD3+ CD8+ Positive Cells", "CD3+ FOXP3+ Positive Cells")) %>%
  #The following condition corresponds to when there was at least one cell type that was not present
  #but due to mislabelling of IF software packages cells can both be classifed as`CD3+ FOXP3+ Positive Cells`
  #and `CD3+ CD8+ Positive Cells` thus we need to check that this did not make a cell count 0. 
  #(is.nan(`Degree of Clustering Permutation`) == TRUE &
            `CD3+ FOXP3+ Positive Cells` > `CD3+ CD8+ Positive Cells`)
  mutate(
    Treg_Abundance_PA = if_else(
      (`CD3+ FOXP3+ Positive Cells` !=
         0 & is.nan(`Degree of Clustering Permutation`) == FALSE)  |
        (
          is.nan(`Degree of Clustering Permutation`) == TRUE &
            `CD3+ FOXP3+ Positive Cells` > `CD3+ CD8+ Positive Cells`
        ),
      'P',
      'A'
    ),
    Cytotoxic_Abundance_PA = if_else((`CD3+ CD8+ Positive Cells` !=
                                        0 & is.nan(`Degree of Clustering Permutation`) == FALSE)  |
                                       (
                                         is.nan(`Degree of Clustering Permutation`) == TRUE &
                                           `CD3+ FOXP3+ Positive Cells` < `CD3+ CD8+ Positive Cells`
                                       ),
                                     'P',
                                     'A'
    ),
    Presence_anchor = if_else(anchor == "CD3+ FOXP3+ Positive Cells", Treg_Abundance_PA, Cytotoxic_Abundance_PA),
    Presence_counted = if_else(counted == "CD3+ FOXP3+ Positive Cells", Treg_Abundance_PA, Cytotoxic_Abundance_PA)
  ) %>%
  inner_join(fold_assignment)

pairs = data.frame(anchor = c("CD3+ CD8+ Positive Cells","CD3+ FOXP3+ Positive Cells"),
              counted = c("CD3+ FOXP3+ Positive Cells","CD3+ CD8+ Positive Cells"))
best = data.frame()

for(i in 1:nrow(pairs)){
#define the possible cut points for each marker
cuts = quantile(final_data %>% 
                  filter(anchor == pairs$anchor[i],
                         counted == pairs$counted[i]) %>%
                  select(`Degree of Clustering Permutation`) %>%
                  unlist(),
                na.rm = TRUE, probs = seq(0.2,0.8,0.05)) %>%
  data.frame(anchor = pairs$anchor[i], counted = pairs$counted[i], 
             left_out = i, 'cut' = .)

for(fold_lo in 1:10){
for(b in 1:nrow(cuts)){
#Construct model leaving one fold out and for a single marker 
tmp = final_data %>% 
  filter(fold != fold_lo,
         anchor == pairs$anchor[i],
         counted == pairs$counted[i]) %>%
  mutate(Spatial_HL = case_when(Presence_anchor != 'P' | Presence_counted != 'P' | 
                                  is.nan(`Degree of Clustering Permutation`) ~ 'N',
                                Presence_anchor == 'P' & Presence_counted == 'P' & 
                                  `Degree of Clustering Permutation` >= cuts$cut[b] ~ 'H',
                                Presence_anchor == 'P' & Presence_counted == 'P' & 
                                  `Degree of Clustering Permutation` < cuts$cut[b] ~ 'L'),
         final_group = paste0(Presence_anchor, Presence_counted, Spatial_HL))

#Full model (spatial + clinical)
model_lrt_1 = coxph(Surv(timelastfu_new, as.numeric(vitalstatus_new)) ~ refage + stage + final_group + cluster(suid), 
                    data = tmp)
                    #Reduced model (clinical only)
model_lrt_2 = coxph(Surv(timelastfu_new, as.numeric(vitalstatus_new)) ~ refage + stage + cluster(suid), 
                    data = tmp)
                    
#compute the log ratio test p-value
cuts$lrt_log_lik_full[b] = summary(model_lrt_1)$loglik[2]
cuts$lrt_log_lik_clinical[b] = summary(model_lrt_2)$loglik[2]
cuts$lrt[b] = -2*log(cuts$lrt_log_lik_full[b] - cuts$lrt_log_lik_clinical[b])
cuts$big_group[b] = all(table(tmp$final_group)[4:5] >= 10) #check if there are more than 10 images between the two groups
cuts$n_AAN[b] = table(tmp$final_group)['AAN']
cuts$n_APN[b] = table(tmp$final_group)['APN']
cuts$n_PAN[b] = table(tmp$final_group)['PAN']
cuts$n_PPL[b] = table(tmp$final_group)['PPL']
cuts$n_PPH[b] = table(tmp$final_group)['PPH']
cuts$total[b] = sum(table(tmp$final_group))
}
  
  best = plyr::rbind.fill(best,
                          cuts %>% 
                            filter(!is.na(lrt)) %>%
                            arrange(lrt) %>%
                            head(1) %>%
                            mutate(left_out = fold_lo,
                                   anchor = pairs$anchor[i],
                                   counted = pairs$counted[i]))
}
}

#Getting the model summary for each marker
final_cut = best %>% group_by(anchor, counted) %>%
  filter(big_group == TRUE) %>%
  select(cut) %>%
  summarize_all(~median(.))


#Use the cutpoint that lead to the biggest change in log liklihood (after adding the spatial information) 
model = lapply(1:nrow(final_cut), function(i){
  tmp = final_data %>% 
    filter(anchor == final_cut$anchor[i],
           counted == final_cut$counted[i]) %>%
    mutate(Spatial_HL = case_when(Presence_anchor != 'P' | Presence_counted != 'P' | 
                                    is.nan(`Degree of Clustering Permutation`) ~ 'N',
                                  Presence_anchor == 'P' & Presence_counted == 'P' & 
                                    `Degree of Clustering Permutation` >= final_cut$cut[i] ~ 'H',
                                  Presence_anchor == 'P' & Presence_counted == 'P' & 
                                    `Degree of Clustering Permutation` < final_cut$cut[i] ~ 'L'),
           final_group = paste0(Presence_anchor, Presence_counted, Spatial_HL)) %>%
    data.frame(check.names = FALSE)
  
  model = coxph(Surv(timelastfu_new, as.numeric(vitalstatus_new)) ~ refage + stage + final_group + cluster(suid), 
                data = tmp)
  model_2 = coxph(Surv(timelastfu_new, as.numeric(vitalstatus_new)) ~ refage + stage + cluster(suid), 
                  data = tmp)
  lrt_log_lik_full = summary(model)$loglik[2]
  lrt_log_lik_clinical = summary(model_2)$loglik[2]
  lrt = 2*(lrt_log_lik_full - lrt_log_lik_clinical)
  print(c(i,lrt, 1-pchisq(lrt,4)))
  
  plot = ggadjustedcurves(model, data = tmp, variable = 'final_group', size = 0.75) +   
    ggtitle(paste0('Anchor = ', pairs$anchor[i], ', Counted = ', pairs$counted[i])) + 
    labs(col="Group") + xlab('Time (Days)') + 
    theme(legend.title = element_text(size=20),  legend.text = element_text(size = 18))
  print(plot)
  return(list(model = model, plot = plot, 
              n = table(tmp$final_group)))
}
)



#Getting the model summary for each marker
final_models_output_tma = lapply(1:nrow(final_cut), function(a){
  model_fits = c(anchor = final_cut$anchor[a],
                 counted = final_cut$counted[a],
                 round(coefficients(summary(model[[a]]$model))[4:7, 'exp(coef)'],4),
                 round(coefficients(summary(model[[a]]$model))[4:7, 'Pr(>|z|)'], 4),
                 model[[a]]$n)
  return(model_fits)
}) %>% plyr::ldply()
colnames(final_models_output_tma)[7:10] = paste0(colnames(final_models_output_tma)[7:10], '_p.val')
