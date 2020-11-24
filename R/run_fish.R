# library(DT)
library(glue)
# ori_variable = 'ever_severe'
run_fish <- function(ori_variable = c('ever_severe', 'death', 'readmitted')){
  ori_variable <- match.arg(ori_variable)

  x_labs <- switch(
    ori_variable,
    ever_severe = 'Number of ever severe patients',
    death = 'Number of patients who died',
    readmitted = 'Number of patients who were readmitted'
  )

  severe_df <- nat_hist_df %>%
    rename('new_var' = !!sym(ori_variable))

  n_severe <- severe_df %>% pull(new_var) %>% sum()
  n_non_severe <- severe_df %>% pull(new_var) %>% `!` %>% sum()

  contingency_df <- neuro_patients %>%
    left_join(severe_df, by = 'patient_num') %>%
    group_by(concept_code) %>%
    summarise(num_patients_ever_severe_icd1 = sum(new_var),
              num_patients_never_severe_icd1 = sum(!new_var),
              num_patients_ever_severe_icd0 = n_severe - num_patients_ever_severe_icd1,
              num_patients_never_severe_icd0 = n_non_severe - num_patients_never_severe_icd1,
              .groups = 'keep')

  nested_obs_exp <- nest(contingency_df)

  fish_obs_exp <- nested_obs_exp %>%
    mutate(fish = map(data, my_fish_warning)) %>%
    dplyr::select(-data) %>%
    unnest(cols = c(fish)) %>%
    mutate(upper = ifelse(is.na(upper), Inf, upper)) %>%
    mutate(
      lestimate = log2(estimate),
      llower = log2(lower),
      lupper = log2(upper),
      lci = paste0('(', round(llower, 2), ', ',
                   round(lupper, 2), ')'),
    )
  fish_obs_exp$p_holm <-  p.adjust(fish_obs_exp$p_value, 'holm')
  fish_obs_exp$p_fdr <- p.adjust(fish_obs_exp$p_value, 'BH')

  fish_tab <- fish_obs_exp  %>%
    left_join(contingency_df, by = 'concept_code') %>%
    mutate(Observed = num_patients_ever_severe_icd1,
           Expected = num_patients_never_severe_icd1/n_non_severe*n_severe,
           over_sev = Observed - Expected,
           drr = round((estimate - 1)*100, 0),
           lower_drr = (lower - 1)*100,
           upper_drr = (upper - 1)*100,
           ci_drr = paste0('(', round(lower_drr, 0), ', ',
                           round(upper_drr, 0), ')')) %>%
    arrange(p_fdr)

  # fish_tab %>%
  #   rename('ICD' = 'concept_code',
  #          'Observed - Expected' = 'over_sev',
  #          'Enrichment' = 'estimate',
  #          'Log2(enrichment)' = 'lestimate',
  #          '95% Confidence interval' = 'lci',
  #          'p_raw' = 'p_value') %>%
  #   datatable(rownames = FALSE, filter = 'top') %>%
  #   formatRound(c('Observed', 'Expected', 'Observed - Expected',
  #                 'Enrichment', 'Log2(enrichment)'), 1) %>%
  #   formatSignif(c('p_raw', 'p_holm', 'p_fdr'), 3) %>%
  #   {.}

  fish_tab %>%
    write_csv(glue('results/icd10_{ori_variable}_raw.csv'))

  fish_tab %>%
    mutate(lestimate = round(lestimate, 2),
           p_fdr = format(p_fdr, digits = 2)) %>%
    select(concept_code, warning, lestimate, drr, ci_drr, p_fdr) %>%
    write_csv(glue('results/icd10_{ori_variable}_show.csv'))

  filtered_obs_exp <- fish_tab %>%
    mutate(
      distance_to_null = case_when(
        lower > 1 ~ lower - 1,
        TRUE ~ upper - 2
      ),
      presentation = case_when(
        lower > 1 & p_fdr < 0.05 ~ '#d46780',
        upper < 1 & p_fdr < 0.05 ~ '#798234',
        TRUE ~ 'grey20'
      )) %>%
    left_join(neuro_icds_10, by =c('concept_code' = 'icd')) %>%
    mutate(full_icd = paste0(`ICD-10 Description`, ' (', concept_code, ')') %>%
             as.factor() %>% fct_reorder(Expected)) %>%
    arrange(Expected) %>%
    {.}

  plot_obs_exp <- filtered_obs_exp %>%
    select(concept_code, full_icd, lestimate, llower, lupper,
           presentation, over_sev, Observed, Expected) %>%
    pivot_longer(- c(concept_code, full_icd, presentation, over_sev), names_to = 'type') %>%
    mutate(subtype = ifelse(type == 'Expected' | type == 'Observed',
                            'Sqrt(number of patients)',
                            'Log2 enrichment, 95% CI')) %>%
    pivot_wider(names_from = type) %>%
    mutate(presentation = as.factor(presentation),
           concept_code = fct_reorder(concept_code, Observed))

  plot_obs_exp_right <- plot_obs_exp %>%
    filter(subtype == 'Sqrt(number of patients)')

  plot_obs_exp_left <- plot_obs_exp %>%
    filter(subtype != 'Sqrt(number of patients)')

  title <- cowplot::ggdraw() +
    cowplot::draw_label(ori_variable, x = 0, hjust = 0) +
    theme(plot.margin = margin(0, 0, 0, 7))

  enrichment_plot_after <- cowplot::plot_grid(
    title, plot_enrich(
      plot_obs_exp_left,
      plot_obs_exp_right,
      filtered_obs_exp,
      xlab = x_labs,
      nudge = 2),
    ncol = 1,
    rel_heights = c(0.05, 1)
  )
  print(enrichment_plot_after)

  ggsave(glue('figs/icd_{ori_variable}_after.png'),
         enrichment_plot_after, width = 9.5, height = 3.5)
}
