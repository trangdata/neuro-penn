# scores_unique %>%
#   ggplot(aes(x = charlson_score, fill = neuro_post)) +
#   geom_histogram(binwidth = 1) +
#   NULL
#
# index_scores %>%
#   pivot_longer(MI:HIV) %>%
#   filter(value != 0) %>%
#   left_join(comorb_names, c('name' = 'Abbreviation')) %>%
#   mutate(fullname = glue::glue('{Comorbidity} ({name})')) %>%
#   ggplot(aes(x = charlson_score, y = fct_reorder(fullname, charlson_score))) +
#   geom_boxplot() +
#   labs(y = NULL)
#
# comorb_long %>%
#   ggplot(aes(
#     y = reorder(name, - charlson_score, na.rm = TRUE),
#     x = charlson_score, fill = neuro_post)) +
#   geom_boxplot() +
#   labs(x = 'Charlson score', y = NULL) +
#   scale_fill_manual(values = c("gray", "slateblue"),
#                     guide = guide_legend(reverse = TRUE)) +
#   NULL

### By grouped neuro diagnosis and severity

# How would we interpret this figure?
#
#
# scores_neuro %>%
#   ggplot(aes(
#     y = reorder(`Neurological Disease Category`, - charlson_score, na.rm = TRUE),
#     x = charlson_score, fill = severe)) +
#   geom_boxplot() +
#   labs(x = 'Charlson score', y = NULL) +
#   scale_fill_manual(name = NULL,
#                     values = c("gray", "darkorange"),
#                     guide = guide_legend(reverse = TRUE)) +
#   NULL


### Median and IQR

# aggregate(charlson_score ~ neuro_post, FUN = "median", scores_unique)
# aggregate(charlson_score ~ neuro_post, FUN = "summary", scores_unique)

