concat <- function(x, y){
  paste0(x, ' (', round(y, 3)*100, '%)')
}

concat_median <- function(med, mi, ma){
  paste0(med, ' [', mi, ', ', ma, ']')
}

concat_mean <- function(mea, s, acc = 0){
  paste0(round(mea, acc), ' (', round(s, acc), ')')
}

count_stats <- function(df, count_var, neg_var, ...) {
  # summary statistics for survival/severity status
  # count values are obfuscated

  count_var <- sym(count_var)
  neg_var <- sym(neg_var)
  Count_var <- sym(stringr::str_to_title(count_var))

  df %>%
    select(neuro_post, count_var) %>%
    group_by(neuro_post) %>%
    summarise(!!neg_var := sum(!!count_var == 0),
              Total = n(),
              .groups = 'drop') %>%
    blur_it(c(neg_var, 'Total'), ...) %>%
    mutate(!!count_var := Total - !!neg_var) %>%
    transmute(
      neuro_post,
      !!neg_var := concat(!!neg_var, !!neg_var/Total),
      !!Count_var := concat(!!count_var, !!count_var/Total)) %>%
    pivot_longer(- neuro_post) %>%
    pivot_wider(names_from = neuro_post, values_from = value)
}

continuous_stats <- function(df, cont_var, name, ...) {
  # summary statistics for length of stay
  # count values are obfuscated
  med <- sym(paste('Median', name, '[Min, Max]'))
  mea <- sym(paste('Mean', name, '(SD)'))
  cont_var <- sym(cont_var)

  df %>%
    select(neuro_post, cont_var) %>%
    group_by(neuro_post) %>%
    summarise(median_var = median(!!cont_var, na.rm = TRUE),
              min_var = min(!!cont_var, na.rm = TRUE),
              max_var = max(!!cont_var, na.rm = TRUE),
              mean_var = mean(!!cont_var, na.rm = TRUE),
              sd_var = sd(!!cont_var, na.rm = TRUE),
              .groups = 'drop') %>%
    transmute(
      neuro_post,
      !!med := concat_median(median_var, min_var, max_var),
      !!mea := concat_mean(mean_var, sd_var)) %>%
    pivot_longer(- neuro_post) %>%
    pivot_wider(names_from = neuro_post, values_from = value)
}

demo_stats <- function(df, var, ...){
  svar <- sym(var)
  df %>%
    group_by(!!svar) %>%
    count(neuro_post, name = 'n_var') %>%
    as.data.frame() %>%
    blur_it('n_var', ...) %>%
    group_by(neuro_post) %>%
    mutate(both_neuro = sum(n_var)) %>%
    ungroup() %>%
    mutate(prop = n_var/both_neuro,
           pres = concat(n_var, n_var/both_neuro)) %>%
    pivot_wider(- c(n_var, both_neuro), names_from = neuro_post,
                values_from = c(n_var, prop, pres)) %>%
    mutate(variable = paste(var, !!svar, sep = '.')) %>%
    replace_na(list(pres_no_neuro_cond = '0 (0%)',
                    pres_neuro_cond = '0 (0%)')) %>%
    select(- !!svar)
}


blur_it <- function(df, vars, blur_abs, mask_thres){
  # Obfuscate count values.
  # If blurring range is +/-3, or blur_abs = 3,
  # the count receive a small addition of a random number from -3 to 3.
  # If a count is less than mask_thres, set that count to 0.

  for (var in vars){
    var <- sym(var)
    blur_vec <- sample(seq(- blur_abs, blur_abs), nrow(df), replace = TRUE)
    df <- df %>%
      mutate(!!var := !!var + blur_vec,
             !!var := ifelse(abs(!!var) < mask_thres, 0, !!var))
  }
  df
}
