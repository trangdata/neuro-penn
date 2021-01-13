library(icd)
library(tidyverse)

# The below was copied from the Comorbidity_Mapping.Rmd file in order to support sourcing it for other analyses.
# df = obs_raw # (4ce LocalPatientObservations.csv)
# comorb_names <- get_charlson_names() or get_elix_names()
# t1 <- earliest time point to consider comorbidities
# t2 <- latest time point to consider comorbidities
# t1 = -365; t2 = -1;
# example above will map all codes up to a year prior but before admission (admission = day 0)
# num_days_prior_admission = -365 indicates that we consider all codes up to a year prior to the first COVID admission as comorbidities
# day_of
# map_type = 'charlson', 'elixhauser' - where charlson will be scored with quan-deyo
# truncate = TRUE # indicates we are using ICD codes truncated to the first 3 characters; set FALSE if you have full ICD codes
concat <- function(x, y){
  paste0(x, ' (', round(y, 3)*100, '%)')
}

map_char_elix_codes <- function(df, comorb_names, t1, t2, map_type, truncate = TRUE) {

  df <- df %>%
    filter(concept_type %in% c("DIAG-ICD10", "DIAG-ICD9"),
           days_since_admission >= t1 & days_since_admission <= t2)

  # Create separate df frames for ICD9 and 10 Codes
  # icd package does not support simultaneous processing of both ICD code types
  # we will recombine after the initial processing
  icd10 <- df %>%
    select(-c(days_since_admission, value)) %>%
    filter(concept_type == "DIAG-ICD10") %>%
    distinct()

  icd9 <- df %>%
    select(-c(days_since_admission, value)) %>%
    filter(concept_type == "DIAG-ICD9") %>%
    distinct()

  # select map_type for analysis
  #if (map_type == "charlson") {
  #  icd10_comorb_map = icd10_map_charlson
  #  icd9_comorb_map = icd9_map_charlson
  #}

  # quan prefixes revised charlson & elixhauser mapping
  if (map_type == "charlson") {
    icd10_comorb_map = icd10_map_quan_deyo
    icd9_comorb_map = icd9_map_quan_deyo
  }
  if (map_type == "elixhauser") {
    icd10_comorb_map = icd10_map_quan_elix
    icd9_comorb_map = icd9_map_quan_elix
  }

  ## Because the 4CE has truncated ICD codes, we will also truncate the icd package index maps
  # Function to select first 3 characters of the ICD Code in all lists of the index map
  if (truncate == TRUE) {
    icd10_comorb_map <- lapply(icd10_comorb_map, first_3)
    icd9_comorb_map <- lapply(icd9_comorb_map, first_3)
  }

  if (truncate == FALSE) {
    # convert icd code to short format (without decimals) to facilitate mapping
    # where diagnosis code is you non-truncated icd column
    icd10 <- icd10 %>%
      mutate(diagnosis_code = decimal_to_short(diagnosis_code)) %>%
      select(-concept_code) %>%
      rename(concept_code = diagnosis_code)

    icd9 <- icd9 %>%
      mutate(diagnosis_code = decimal_to_short(diagnosis_code)) %>%
      select(-concept_code) %>%
      rename(concept_code = diagnosis_code)

  }

  # perform the mapping
  icd10_map <-
    icd10_comorbid(
      icd10,
      map = icd10_comorb_map,
      return_df = TRUE,
      visit_name = "patient_num",
      return_binary = TRUE,
    )

  icd9_map <-
    icd9_comorbid(
      icd9,
      map = icd9_comorb_map,
      icd_name = "concept_code",
      return_df = TRUE,
      visit_name = "patient_num",
      return_binary = TRUE
    )

  # If multiple rows due to a patient having both ICD 9 and 10 codes, we will take the max of the column
  # This will allow us to capture the 1s indicating that the comorbidity is present
  # the try wrapper is important in cases there is not an instance of a specific comorbidity in the data - try silences errors
  icd_map <- rbind(icd9_map, icd10_map) %>% # results of 9 and 10 mapping
    group_by(patient_num) %>%
    summarise(
      across(everything(), ~ try(max(.x), silent = TRUE)),
      .groups = 'drop'
    )

  ## Calculate Index Scores
  if (map_type == 'charlson') {
    charlson_score <- charlson_from_comorbid(
      icd_map,
      visit_name = "patient_num",
      scoring_system = "charlson",
      hierarchy = TRUE
    ) %>%
      data.frame(charlson_score = .) %>%
      tibble::rownames_to_column("patient_num")
  }

  # need to check If I've done this correctly - it seems to be that their are different versions
  # of the elixhuaser mapping and that in one HTN is combined. I think this is the version needed
  # to run the van_walraven_from_comorb function
  if (map_type == "elixhauser") {
    # combine hypertension into one category
    icd_map <- icd_map %>%
      mutate(HTN = pmax(HTN, HTNcx, na.rm = TRUE)) %>%
      select(-HTNcx)

    van_walraven_score <- van_walraven_from_comorbid(
      icd_map,
      visit_name = 'patient_num',
      hierarchy = TRUE
    ) %>%
      data.frame(van_walraven_score = .) %>%
      tibble::rownames_to_column("patient_num")
  }

  if(map_type == 'charlson') {
    index_scores <- icd_map %>%
      full_join(charlson_score, by = "patient_num") %>%
      arrange(desc(charlson_score)) %>%
      select(
        patient_num,
        charlson_score,
        everything()
      )
  } else {
    index_scores <- icd_map %>%
      full_join(van_walraven_score, by = "patient_num") %>%
      arrange(desc(van_walraven_score)) %>%
      select(
        patient_num,
        van_walraven_score,
        everything()
      )
  }


  # Identify the specific codes that mapped
  # unlist the charlson mapping lists
  icd10$concept_code <- as.character(icd10$concept_code)
  icd9$concept_code <- as.character(icd9$concept_code)

  icd10_map <-
    map_df(icd10_comorb_map, ~ as.data.frame(.x), .id = "name") %>%
    `colnames<-`(c("Abbreviation", "concept_code")) %>%
    mutate(concept_code = as.character(concept_code)) %>%
    distinct() %>%
    # merge the mapping dataframe to the patient level ICD codes
    # this will return all comorbidities that mapped to our patient data
    inner_join(icd10, by = "concept_code")

  icd9_map <-
    map_df(icd9_comorb_map, ~ as.data.frame(.x), .id = "name") %>%
    `colnames<-`(c("Abbreviation", "concept_code")) %>%
    mutate(concept_code = as.character(concept_code)) %>%
    distinct() %>%
    inner_join(icd9, by = "concept_code")

  # explain_codes will add additional information regarding the code name
  # add if statements in order to handle sites that only have ICD 9 or 10 codes but not both
  if (nrow(icd10_map) > 0) {
    icd10_mapped_table <- unique(icd10_map$concept_code) %>%
      explain_table() %>%
      left_join(icd10_map, ., by = c("concept_code" = "code")) %>%
      select(patient_num, concept_code, Abbreviation, long_desc) %>%
      distinct()
  }

  if (nrow(icd9_map) > 0) {
    icd9_mapped_table <- unique(icd9_map$concept_code) %>%
      explain_table() %>%
      left_join(icd9_map, ., by = c("concept_code" = "code")) %>%
      select(patient_num, concept_code, Abbreviation, long_desc) %>%
      distinct()
  }

  # Bind both ICD 9 and 10 code tables together
  mapped_codes_table <-
    icd10_mapped_table %>%
    try(bind_rows(icd9_mapped_table), silent = TRUE) %>%
    # calculate how many patients had each unique comorbidity/concept_code
    count(concept_code, Abbreviation, long_desc, name = 'n_patients') %>%
    arrange(desc(n_patients))

  map_results <-
    list(
      index_scores = index_scores,
      mapped_codes_table = mapped_codes_table
    )

  map_results

}

# where df takes in the matrix for the initial mapping

get_table1 <- function(
  df,
  comorbidities = comorb_name_df,
  comorbidities_map = comorbidities$Abbreviation
)
  {
  df %>%
    select(all_of(comorbidities_map)) %>%
    colSums() %>%
    data.frame(n_patients = .) %>%
    rownames_to_column("Abbreviation") %>%
    right_join(comorbidities, ., by = "Abbreviation")
}

process_tables <- function(index_scores, ...) {
  get_table1(index_scores %>% filter(patient_num %in% neuro_pt_post), ...) %>%
    rename('n_neuro_pats' = n_patients) %>%
    left_join(get_table1(index_scores, ...),
              by = c("Comorbidity", "Abbreviation"))
}

get_charlson_names <- function(){
  data.frame(
    Comorbidity = do.call(rbind, names_charlson),
    Abbreviation = do.call(rbind, names_charlson_abbrev))
}

get_quan_elix_names <- function(){
  data.frame(
    Comorbidity = do.call(rbind, names_quan_elix),
    Abbreviation = do.call(rbind, names_quan_elix_abbrev))
}

first_3 <- function(x) {
  # retain first 3 characters of the ICD code
  substr(x, 1, 3) %>% unique()
}

concat_median <- function(med, mi, ma){
  paste0(med, ' [', mi, ', ', ma, ']')
}

concat_mean <- function(mea, s, acc = 0){
  paste0(round(mea, acc), ' (', round(s, acc), ')')
}

severity_stats <- function(df, neuro_cond, ...) {
  # summary statistics for severity status
  # count values are obfuscated
  df %>%
    select(neuro_post, time_severe = `Time to severity onset (days)`) %>%
    group_by(neuro_post) %>%
    summarise(median_time = median(time_severe, na.rm = TRUE),
              min_time = min(time_severe, na.rm = TRUE),
              max_time = max(time_severe, na.rm = TRUE),
              mean_time = mean(time_severe, na.rm = TRUE),
              sd_time = sd(time_severe, na.rm = TRUE),
              non_severe = sum(is.na(time_severe)),
              Total = n(),
              .groups = 'drop') %>%
    blur_it(c('non_severe', 'Total'), ...) %>%
    mutate(severe = Total - non_severe) %>%
    transmute(
      neuro_post,
      Nonsevere = concat(non_severe, non_severe/Total),
      Severe = concat(severe, severe/Total),
      `Median time to severity onset [Min, Max] (days)` = concat_median(median_time, min_time, max_time),
      `Mean time to severity onset (SD) (days)` = concat_mean(mean_time, sd_time)) %>%
    pivot_longer(-neuro_post) %>%
    pivot_wider(names_from = neuro_post, values_from = value) %>%
    column_to_rownames('name')
}

survival_stats <- function(df, neuro_cond, ...) {
  # summary statistics for survival status
  # count values are obfuscated
  df %>%
    select(neuro_post, time_death = `Time to death (days)`) %>%
    group_by(neuro_post) %>%
    summarise(median_time = median(time_death, na.rm = TRUE),
              min_time = min(time_death, na.rm = TRUE),
              max_time = max(time_death, na.rm = TRUE),
              mean_time = mean(time_death, na.rm = TRUE),
              sd_time = sd(time_death, na.rm = TRUE),
              alive = sum(is.na(time_death)),
              Total = n(),
              .groups = 'drop') %>%
    blur_it(c('alive', 'Total'), ...) %>%
    mutate(deceased = Total - alive) %>%
    transmute(
      neuro_post,
      Alive = concat(alive, alive/Total),
      Deceased = concat(deceased, deceased/Total),
      `Median time to death [Min, Max] (days)` = concat_median(median_time, min_time, max_time),
      `Mean time to death (SD) (days)` = concat_mean(mean_time, sd_time)) %>%
    pivot_longer(-neuro_post) %>%
    pivot_wider(names_from = neuro_post, values_from = value) %>%
    column_to_rownames('name')
}

blur_it <- function(df, vars, blur_abs, mask_thres){
  # Obfuscate count values.
  # If blurring range is +/-3, or blur_abs = 3,
  # the count receive a small addition of a random number from -3 to 3.
  # If a count is smaller than mask_thres, set that count to 0.

  for (var in vars){
    var <- sym(var)
    blur_vec <- sample(seq(- blur_abs, blur_abs), nrow(df), replace = TRUE)
    df <- df %>%
      mutate(!!var := ifelse(!!var < mask_thres, 0, !!var  + blur_vec))
  }
  df
}
