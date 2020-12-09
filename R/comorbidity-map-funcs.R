library(icd)
library(tidyverse)
library(tableone)

# The below was copied from the Comorbidity_Mapping.Rmd file in order to support sourcing it for other analyses.
# data = obs_raw
# comorb_names <- get_comorb_names()

map_charlson_codes <- function(data, comorb_names, num_days_since_admission, truncate = TRUE) {
  data <- data %>%
    filter(concept_type %in% c("DIAG-ICD10", "DIAG-ICD9"),
           # filter for diagnoses prior to admission
           days_since_admission <= num_days_since_admission)

  # Create separate data frames for ICD9 and 10 Codes
  # icd package does not support simultaneous processing of both ICD code types
  # we will recombine after the initial processing
  icd10 <- data %>%
    select(-c(days_since_admission, value)) %>%
    filter(concept_type == "DIAG-ICD10") %>%
    distinct()

  icd9 <- data %>%
    select(-c(days_since_admission, value)) %>%
    filter(concept_type == "DIAG-ICD9") %>%
    distinct()

  ## Because the 4CE has truncated ICD codes, we will also truncate the icd package index maps
  # Function to select first 3 characters of the ICD Code in all lists of the index map
  if (truncate == TRUE) {
    icd10_map_charlson <- lapply(icd10_map_charlson, first_3)
    icd9_map_charlson <- lapply(icd9_map_charlson, first_3)
  }

  if (truncate == FALSE) {
    # convert icd code to short format (without decimals)
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
      map = icd10_map_charlson,
      icd_name = "concept_code",
      return_df = TRUE,
      visit_name = "patient_num",
      return_binary = TRUE
    )

  icd9_map <-
    icd9_comorbid(
      icd9,
      map = icd9_map_charlson,
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
    ) %>%
    arrange(as.numeric(patient_num))

  ## Calculate Index Scores
  charlson_score <- charlson_from_comorbid(
    icd_map,
    visit_name = "patient_num",
    scoring_system = "charlson",
    hierarchy = TRUE
  ) %>%
    data.frame(charlson_score = .) %>%
    tibble::rownames_to_column("patient_num")

  quan_score <- charlson_from_comorbid(
    icd_map,
    visit_name = "patient_num",
    scoring_system = "quan",
    hierarchy = TRUE
  ) %>%
    data.frame(quan_score = .) %>%
    tibble::rownames_to_column("patient_num")

  index_scores <- icd_map %>%
    full_join(charlson_score, by = "patient_num") %>%
    full_join(quan_score, by = "patient_num") %>%
    # calculate difference between Charlson and Quan scores
    mutate(score_diff = charlson_score - quan_score,
           abs_score_diff = abs(score_diff)) %>%
    arrange(desc(charlson_score)) %>%
    select(
      patient_num,
      charlson_score,
      quan_score,
      score_diff,
      abs_score_diff,
      everything()
    )

  # Identify the specific codes that mapped
  # unlist the charlson mapping lists
  icd10$concept_code <- as.character(icd10$concept_code)
  icd9$concept_code <- as.character(icd9$concept_code)

  icd10_map <-
    map_df(icd10_map_charlson, ~ as.data.frame(.x), .id = "name") %>%
    `colnames<-`(c("Abbreviation", "concept_code")) %>%
    filter(!concept_code %in% comorb_names$Abbreviation) %>%
    mutate(concept_code = as.character(concept_code)) %>%
    distinct() %>%
    # merge the mapping dataframe to the patient level ICD codes
    # this will return all comorbidities that mapped to our patient data
    inner_join(icd10, by = "concept_code")

  icd9_map <-
    map_df(icd9_map_charlson, ~ as.data.frame(.x), .id = "name") %>%
    `colnames<-`(c("Abbreviation", "concept_code")) %>%
    filter(!concept_code %in% comorb_names$Abbreviation) %>%
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
    {if (nrow(icd9_map) > 0) bind_rows(icd9_mapped_table) else .} %>%
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

get_table1 <- function(
  df, comorbidities = comorb_names$Abbreviation
){
  df %>%
    select(all_of(comorbidities)) %>%
    colSums() %>%
    data.frame(n_patients = .) %>%
    rownames_to_column("Abbreviation") %>%
    right_join(comorb_names, ., by = "Abbreviation")
}

get_comorb_names <- function(){
  data.frame(
    Comorbidity = do.call(rbind, names_charlson),
    Abbreviation = do.call(rbind, names_charlson_abbrev))
}

first_3 <- function(x) {
  # retain first 3 characters of the ICD code
  substr(x, 1, 3) %>% unique()
}
