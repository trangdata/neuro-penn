library(icd)
library(tidyverse)
library(tableone)


first_3 <- function(x) {
  # retain first 3 characters of the ICD code
  substr(x, 1, 3) %>% unique()
}

# The below was copied from the Comorbidity_Mapping.Rmd file in order to support sourcing it for other analyses.
# data = obs_raw
map_charlson_codes <- function(data) {
  data <- data %>%
    filter(concept_type %in% c("DIAG-ICD10", "DIAG-ICD9"),
           # filter for diagnoses prior to admission
           days_since_admission < 0)

  # Create separate data frames for ICD9 and 10 Codes
  # icd package does not support simultaneous processing of both ICD code types
  # we will recombine after the initial processing
  icd10 <- data %>%
    filter(concept_type == "DIAG-ICD10") %>%
    distinct()

  icd9 <- data %>%
    filter(concept_type == "DIAG-ICD9") %>%
    distinct()

  ## Because the 4CE has truncated ICD codes, we will also truncate the icd package index maps
  # Function to select first 3 characters of the ICD Code in all lists of the index map
  icd10_map_charlson3 <- lapply(icd10_map_charlson, first_3)
  icd9_map_charlson3 <- lapply(icd9_map_charlson, first_3)

  # perform the mapping
  icd10_map <-
    icd10_comorbid(
      icd10,
      map = icd10_map_charlson3,
      icd_name = "concept_code",
      return_df = TRUE,
      visit_name = "patient_num",
      return_binary = TRUE
    )

  icd9_map <-
    icd9_comorbid(
      icd9,
      map = icd9_map_charlson3,
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

    # replace abbreviations with full comorbidity name
  comorb_list_names <-
    map_df(
      names_charlson, ~ as.data.frame(.x),
      .id  =  "code"
    ) %>%
    full_join(
      map_df(names_charlson_abbrev, ~ as.data.frame(.x),
             .id="code"),
      by = "code"
    ) %>%
    select(- code) %>%
    `colnames<-`(c('Name', 'Comorbidity'))

  table1 <- icd_map %>%
    select(- patient_num) %>%
    pivot_longer(everything(), names_to = 'Comorbidity') %>%
    group_by(Comorbidity) %>%
    summarise(Patients = sum(value), .groups = 'drop') %>%
    left_join(comorb_list_names, by = "Comorbidity") %>%
    select(Comorbidity = Name, Patients) %>%
    mutate(percent = round(Patients / nrow(icd_map) * 100, 2)) %>%
    arrange(
      desc(Patients)
    )

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

  ## Identify the specific codes that mapped

  # unlist the charlson mapping lists
  icd10_map <-
    map_df(icd10_map_charlson3, ~ as.data.frame(.x), .id = "name") %>%
    `colnames<-`(c("Comorbidity", "concept_code")) %>%
    filter(!concept_code %in% comorb_list_names$Comorbidity) %>%
    distinct() %>%
    # merge the mapping dataframe to the patient level ICD codes
    # this will return all comorbidities that mapped to our patient data
    inner_join(icd10, by = "concept_code")

  icd9_map <-
    map_df(icd9_map_charlson3, ~ as.data.frame(.x), .id = "name") %>%
    `colnames<-`(c("Comorbidity", "concept_code")) %>%
    filter(!concept_code %in% comorb_list_names$Comorbidity) %>%
    distinct() %>%
    inner_join(icd9, by = "concept_code")

  # explain_codes will add additional information regarding the code name
  # add if statements in order to handle sites that only have ICD 9 or 10 codes but not both

  if (nrow(icd10_map) > 0) {
    icd10_mapped_table <- explain_table(icd10_map$concept_code) %>%
      distinct() %>%
      left_join(icd10_map, ., by = c("concept_code" = "code")) %>%
      select(patient_num, concept_code, Comorbidity, long_desc) %>%
      distinct()
  }

  if (nrow(icd9_map) > 0) {
    icd9_mapped_table <- explain_table(icd9_map$concept_code) %>%
      distinct() %>%
      left_join(icd9_map, ., by = c("concept_code" = "code")) %>%
      select(patient_num, concept_code, Comorbidity, long_desc) %>%
      distinct()
  }

  # Bind both ICD 9 and 10 code tables together
  mapped_codes_table <-
    rbind(try(icd10_mapped_table, icd9_mapped_table), silent = TRUE)

  # calculate how many patients had each unique Comorbidity/concept_code
  mapped_codes_table <- mapped_codes_table %>%
    add_count(concept_code, Comorbidity, long_desc, name = 'n_patients') %>%
    group_by(long_desc, n_patients) %>%
    arrange(desc(n_patients)) %>%
    select(- patient_num) %>%
    distinct()

  map_results <-
    list(
      comorb_list_names = comorb_list_names,
      icd_map = icd_map,
      index_scores = index_scores,
      table1 = table1,
      mapped_codes_table = mapped_codes_table
    )

  map_results

}
