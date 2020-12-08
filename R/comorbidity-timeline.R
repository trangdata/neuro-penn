# Comorbidity timesline
library(tidyverse)
library(ggplot2)
library(scales)

comorb_time <- function(data, mapped_codes, comorbidity, day_start){

  # merge by mapped comorbidity codes
  data <- merge(data, mapped_codes, by = "concept_code")

  # sum unique counts
  data <- data %>%
    select(patient_num, Abbreviation, days_since_admission) %>%
    group_by(Abbreviation, days_since_admission) %>%
    distinct() %>% # patients may have more than one code for a comorbidity
    mutate(count = n()) %>%
    distinct(Abbreviation, days_since_admission, count) %>%
    ungroup()

  #count number of unique occurrences by patient
  t_plot <- ggplot(data %>% filter(Abbreviation == comorbidity,
                         days_since_admission >= day_start),
         aes(x = days_since_admission, y = count)) +
    geom_line(color = "tomato", size = 1) +
    geom_point(color = "tomato", size = 0.5) +
    scale_x_continuous(name ="Days Since Initial COVID-19 Admission",
                     breaks = pretty_breaks(n = 10)) +
    #expand_limits(x = 0, y = 0) +
    ggtitle(paste(comorbidity))

  t_plot <- t_plot + theme(plot.title = element_text(hjust = 0.5))

  return(t_plot)

}
