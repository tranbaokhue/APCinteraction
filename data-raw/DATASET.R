# Load required packages
library(tidyverse)

# Read and process the dataset
FirstNullAPCSSA <- read_csv("data-raw/NullDistributionAPCSSA.csv") %>%
  mutate(
    SD_CRA = sqrt(V_CRA),
    SD_RCA = sqrt(V_RCA)
  ) %>%
  select(I, J, K, E_CRA, SD_CRA, E_RCA, SD_RCA)

# Save the cleaned dataset internally
usethis::use_data(FirstNullAPCSSA, internal = TRUE, overwrite = TRUE)
