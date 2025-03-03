# Load required packages
library(tidyverse)

# Read and process the dataset for APCSSA
FirstNullAPCSSM <- read_csv("data-raw/NullDistributionAPCSSM.csv") %>%
  mutate(
    SD_CRM = sqrt(V_CRM),
    SD_RCM = sqrt(V_RCM)
  ) %>%
  select(I, J, K, E_CRM, SD_CRM, E_RCM, SD_RCM)

