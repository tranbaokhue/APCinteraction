# Source individual dataset scripts to load them into the environment
source("data-raw/Null_APCSSA.R")
source("data-raw/Null_APCSSM.R")

# Save all internal datasets together in sysdata.rda
usethis::use_data(nullAPCXXA_summary, nullAPCXXM_summary, internal = TRUE, overwrite = TRUE)
