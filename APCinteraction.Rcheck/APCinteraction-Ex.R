pkgname <- "APCinteraction"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "APCinteraction-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('APCinteraction')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("APCSSA")
### * APCSSA

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: APCSSA
### Title: APCSSA
### Aliases: APCSSA

### ** Examples

# Set the seed for reproducibility
set.seed(206)

# Parameters - number of levels for factors A and B, and replications
nA <- 3; nB <- 3; nrep <- 3

# Generate levels and create the full design
A_vals <- 1:nA
B_vals <- 1:nB
design <- expand.grid(A = A_vals, B = B_vals)
design <- design[rep(seq_len(nrow(design)), each = nrep), ]

# Compute response as interaction (product of numeric factor levels)
mu <- A_vals[design$A] * B_vals[design$B]
value <- mu + rnorm(length(mu))  # Add normal noise

# Assemble the data
data <- data.frame(
  value = value,
  A = factor(design$A),
  B = factor(design$B)
)

# Run the APCSSA test
APCSSA(value ~ A + B, data = data)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("APCSSA", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("APCSSM")
### * APCSSM

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: APCSSM
### Title: APCSSM
### Aliases: APCSSM

### ** Examples

# Set seed for reproducibility
set.seed(36)

# Parameters - number of levels for factos A and B, replications, and interaction effect
nA <- 3; nB <- 4; nrep <- 2; c <- 1.25

# Generate levels and create the full design
A_vals <- seq(-2, 2, length.out = nA)
B_vals <- seq(-1.5, 1.5, length.out = nB)
design <- expand.grid(A = 1:nA, B = 1:nB)
design <- design[rep(seq_len(nrow(design)), each = nrep), ]

# Create the specific interaction matrix: alternate ±c in top rows
specInt <- matrix(0, nA, nB)
specInt[1:2, 1:nB] <- matrix(c(c, -c, -c, c), nrow = 2, byrow = TRUE)

# Compute response
mu <- A_vals[design$A] + B_vals[design$B] + specInt[cbind(design$A, design$B)]
value <- mu + rt(length(mu), df = 1)

# Assemble the data
data <- data.frame(
  value = value,
  A = factor(design$A, labels = round(A_vals, 2)),
  B = factor(design$B, labels = round(B_vals, 2))
)

# With Cauchy errors, we opt for APCSSM to check for interaction
APCSSM(value ~ A + B, data = data)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("APCSSM", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("save_null")
### * save_null

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: save_null
### Title: Save null distribution to working directory
### Aliases: save_null

### ** Examples

## Not run: 
##D # After running simulation
##D sim_nullAPCSSA(2, 2, 2)
##D 
##D # Save to working directory
##D save_null("APCSSA", 2, 2, 2)
##D 
##D # Or save to custom location
##D save_null("APCSSA", 2, 2, 2, path = "~/my_nulls/")
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("save_null", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sim_nullAPCSSA")
### * sim_nullAPCSSA

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sim_nullAPCSSA
### Title: Simulate the null distribution for APCSSA
### Aliases: sim_nullAPCSSA

### ** Examples

## Not run: 
##D sim_nullAPCSSA(2, 2, 2, 5000)  # This should take only a few seconds
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sim_nullAPCSSA", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sim_nullAPCSSM")
### * sim_nullAPCSSM

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sim_nullAPCSSM
### Title: Simulate the null distribution for APCSSM
### Aliases: sim_nullAPCSSM

### ** Examples

## Not run: 
##D sim_nullAPCSSM(2, 2, 2, 5000)  # This should take only a few seconds
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sim_nullAPCSSM", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
