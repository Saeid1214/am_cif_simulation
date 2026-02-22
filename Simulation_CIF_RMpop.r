#CIF made using RM  parents (correlated parental liabilities)

rm(list = ls())
library(MASS)
library(survival)
library(polycor)



# ============================================================
# Helper functions
# ============================================================

# Truncated exponential sampler: Exp(rate) conditional on T <= max_time
rtruncexp <- function(n, rate, max_time) {
  u <- runif(n)
  -log(1 - u * (1 - exp(-rate * max_time))) / rate
}

# CIF/KM simulation in a Random Mating (RM) population
simulate_CIF_RM <- function(n, heritability, prevalence, followup_time, rate = 1/20) {

  mother_dx <- rbinom(n, 1, prevalence)
  father_dx <- rbinom(n, 1, prevalence)

  environmental_risk <- sqrt(1 - heritability) * rnorm(n)
  genetic_risk       <- sqrt(heritability) * (mother_dx + father_dx) / 2
  liability          <- genetic_risk + environmental_risk

  thr_child      <- qnorm(1 - prevalence)
  disease_status <- as.integer(liability > thr_child)

  onset_time <- rtruncexp(n, rate = rate, max_time = followup_time)
  time   <- ifelse(disease_status == 1, onset_time, followup_time)
  status <- disease_status

  dat <- data.frame(time = time, status = status, mother_dx = mother_dx, father_dx = father_dx)

  dat$group <- with(dat,
                    ifelse(mother_dx == 0 & father_dx == 0, "none",
                    ifelse(mother_dx == 1 & father_dx == 0, "mother",
                    ifelse(mother_dx == 0 & father_dx == 1, "father", "both"))))

  fits <- lapply(split(dat, dat$group), function(d) survfit(Surv(time, status) ~ 1, data = d))

  list(
    km_fit_none   = fits$none,
    km_fit_mother = fits$mother,
    km_fit_father = fits$father,
    km_fit_both   = fits$both
  )
}

# Build risk lookup table: risk(group, age) = 1 - S(age) from KM
make_risk_lookup <- function(cif_results, ages) {

  ages <- sort(unique(ages))

  risk_lookup <- expand.grid(
    group     = c("none", "mother", "father", "both"),
    child_age = ages,
    stringsAsFactors = FALSE
  )

  get_risk_from_fit <- function(fit, ages) {
    s <- summary(fit, times = ages)$surv
    1 - s
  }

  risk_lookup$risk <- NA_real_
  risk_lookup$risk[risk_lookup$group == "none"]   <- get_risk_from_fit(cif_results$km_fit_none,   ages)
  risk_lookup$risk[risk_lookup$group == "mother"] <- get_risk_from_fit(cif_results$km_fit_mother, ages)
  risk_lookup$risk[risk_lookup$group == "father"] <- get_risk_from_fit(cif_results$km_fit_father, ages)
  risk_lookup$risk[risk_lookup$group == "both"]   <- get_risk_from_fit(cif_results$km_fit_both,   ages)

  risk_lookup
}

# Robust polychoric summary helper (safe for rare/degenerate binary outcomes)
polychor_summary <- function(x, y) {

  # If either variable is constant, correlation is undefined
  if (length(unique(x[!is.na(x)])) < 2 || length(unique(y[!is.na(y)])) < 2) {
    return(data.frame(rho = NA_real_, se = NA_real_, prev = mean(y, na.rm = TRUE)))
  }

  pc <- tryCatch(
    polycor::polychor(x, y, std.err = FALSE),
    error = function(e) NULL
  )

  if (is.null(pc)) {
    return(data.frame(rho = NA_real_, se = NA_real_, prev = mean(y, na.rm = TRUE)))
  }

  # polychor can return either:
  # - a numeric correlation (atomic)
  # - a list with $rho (and maybe $var)
  if (is.list(pc)) {
    rho <- unname(pc$rho)
  } else {
    rho <- unname(pc)
  }

  if (!is.finite(rho)) rho <- NA_real_

  data.frame(
    rho  = rho,
    se   = NA_real_,                      # std.err disabled -> return NA
    prev = mean(y, na.rm = TRUE)
  )
}
# percentage difference: (x - y) / y * 100 ; y is reference
perc_diff <- function(x, y) {
  if (any(y == 0, na.rm = TRUE)) stop("Reference value y contains zero — cannot divide by zero.")
  ((x - y) / y) * 100
}

# ============================================================
# Simulation settings
# ============================================================

n_sim         <- 100000   # CIF generation sample size
h2            <- 0.50
prev          <- 0.10
followup_time <- 50
rate_onset    <- 1/20

am            <- 0.30
n_couples     <- 100000

N             <- 1000     # outer replicates
P             <- 100      # permutations per replicate

# Preallocate results (rho, se, prev percent-diff vs observed)
out_final <- matrix(NA_real_, nrow = N, ncol = 3)
colnames(out_final) <- c("rho", "se", "prev")

# ============================================================
# Main loop
# ============================================================

for (j in 1:N) {

  # 1) CIF from RM population
  cif_results <- simulate_CIF_RM(
    n = n_sim,
    heritability = h2,
    prevalence = prev,
    followup_time = followup_time,
    rate = rate_onset
  )

  # 2) Generate couples with AM on parental liability
  Sigma <- matrix(c(1, am,
                    am, 1), nrow = 2)
  liab <- MASS::mvrnorm(n = n_couples, mu = c(0, 0), Sigma = Sigma)

  thr_parent <- qnorm(1 - prev)

  couples <- data.frame(
    couple_id = 1:n_couples,
    mother_dx = as.integer(liab[, 1] > thr_parent),
    father_dx = as.integer(liab[, 2] > thr_parent)
  )

  # Child observed ages (1..30)
  couples$child_age <- sample(1:30, size = nrow(couples), replace = TRUE)

  # Parental diagnostic group
  couples$group <- with(couples,
                        ifelse(mother_dx == 0 & father_dx == 0, "none",
                        ifelse(mother_dx == 1 & father_dx == 0, "mother",
                        ifelse(mother_dx == 0 & father_dx == 1, "father", "both"))))

  # 3) Compute child risks (lookup) and generate children
  risk_lookup <- make_risk_lookup(cif_results, couples$child_age)

  couples <- merge(
    couples, risk_lookup,
    by = c("group", "child_age"),
    all.x = TRUE, sort = FALSE
  )
  names(couples)[names(couples) == "risk"] <- "child_risk"

  couples$child_dx <- rbinom(nrow(couples), 1, couples$child_risk)

  # Observed (non-permuted) summary
  obs <- polychor_summary(couples$mother_dx, couples$child_dx)
  obs_vec <- unlist(obs[1, c("rho", "se", "prev")])

  # 4) Permutations: shuffle fathers, resimulate children, summarize
  perm_out <- matrix(NA_real_, nrow = P, ncol = 3)
  colnames(perm_out) <- c("rho", "se", "prev")

  for (i in 1:P) {

    Rcouples <- couples
    Rcouples$Rfather_dx <- sample(Rcouples$father_dx)

    Rcouples$Rgroup <- with(Rcouples,
                            ifelse(mother_dx == 0 & Rfather_dx == 0, "none",
                            ifelse(mother_dx == 1 & Rfather_dx == 0, "mother",
                            ifelse(mother_dx == 0 & Rfather_dx == 1, "father", "both"))))

    Ncouples <- merge(
      Rcouples,
      risk_lookup,
      by.x = c("Rgroup", "child_age"),
      by.y = c("group", "child_age"),
      all.x = TRUE,
      sort = FALSE
    )
    names(Ncouples)[names(Ncouples) == "risk"] <- "Rchild_risk"

    Ncouples$Rchild_dx <- rbinom(nrow(Ncouples), 1, Ncouples$Rchild_risk)

    p1 <- polychor_summary(Ncouples$mother_dx, Ncouples$Rchild_dx)
    perm_out[i, ] <- unlist(p1[1, c("rho", "se", "prev")])
  }

  perm_mean <- colMeans(perm_out, na.rm = TRUE)

  # 5) Percent difference: permuted vs observed (reference = observed)
  out_final[j, ] <- perc_diff(perm_mean, obs_vec)
}

# ============================================================
# Output
# ============================================================

out_final_df <- as.data.frame(out_final)

print(head(out_final_df))
print(colMeans(out_final_df, na.rm = TRUE))

out_file <- paste0("simulation_RM1000_prev_", prev, "_h2_", h2, ".txt")
write.table(out_final_df, out_file, col.names = TRUE, row.names = FALSE, quote = FALSE)

