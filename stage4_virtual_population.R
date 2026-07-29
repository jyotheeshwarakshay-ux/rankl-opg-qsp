## ============================================================
## Stage 4: virtual population on the Stage 3 PK-coupled model
##
## Six parameters are varied across virtual patients, chosen for direct
## physiological interpretability (not an exhaustive sensitivity sweep):
##
##   PK side (drug exposure):
##     ka  - absorption rate:  governs how fast peak concentration is reached
##     Vc  - central volume:   governs peak concentration magnitude
##     CL  - clearance:        governs overall elimination / exposure (AUC)
##
##   Bone side (remodeling balance):
##     D_R - responding-osteoblast recruitment rate (formation drive)
##     D_C - osteoclast differentiation rate         (resorption drive)
##     k_B - active-osteoblast apoptosis rate         (osteoblast turnover)
##
## All six sampled independently, lognormal around the typical (Stage 3)
## value, CV = 30% -- a single modest, illustrative CV, NOT the literature
## IIVs from the denosumab-tmdd-qss project (some of which exceed 200%
## and are themselves flagged there as poorly identified).
##
## MANDATORY plausibility filter (before any dosing): simulate each
## candidate's UNTREATED model (I_O = 0, no dosing) from the population
## baseline state for 3000 days and require R, B, C to be finite,
## positive, and converged to a stable steady state throughout. Only
## patients that pass are dosed.
##
## This run: N = 100 candidates, to confirm the filter and outputs are
## sensible before scaling up.
## ============================================================

library(rxode2)
library(ggplot2)

set.seed(42)

## ---- Combined model (identical to stage3_pk_coupling.R) ----------------
mod4 <- rxode2({
  disc   <- (Ctot - Rtot - Kss)^2 + 4 * Kss * Ctot
  discP  <- max(disc, 0)
  C_deno <- max(0.5 * ((Ctot - Rtot - Kss) + sqrt(discP)), 0)
  RC     <- Ctot - C_deno

  d/dt(depot) <- -ka * depot
  d/dt(Ctot)  <- (ka * depot) / Vc - (CL / Vc) * C_deno - (Q / Vc) * C_deno +
                 (Q / Vp) * Cp - kint * RC
  d/dt(Cp)    <- (Q / Vc) * C_deno - (Q / Vp) * Cp
  d/dt(Rtot)  <- ksyn - kdeg * (Rtot - RC) - kint * RC

  I_O <- k_scale * C_deno

  D_B   <- f0 * d_B
  Phi_C <- (C + f0 * C_s) / (C + C_s)

  Pbar  <- I_P / k_P
  P_O   <- S_P / k_P
  P_S   <- k6 / k5
  Phi_P <- (Pbar + P_O) / (Pbar + P_S)

  Phi_L <- ((k3 / k4) * K_L_P * Phi_P * B /
              (1 + (k3 * K / k4) + (k1 / (k2 * kO)) * (I_O + K_O_P * R / Phi_P))) /
           (1 + I_L / r_L)

  d/dt(R) <- D_R * Phi_C - D_B * R / Phi_C
  d/dt(B) <- D_B * R / Phi_C - k_B * B
  d/dt(C) <- D_C * Phi_L - D_A * Phi_C * C
})

## ---- Typical (Stage 3) parameter values ---------------------------------
pk_ka_h   <- 0.0078; pk_Vc <- 1.58; pk_Vp <- 6.06; pk_CL_h <- 0.006
pk_Q_h    <- 0.20;   pk_kint_h <- 0.022; pk_Kss <- 1.56
pk_ksyn_h <- 0.01;   pk_R0 <- 15.23
pk_kdeg_h <- pk_ksyn_h / pk_R0

bone_typical <- c(
  C_s = 0.005, D_A = 0.7, d_B = 0.7, D_C = 0.0021, D_R = 0.0007, f0 = 0.05,
  I_L = 0, I_P = 0, K = 10,
  k1 = 0.01, k2 = 10, k3 = 0.00058, k4 = 0.017, k5 = 0.02, k6 = 3,
  k_B = 0.189, K_L_P = 3000000, kO = 0.35, K_O_P = 200000, k_P = 86,
  r_L = 1000, S_P = 250
)

k_scale <- 100

baseline_inits_bone <- c(R = 0.0007734, B = 0.0007282, C = 0.0009127)

## ---- Sample the virtual population --------------------------------------
cv    <- 0.30
omega <- sqrt(log(1 + cv^2))

n_candidates <- 100

varied_names <- c("ka", "Vc", "CL", "D_R", "D_C", "k_B")
typical_vals <- c(
  ka = pk_ka_h * 24, Vc = pk_Vc, CL = pk_CL_h * 24,
  D_R = unname(bone_typical["D_R"]), D_C = unname(bone_typical["D_C"]), k_B = unname(bone_typical["k_B"])
)

eta <- matrix(
  rnorm(n_candidates * length(varied_names), 0, omega),
  nrow = n_candidates, ncol = length(varied_names),
  dimnames = list(NULL, varied_names)
)

patient_mult <- sweep(exp(eta), 2, typical_vals, `*`)
patient_df <- as.data.frame(patient_mult)
patient_df$id <- seq_len(n_candidates)

cat("Sampled", n_candidates, "candidate virtual patients.\n")
cat("Varied parameters (lognormal, CV =", cv, "):", paste(varied_names, collapse = ", "), "\n\n")

## ---- Full parameter vector for one patient -------------------------------
build_pars <- function(row) {
  c(
    ka = row$ka, Vc = row$Vc, Vp = pk_Vp, CL = row$CL, Q = pk_Q_h * 24,
    kint = pk_kint_h * 24, Kss = pk_Kss, ksyn = pk_ksyn_h * 24, kdeg = pk_kdeg_h * 24,
    k_scale = k_scale,
    C_s = unname(bone_typical["C_s"]), D_A = unname(bone_typical["D_A"]), d_B = unname(bone_typical["d_B"]),
    D_C = row$D_C, D_R = row$D_R, f0 = unname(bone_typical["f0"]),
    I_L = 0, I_P = 0, K = unname(bone_typical["K"]),
    k1 = unname(bone_typical["k1"]), k2 = unname(bone_typical["k2"]), k3 = unname(bone_typical["k3"]),
    k4 = unname(bone_typical["k4"]), k5 = unname(bone_typical["k5"]), k6 = unname(bone_typical["k6"]),
    k_B = row$k_B, K_L_P = unname(bone_typical["K_L_P"]), kO = unname(bone_typical["kO"]),
    K_O_P = unname(bone_typical["K_O_P"]), k_P = unname(bone_typical["k_P"]), r_L = unname(bone_typical["r_L"]),
    S_P = unname(bone_typical["S_P"])
  )
}

## ---- MANDATORY plausibility filter: untreated steady state --------------
## Simulate 3000 days with NO dosing (I_O stays 0 throughout) from the
## population baseline state. Require finite, positive, and converged
## (final state within 0.1% of the state 500 days earlier).
check_untreated_ss <- function(pars_vec, t_end = 3000) {
  inits0 <- c(depot = 0, Ctot = 0, Cp = 0, Rtot = pk_R0,
              R = unname(baseline_inits_bone["R"]), B = unname(baseline_inits_bone["B"]), C = unname(baseline_inits_bone["C"]))
  ev0 <- eventTable()
  ev0$add.sampling(seq(0, t_end, length.out = 400))

  sol0 <- tryCatch(as.data.frame(rxSolve(mod4, pars_vec, ev0, inits = inits0)),
                    error = function(e) NULL)
  if (is.null(sol0)) return(list(pass = FALSE, reason = "solver_error"))

  if (any(!is.finite(sol0$R)) || any(!is.finite(sol0$B)) || any(!is.finite(sol0$C))) {
    return(list(pass = FALSE, reason = "non_finite"))
  }
  if (any(sol0$R <= 0) || any(sol0$B <= 0) || any(sol0$C <= 0)) {
    return(list(pass = FALSE, reason = "non_positive"))
  }

  n <- nrow(sol0)
  idx_check <- which.min(abs(sol0$time - t_end * (5 / 6)))
  final <- c(sol0$R[n], sol0$B[n], sol0$C[n])
  earlier <- c(sol0$R[idx_check], sol0$B[idx_check], sol0$C[idx_check])
  rel_change <- abs(final - earlier) / final

  if (any(rel_change > 0.001)) {
    return(list(pass = FALSE, reason = "not_converged"))
  }

  list(pass = TRUE, reason = "ok", R_ss = final[1], B_ss = final[2], C_ss = final[3])
}

cat("Running mandatory untreated-steady-state filter on all", n_candidates, "candidates...\n")

filter_results <- lapply(seq_len(n_candidates), function(i) {
  pars_vec <- build_pars(patient_df[i, ])
  res <- check_untreated_ss(pars_vec)
  res$id <- i
  res
})

pass_flags <- sapply(filter_results, function(r) r$pass)
reasons    <- sapply(filter_results, function(r) r$reason)

cat("\n=== Plausibility filter report ===\n")
cat("Candidates generated:", n_candidates, "\n")
cat("Passed               :", sum(pass_flags), "\n")
cat("Excluded              :", sum(!pass_flags), "\n")
if (any(!pass_flags)) {
  cat("\nExclusion reasons:\n")
  print(table(reasons[!pass_flags]))
} else {
  cat("(no exclusions)\n")
}

## ---- Dose only the patients that passed ----------------------------------
passed_ids <- which(pass_flags)
cat("\nDosing the", length(passed_ids), "patients that passed the filter (60 mg SC q6mo x 3)...\n")

MW_kDa <- 147; dose_mg <- 60
dose_nmol <- (dose_mg * 1e6) / (MW_kDa * 1e3)
month_d <- 30.44
dose_times_d <- c(0, 6, 12) * month_d

times <- seq(0, 700, length.out = 1400)

dosed_results <- lapply(passed_ids, function(i) {
  row <- patient_df[i, ]
  ss  <- filter_results[[i]]
  pars_vec <- build_pars(row)

  inits <- c(depot = 0, Ctot = 0, Cp = 0, Rtot = pk_R0,
             R = ss$R_ss, B = ss$B_ss, C = ss$C_ss)

  ev <- eventTable()
  for (dt in dose_times_d) ev$add.dosing(dose = dose_nmol, start.time = dt, dosing.to = "depot")
  ev$add.sampling(times)

  d <- as.data.frame(rxSolve(mod4, pars_vec, ev, inits = inits))
  d$id <- i
  d$C_ss_baseline <- ss$C_ss
  d$R_ss_baseline <- ss$R_ss
  d$B_ss_baseline <- ss$B_ss
  d
})

pop_df <- do.call(rbind, dosed_results)

## ---- Reference "typical" patient (all etas = 0, i.e. Stage 3) -----------
typical_row <- as.data.frame(as.list(typical_vals))
typical_pars <- build_pars(typical_row)
typical_inits <- c(depot = 0, Ctot = 0, Cp = 0, Rtot = pk_R0,
                    R = unname(baseline_inits_bone["R"]), B = unname(baseline_inits_bone["B"]), C = unname(baseline_inits_bone["C"]))
ev_typ <- eventTable()
for (dt in dose_times_d) ev_typ$add.dosing(dose = dose_nmol, start.time = dt, dosing.to = "depot")
ev_typ$add.sampling(times)
typical_sol <- as.data.frame(rxSolve(mod4, typical_pars, ev_typ, inits = typical_inits))

## ---- Summary plot 1: population median + 10th-90th percentile ribbon ----
summarize_state <- function(df, state) {
  agg <- aggregate(df[[state]], by = list(time = df$time), FUN = function(x) {
    c(median = median(x), lo = quantile(x, 0.10), hi = quantile(x, 0.90))
  })
  data.frame(time = agg$time, median = agg$x[, 1], lo = agg$x[, 2], hi = agg$x[, 3], state = state)
}

ribbon_df <- do.call(rbind, lapply(c("R", "B", "C"), summarize_state, df = pop_df))
ribbon_df$state <- factor(ribbon_df$state, levels = c("R", "B", "C"))

typ_long <- do.call(rbind, lapply(c("R", "B", "C"), function(s) {
  data.frame(time = typical_sol$time, value = typical_sol[[s]], state = s)
}))
typ_long$state <- factor(typ_long$state, levels = c("R", "B", "C"))

state_labels <- c(
  R = "R: Responding osteoblasts",
  B = "B: Active osteoblasts",
  C = "C: Active osteoclasts"
)

p1 <- ggplot(ribbon_df, aes(time, median)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.25) +
  geom_line(color = "steelblue", linewidth = 0.9) +
  geom_line(data = typ_long, aes(time, value), color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = dose_times_d, linetype = "dotted", color = "grey40") +
  facet_wrap(~state, scales = "free_y", ncol = 1, labeller = as_labeller(state_labels)) +
  labs(
    title = paste0("Stage 4: virtual population bone response (n = ", length(passed_ids), " plausible patients)"),
    subtitle = "Shaded band = 10th-90th percentile, solid = population median, dashed = Stage 3 typical patient. Dotted lines = doses.",
    x = "Time (days)", y = NULL
  ) +
  theme_minimal()

print(p1)
ggsave("stage4_population_response.png", p1, width = 8, height = 9, dpi = 150)

## ---- Summary plot 2: distribution of peak C suppression (first dose) ----
peak_suppression <- sapply(passed_ids, function(i) {
  d <- pop_df[pop_df$id == i, ]
  base_C <- unique(d$C_ss_baseline)
  min_C_first_dose <- min(d$C[d$time <= 30])
  100 * (base_C - min_C_first_dose) / base_C
})

supp_df <- data.frame(id = passed_ids, peak_suppression_pct = peak_suppression)

cat("\n=== Peak osteoclast (C) suppression after first dose, across", length(passed_ids), "patients ===\n")
print(summary(supp_df$peak_suppression_pct))

p2 <- ggplot(supp_df, aes(peak_suppression_pct)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "white") +
  geom_vline(xintercept = median(supp_df$peak_suppression_pct), linetype = "dashed", color = "black") +
  labs(
    title = "Stage 4: distribution of peak osteoclast (C) suppression after dose 1",
    subtitle = paste0("n = ", length(passed_ids), " plausible patients. Dashed line = median (",
                       round(median(supp_df$peak_suppression_pct), 1), "%)"),
    x = "Peak % suppression of C relative to patient's own untreated baseline",
    y = "Number of patients"
  ) +
  theme_minimal()

print(p2)
ggsave("stage4_suppression_distribution.png", p2, width = 7, height = 5, dpi = 150)

cat("\n=== Stage 4 (N =", n_candidates, "candidates) complete ===\n")
