## ============================================================
## Stage 4 (continued): dosed-response across the virtual population,
## and a variance decomposition of what actually drives the spread.
##
## This is an illustrative virtual population -- a wide physiological
## range of parameter values chosen for interpretability, NOT calibrated
## to real patient data. Same 100-patient population and CV=30% already
## validated (0/100 excluded by the mandatory untreated-steady-state
## filter; see the prior diagnostic run for why that's a real finding,
## not a weak filter).
##
## PART 1: dose all 100 plausible patients (60 mg SC q6mo x 3, as in
## Stage 3) and characterize the response distribution -- peak
## osteoclast (C) suppression, and how completely C recovers before
## each re-dose.
##
## PART 2 (the key analysis): decompose where the response variability
## comes from. Three scenarios, built from the SAME underlying random
## draws (same seed, same per-patient z-scores), so only the source of
## variation differs:
##   - Full:      PK (ka, Vc, CL) and bone (D_R, D_C, k_B) both vary
##   - Bone-only: bone parameters vary per patient; PK fixed at typical
##   - PK-only:   PK parameters vary per patient; bone fixed at typical
## Comparing the response spread across these three isolates whether
## baseline-equilibrium differences (bone) or drug-exposure differences
## (PK) dominate the variability in denosumab response.
## ============================================================

library(rxode2)
library(ggplot2)

set.seed(42)

## ---- Combined model (identical to stage3_pk_coupling.R / stage4_virtual_population.R) ----
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

## ---- Sample the population: same draws underlie all three scenarios -----
cv    <- 0.30
omega <- sqrt(log(1 + cv^2))
n_candidates <- 100

varied_names <- c("ka", "Vc", "CL", "D_R", "D_C", "k_B")
typical_vals <- c(
  ka = pk_ka_h * 24, Vc = pk_Vc, CL = pk_CL_h * 24,
  D_R = unname(bone_typical["D_R"]), D_C = unname(bone_typical["D_C"]), k_B = unname(bone_typical["k_B"])
)

z <- matrix(
  rnorm(n_candidates * length(varied_names), 0, 1),
  nrow = n_candidates, ncol = length(varied_names),
  dimnames = list(NULL, varied_names)
)
eta <- z * omega
patient_mult <- as.data.frame(sweep(exp(eta), 2, typical_vals, `*`))
patient_mult$id <- seq_len(n_candidates)

## Three scenarios, same underlying patients, different parameters "frozen"
patient_df_full <- patient_mult

patient_df_bone_only <- patient_mult
patient_df_bone_only$ka <- typical_vals[["ka"]]
patient_df_bone_only$Vc <- typical_vals[["Vc"]]
patient_df_bone_only$CL <- typical_vals[["CL"]]

patient_df_pk_only <- patient_mult
patient_df_pk_only$D_R <- typical_vals[["D_R"]]
patient_df_pk_only$D_C <- typical_vals[["D_C"]]
patient_df_pk_only$k_B <- typical_vals[["k_B"]]

## ---- Helper functions -----------------------------------------------------
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

check_untreated_ss <- function(pars_vec, t_end = 3000) {
  inits0 <- c(depot = 0, Ctot = 0, Cp = 0, Rtot = pk_R0,
              R = unname(baseline_inits_bone["R"]), B = unname(baseline_inits_bone["B"]), C = unname(baseline_inits_bone["C"]))
  ev0 <- eventTable()
  ev0$add.sampling(seq(0, t_end, length.out = 400))

  sol0 <- tryCatch(as.data.frame(rxSolve(mod4, pars_vec, ev0, inits = inits0)), error = function(e) NULL)
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
  if (any(rel_change > 0.001)) return(list(pass = FALSE, reason = "not_converged"))
  list(pass = TRUE, reason = "ok", R_ss = final[1], B_ss = final[2], C_ss = final[3])
}

MW_kDa <- 147; dose_mg <- 60
dose_nmol <- (dose_mg * 1e6) / (MW_kDa * 1e3)
month_d <- 30.44
dose_times_d <- c(0, 6, 12) * month_d
times <- seq(0, 700, length.out = 1400)

## Runs the full filter -> dose -> metrics pipeline for one parameter scenario.
run_scenario <- function(patient_df, label) {
  filt <- lapply(seq_len(nrow(patient_df)), function(i) {
    pars_vec <- build_pars(patient_df[i, ])
    res <- check_untreated_ss(pars_vec)
    res$id <- i
    res
  })
  pass_flags <- sapply(filt, function(r) r$pass)
  passed_ids <- which(pass_flags)

  dosed <- lapply(passed_ids, function(i) {
    row <- patient_df[i, ]; ss <- filt[[i]]
    pars_vec <- build_pars(row)
    inits <- c(depot = 0, Ctot = 0, Cp = 0, Rtot = pk_R0, R = ss$R_ss, B = ss$B_ss, C = ss$C_ss)
    ev <- eventTable()
    for (dt in dose_times_d) ev$add.dosing(dose = dose_nmol, start.time = dt, dosing.to = "depot")
    ev$add.sampling(times)
    d <- as.data.frame(rxSolve(mod4, pars_vec, ev, inits = inits))
    d$id <- i; d$C_ss_baseline <- ss$C_ss; d$R_ss_baseline <- ss$R_ss; d$B_ss_baseline <- ss$B_ss
    d
  })
  pop_df <- do.call(rbind, dosed)

  metrics <- do.call(rbind, lapply(passed_ids, function(i) {
    d <- pop_df[pop_df$id == i, ]
    base_C <- unique(d$C_ss_baseline)
    peak_supp <- 100 * (base_C - min(d$C[d$time <= 30])) / base_C
    idx_trough <- which.min(abs(d$time - (dose_times_d[2] - 0.1)))
    trough_recovery <- 100 * d$C[idx_trough] / base_C
    data.frame(id = i, peak_suppression_pct = peak_supp, trough_recovery_pct = trough_recovery)
  }))

  cat(sprintf("[%s] passed %d / %d\n", label, length(passed_ids), nrow(patient_df)))

  list(label = label, pop_df = pop_df, metrics = metrics, n_passed = length(passed_ids))
}

cat("Running dosed-response pipeline for 3 scenarios (Full, Bone-only, PK-only)...\n\n")
res_full      <- run_scenario(patient_df_full,      "Full (PK + bone vary)")
res_bone_only <- run_scenario(patient_df_bone_only, "Bone-only (PK fixed)")
res_pk_only   <- run_scenario(patient_df_pk_only,   "PK-only (bone fixed)")

## =============================================================================
## PART 1: response distribution across the population (Full scenario)
## =============================================================================
cat("\n=== PART 1: dosed-response distribution across the population (Full scenario) ===\n")
cat("\nPeak osteoclast (C) suppression after dose 1:\n")
print(summary(res_full$metrics$peak_suppression_pct))
cat(sprintf("  SD = %.2f%%\n", sd(res_full$metrics$peak_suppression_pct)))

cat("\nTrough recovery of C (%% of baseline restored) just before dose 2:\n")
print(summary(res_full$metrics$trough_recovery_pct))
cat(sprintf("  SD = %.2f%%\n", sd(res_full$metrics$trough_recovery_pct)))

## ---- Reference "typical" patient (all etas = 0) --------------------------
typical_row <- as.data.frame(as.list(typical_vals))
typical_pars <- build_pars(typical_row)
typical_inits <- c(depot = 0, Ctot = 0, Cp = 0, Rtot = pk_R0,
                    R = unname(baseline_inits_bone["R"]), B = unname(baseline_inits_bone["B"]), C = unname(baseline_inits_bone["C"]))
ev_typ <- eventTable()
for (dt in dose_times_d) ev_typ$add.dosing(dose = dose_nmol, start.time = dt, dosing.to = "depot")
ev_typ$add.sampling(times)
typical_sol <- as.data.frame(rxSolve(mod4, typical_pars, ev_typ, inits = typical_inits))

## ---- Plot 1: population ribbon (median + 10th-90th pct) over time --------
summarize_state <- function(df, state) {
  agg <- aggregate(df[[state]], by = list(time = df$time), FUN = function(x) {
    c(median = median(x), lo = quantile(x, 0.10), hi = quantile(x, 0.90))
  })
  data.frame(time = agg$time, median = agg$x[, 1], lo = agg$x[, 2], hi = agg$x[, 3], state = state)
}
ribbon_df <- do.call(rbind, lapply(c("R", "B", "C"), summarize_state, df = res_full$pop_df))
ribbon_df$state <- factor(ribbon_df$state, levels = c("R", "B", "C"))

typ_long <- do.call(rbind, lapply(c("R", "B", "C"), function(s) {
  data.frame(time = typical_sol$time, value = typical_sol[[s]], state = s)
}))
typ_long$state <- factor(typ_long$state, levels = c("R", "B", "C"))

state_labels <- c(R = "R: Responding osteoblasts", B = "B: Active osteoblasts", C = "C: Active osteoclasts")

p1 <- ggplot(ribbon_df, aes(time, median)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.25) +
  geom_line(color = "steelblue", linewidth = 0.9) +
  geom_line(data = typ_long, aes(time, value), color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = dose_times_d, linetype = "dotted", color = "grey40") +
  facet_wrap(~state, scales = "free_y", ncol = 1, labeller = as_labeller(state_labels)) +
  labs(
    title = paste0("Stage 4: dosed-response across the population (n = ", res_full$n_passed, ")"),
    subtitle = "Shaded band = 10th-90th percentile, solid = median, dashed = typical patient. Illustrative population, not calibrated.",
    x = "Time (days)", y = NULL
  ) +
  theme_minimal()
print(p1)
ggsave("stage4_dosed_population_response.png", p1, width = 8, height = 9, dpi = 150)

## ---- Plot 2: peak suppression and trough recovery histograms -------------
p2a <- ggplot(res_full$metrics, aes(peak_suppression_pct)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "white") +
  geom_vline(xintercept = median(res_full$metrics$peak_suppression_pct), linetype = "dashed") +
  labs(title = "Peak osteoclast (C) suppression after dose 1",
       subtitle = paste0("median = ", round(median(res_full$metrics$peak_suppression_pct), 1), "%"),
       x = "Peak % suppression of C vs. own baseline", y = "Patients")

p2b <- ggplot(res_full$metrics, aes(trough_recovery_pct)) +
  geom_histogram(bins = 20, fill = "darkorange", color = "white") +
  geom_vline(xintercept = median(res_full$metrics$trough_recovery_pct), linetype = "dashed") +
  labs(title = "C recovery just before dose 2",
       subtitle = paste0("median = ", round(median(res_full$metrics$trough_recovery_pct), 1), "% of baseline"),
       x = "% of baseline C restored", y = "Patients")

print(p2a); print(p2b)
ggsave("stage4_dosed_peak_suppression.png", p2a, width = 6, height = 5, dpi = 150)
ggsave("stage4_dosed_trough_recovery.png", p2b, width = 6, height = 5, dpi = 150)

## =============================================================================
## PART 2: variance decomposition -- what drives the spread?
## =============================================================================
cat("\n=== PART 2: variance decomposition -- bone vs. PK contribution to spread ===\n")

decomp_df <- rbind(
  data.frame(scenario = "Full (both vary)", res_full$metrics),
  data.frame(scenario = "Bone-only (PK fixed)", res_bone_only$metrics),
  data.frame(scenario = "PK-only (bone fixed)", res_pk_only$metrics)
)
decomp_df$scenario <- factor(decomp_df$scenario,
                              levels = c("Full (both vary)", "Bone-only (PK fixed)", "PK-only (bone fixed)"))

cat("\nPeak suppression %% -- summary by scenario:\n")
tapply(decomp_df$peak_suppression_pct, decomp_df$scenario, summary)

sd_table <- aggregate(peak_suppression_pct ~ scenario, decomp_df, sd)
names(sd_table) <- c("scenario", "sd_peak_suppression_pct")
var_table <- aggregate(peak_suppression_pct ~ scenario, decomp_df, var)
names(var_table) <- c("scenario", "var_peak_suppression_pct")

cat("\nSD of peak suppression %% by scenario:\n")
print(sd_table, row.names = FALSE)

sd_full      <- sd_table$sd_peak_suppression_pct[sd_table$scenario == "Full (both vary)"]
sd_bone_only <- sd_table$sd_peak_suppression_pct[sd_table$scenario == "Bone-only (PK fixed)"]
sd_pk_only   <- sd_table$sd_peak_suppression_pct[sd_table$scenario == "PK-only (bone fixed)"]

var_full      <- var_table$var_peak_suppression_pct[var_table$scenario == "Full (both vary)"]
var_bone_only <- var_table$var_peak_suppression_pct[var_table$scenario == "Bone-only (PK fixed)"]
var_pk_only   <- var_table$var_peak_suppression_pct[var_table$scenario == "PK-only (bone fixed)"]

cat(sprintf("\nBone-only SD is %.2fx the Full SD (%.2f%% vs %.2f%%)\n", sd_bone_only / sd_full, sd_bone_only, sd_full))
cat(sprintf("PK-only   SD is %.2fx the Full SD (%.2f%% vs %.2f%%)\n", sd_pk_only / sd_full, sd_pk_only, sd_full))
cat(sprintf("\nVariance share (approx, if independent): bone %.1f%%, PK %.1f%% of (var_bone+var_pk)\n",
            100 * var_bone_only / (var_bone_only + var_pk_only),
            100 * var_pk_only / (var_bone_only + var_pk_only)))
cat(sprintf("Additive check: var_bone_only + var_pk_only = %.2f  vs  var_full = %.2f  (ratio %.2f)\n",
            var_bone_only + var_pk_only, var_full, (var_bone_only + var_pk_only) / var_full))

if (sd_bone_only > sd_pk_only) {
  cat(sprintf("\n=> BONE-parameter differences dominate the response spread (%.1fx more spread than PK alone).\n",
              sd_bone_only / sd_pk_only))
} else {
  cat(sprintf("\n=> PK differences dominate the response spread (%.1fx more spread than bone alone).\n",
              sd_pk_only / sd_bone_only))
}

## ---- Plot 3: variance decomposition comparison ---------------------------
p3 <- ggplot(decomp_df, aes(scenario, peak_suppression_pct, fill = scenario)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 1) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 0.8) +
  labs(
    title = "Variance decomposition: what drives spread in peak C suppression?",
    subtitle = "Same 100 patients' underlying random draws; only which parameter group varies differs across panels",
    x = NULL, y = "Peak % suppression of C after dose 1"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
print(p3)
ggsave("stage4_variance_decomposition.png", p3, width = 8, height = 6, dpi = 150)

cat("\n=== Stage 4 dosed-response + variance decomposition complete ===\n")
