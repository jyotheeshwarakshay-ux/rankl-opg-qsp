## ============================================================
## Stage 2: Denosumab (OPG-like RANKL blocker) dose-response
##
## Denosumab mimics the natural function of OPG: it neutralizes RANKL
## before it can bind RANK on osteoclast precursors. In the Lemaire (2004)
## model this is exactly what the existing I_O intervention input
## represents, so denosumab is modeled here purely as I_O > 0.
##
## No core equations or parameter values are changed from stage1_baseline.R.
## Only I_O is driven, across three scenarios simulated on the same axes:
##   1. Baseline        I_O = 0
##   2. Moderate dose    I_O = 2000
##   3. Strong dose      I_O = 16000
## ============================================================

library(rxode2)
library(ggplot2)

## ---- Model (identical to Stage 1 — no equations/parameters changed) ----
mod <- rxode2({
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

## ---- Base parameter set (identical to Stage 1; I_O varied per scenario) ----
base_pars <- c(
  C_s = 0.005, D_A = 0.7, d_B = 0.7, D_C = 0.0021, D_R = 0.0007, f0 = 0.05,
  I_L = 0, I_O = 0, I_P = 0, K = 10,
  k1 = 0.01, k2 = 10, k3 = 0.00058, k4 = 0.017, k5 = 0.02, k6 = 3,
  k_B = 0.189, K_L_P = 3000000, kO = 0.35, K_O_P = 200000, k_P = 86,
  r_L = 1000, S_P = 250
)

## ---- Initial conditions: the untreated baseline steady state -------------
inits <- c(R = 0.0007734, B = 0.0007282, C = 0.0009127)

## ---- Denosumab-like OPG input scenarios -----------------------------------
## I_O doses chosen from the model's own steady-state dose-response
## (scanned numerically): 2000 gives a moderate ~10% drop in osteoclasts,
## 16000 gives a strong ~57% drop.
scenarios <- list(
  "Baseline (I_O = 0)"        = 0,
  "Moderate dose (I_O = 2000)"  = 2000,
  "Strong dose (I_O = 16000)"   = 16000
)

## ---- Simulate long enough to reach the new steady state -------------------
## New-dose transients settle within ~100 days; 1500 days leaves a long
## flat tail to confirm convergence.
times <- seq(0, 1500, length.out = 3000)
ev <- et(times)

sim_list <- lapply(names(scenarios), function(nm) {
  pars <- base_pars
  pars["I_O"] <- scenarios[[nm]]
  d <- as.data.frame(rxSolve(mod, params = pars, events = ev, inits = inits))
  d$scenario <- nm
  d
})
df <- do.call(rbind, sim_list)
df$scenario <- factor(df$scenario, levels = names(scenarios))

## ---- Reshape to long format for faceted plotting --------------------------
long <- do.call(rbind, lapply(c("R", "B", "C"), function(st) {
  data.frame(time = df$time, scenario = df$scenario, state = st, value = df[[st]])
}))
long$state <- factor(long$state, levels = c("R", "B", "C"))

state_labels <- c(
  R = "R: Responding osteoblasts",
  B = "B: Active osteoblasts",
  C = "C: Active osteoclasts"
)

p <- ggplot(long, aes(time, value, color = scenario)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~state, scales = "free_y", ncol = 1, labeller = as_labeller(state_labels)) +
  labs(
    title = "Denosumab (OPG-like RANKL blocker) dose-response",
    subtitle = "Driving I_O only — no other equations or parameters changed",
    x = "Time (days)",
    y = "Cell population (pM)",
    color = "Scenario"
  ) +
  theme_minimal()

print(p)
ggsave("stage2_denosumab.png", p, width = 8, height = 9, dpi = 150)

## ---- Report steady-state values -------------------------------------------
final_time <- max(df$time)
final_df <- df[df$time == final_time, ]

R0 <- inits["R"]; B0 <- inits["B"]; C0 <- inits["C"]

cat("\n=== Steady-state values at t =", final_time, "days ===\n\n")
tbl <- data.frame(
  Scenario = final_df$scenario,
  R = round(final_df$R, 7),
  B = round(final_df$B, 7),
  C = round(final_df$C, 7),
  R_pct_vs_baseline = round(100 * (final_df$R - R0) / R0, 2),
  B_pct_vs_baseline = round(100 * (final_df$B - B0) / B0, 2),
  C_pct_vs_baseline = round(100 * (final_df$C - C0) / C0, 2)
)
print(tbl, row.names = FALSE)
