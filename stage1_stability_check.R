## ============================================================
## Stage 1 stability check: is the baseline steady state a true
## attractor, or just a fixed point we happened to start on?
##
## Method: perturb the initial conditions of R, B, C by +20%
## (individually, and all three together), simulate forward with
## the same untreated model (I_L = I_O = I_P = 0), and confirm the
## trajectories relax back to the published baseline steady state.
## ============================================================

library(rxode2)
library(ggplot2)

## ---- Model (identical to Stage 1) ----------------------------
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

pars <- c(
  C_s = 0.005, D_A = 0.7, d_B = 0.7, D_C = 0.0021, D_R = 0.0007, f0 = 0.05,
  I_L = 0, I_O = 0, I_P = 0, K = 10,
  k1 = 0.01, k2 = 10, k3 = 0.00058, k4 = 0.017, k5 = 0.02, k6 = 3,
  k_B = 0.189, K_L_P = 3000000, kO = 0.35, K_O_P = 200000, k_P = 86,
  r_L = 1000, S_P = 250
)

## ---- Baseline steady state (from the SBML file) --------------
R0 <- 0.0007734
B0 <- 0.0007282
C0 <- 0.0009127

## ---- Perturbation scenarios: +20% on each state, alone and together
scenarios <- list(
  "Baseline (no perturbation)" = c(R = R0,       B = B0,       C = C0),
  "R +20%"                     = c(R = R0 * 1.2, B = B0,       C = C0),
  "B +20%"                     = c(R = R0,       B = B0 * 1.2, C = C0),
  "C +20%"                     = c(R = R0,       B = B0,       C = C0 * 1.2),
  "R, B, C all +20%"           = c(R = R0 * 1.2, B = B0 * 1.2, C = C0 * 1.2)
)

## ---- Simulate each scenario -----------------------------------
## 400 days is well beyond the ~30-70 day relaxation time observed
## for this system, leaving a long flat tail to confirm convergence.
times <- seq(0, 400, length.out = 4000)
ev <- et(times)

sim_list <- lapply(names(scenarios), function(nm) {
  d <- as.data.frame(rxSolve(mod, params = pars, events = ev, inits = scenarios[[nm]]))
  d$scenario <- nm
  d
})
df <- do.call(rbind, sim_list)
df$scenario <- factor(df$scenario, levels = names(scenarios))

## ---- Reshape to long format for faceted plotting ---------------
long <- do.call(rbind, lapply(c("R", "B", "C"), function(st) {
  data.frame(time = df$time, scenario = df$scenario, state = st, value = df[[st]])
}))
long$state <- factor(long$state, levels = c("R", "B", "C"))

ref <- data.frame(
  state = factor(c("R", "B", "C"), levels = c("R", "B", "C")),
  value = c(R0, B0, C0)
)

## ---- Plot -------------------------------------------------------
state_labels <- c(
  R = "R: Responding osteoblasts",
  B = "B: Active osteoblasts",
  C = "C: Active osteoclasts"
)

p <- ggplot(long, aes(time, value, color = scenario)) +
  geom_line(linewidth = 0.7) +
  geom_hline(data = ref, aes(yintercept = value), linetype = "dashed", color = "grey30") +
  facet_wrap(~state, scales = "free_y", ncol = 1, labeller = as_labeller(state_labels)) +
  labs(
    title = "Stability check: return to baseline after +20% perturbation of initial conditions",
    subtitle = "Dashed line = published baseline steady state",
    x = "Time (days)",
    y = "Cell population (pM)",
    color = "Scenario"
  ) +
  theme_minimal()

print(p)
ggsave("stage1_stability_check.png", p, width = 8, height = 9, dpi = 150)

## ---- Quantify the return to equilibrium --------------------------
final_time <- max(df$time)
final_df <- df[df$time == final_time, ]

cat("\n=== Final state at t =", final_time, "days ===\n")
for (i in seq_len(nrow(final_df))) {
  row <- final_df[i, ]
  cat(sprintf(
    "%-28s  R=%.7f (%+6.3f%%)  B=%.7f (%+6.3f%%)  C=%.7f (%+6.3f%%)\n",
    as.character(row$scenario),
    row$R, 100 * (row$R - R0) / R0,
    row$B, 100 * (row$B - B0) / B0,
    row$C, 100 * (row$C - C0) / C0
  ))
}

## Settling time: earliest time after which the relative deviation from
## the baseline steady state stays below a threshold for the rest of
## the simulated horizon (NA = did not settle within the simulated window).
settle_time <- function(time, value, ss, thresh) {
  rel <- abs(value - ss) / ss
  above <- which(rel >= thresh)
  if (length(above) == 0) return(0)
  last <- max(above)
  if (last == length(time)) return(NA_real_)
  time[last + 1]
}

cat("\n=== Settling time (days) to within 1% / 0.1% of baseline SS ===\n")
for (nm in names(scenarios)) {
  d <- df[df$scenario == nm, ]
  st_R_1  <- settle_time(d$time, d$R, R0, 0.01)
  st_B_1  <- settle_time(d$time, d$B, B0, 0.01)
  st_C_1  <- settle_time(d$time, d$C, C0, 0.01)
  st_R_01 <- settle_time(d$time, d$R, R0, 0.001)
  st_B_01 <- settle_time(d$time, d$B, B0, 0.001)
  st_C_01 <- settle_time(d$time, d$C, C0, 0.001)
  cat(sprintf(
    "%-28s  1%%: R=%-5s B=%-5s C=%-5s   |   0.1%%: R=%-5s B=%-5s C=%-5s\n",
    nm, st_R_1, st_B_1, st_C_1, st_R_01, st_B_01, st_C_01
  ))
}
