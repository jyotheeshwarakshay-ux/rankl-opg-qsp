## ============================================================
## Stage 1: Lemaire (2004) bone-remodeling model — untreated baseline
## Source: BIOMD0000000278_url.xml (BioModels)
## Three ODEs (responding osteoblasts R, active osteoblasts B,
## active osteoclasts C). Intervention inputs I_L, I_O, I_P held at 0.
## ============================================================

library(rxode2)
library(ggplot2)

## ---- Model -------------------------------------------------
## Rate laws and algebraic relations transcribed directly from the
## <listOfRules> section of the SBML file (rateRule / assignmentRule).
mod <- rxode2({
  ## assignment rules (algebraic, recomputed at every step)
  D_B   <- f0 * d_B
  Phi_C <- (C + f0 * C_s) / (C + C_s)

  Pbar  <- I_P / k_P
  P_O   <- S_P / k_P
  P_S   <- k6 / k5
  Phi_P <- (Pbar + P_O) / (Pbar + P_S)

  Phi_L <- ((k3 / k4) * K_L_P * Phi_P * B /
              (1 + (k3 * K / k4) + (k1 / (k2 * kO)) * (I_O + K_O_P * R / Phi_P))) /
           (1 + I_L / r_L)

  ## rate rules (state derivatives)
  d/dt(R) <- D_R * Phi_C - D_B * R / Phi_C
  d/dt(B) <- D_B * R / Phi_C - k_B * B
  d/dt(C) <- D_C * Phi_L - D_A * Phi_C * C
})

## ---- Parameters (verbatim from <listOfParameters>) ---------
pars <- c(
  C_s   = 0.005,      # OC-related resorption saturation (pM)
  D_A   = 0.7,        # OC apoptosis rate related term (1/day)
  d_B   = 0.7,        # OB differentiation rate (1/day)
  D_C   = 0.0021,     # OC differentiation rate (1/day)
  D_R   = 0.0007,     # responding OB differentiation rate (1/day)
  f0    = 0.05,       # fraction of C_s contributing at baseline
  I_L   = 0,          # exogenous RANKL input -> 0 = untreated baseline
  I_O   = 0,          # exogenous OPG input   -> 0 = untreated baseline
  I_P   = 0,          # exogenous PTH input    -> 0 = untreated baseline
  K     = 10,         # RANKL/OPG binding ratio constant
  k1    = 0.01,       # OPG production rate
  k2    = 10,         # OPG degradation rate
  k3    = 0.00058,    # RANKL-RANK binding rate
  k4    = 0.017,      # RANKL-RANK dissociation rate
  k5    = 0.02,       # PTH receptor degradation rate
  k6    = 3,          # PTH receptor production rate
  k_B   = 0.189,      # active OB apoptosis rate (1/day)
  K_L_P = 3000000,    # RANKL production rate constant
  kO    = 0.35,       # OPG production rate constant
  K_O_P = 200000,     # OPG production rate constant (PTH-dependent)
  k_P   = 86,         # PTH clearance rate (1/day)
  r_L   = 1000,       # RANKL saturation constant for I_L
  S_P   = 250         # PTH synthesis rate
)

## ---- Initial conditions (from <listOfSpecies>) --------------
inits <- c(
  R = 0.0007734,   # Responding osteoblasts
  B = 0.0007282,   # Active osteoblasts
  C = 0.0009127    # Active osteoclasts
)

## ---- Simulate long enough to confirm/reach steady state -----
## Time unit = days (consistent with rate constants, e.g. k_B ~ 1/day).
## 3000 days (~8 years) is far more than the system's relaxation time.
times <- seq(0, 3000, length.out = 2000)
ev <- et(times)

sim <- rxSolve(mod, params = pars, events = ev, inits = inits)

## ---- Plot ----------------------------------------------------
df <- as.data.frame(sim)
df_long <- reshape(
  df[, c("time", "R", "B", "C")],
  varying = c("R", "B", "C"),
  v.names = "value",
  timevar = "state",
  times = c("R", "B", "C"),
  direction = "long"
)

p <- ggplot(df_long, aes(time, value, color = state)) +
  geom_line(linewidth = 1) +
  labs(
    title = "Lemaire (2004) bone-remodeling model — untreated baseline",
    x = "Time (days)",
    y = "Cell population (pM)",
    color = "State"
  ) +
  theme_minimal()

print(p)
ggsave("stage1_baseline.png", p, width = 7, height = 5, dpi = 150)

## ---- Report steady-state values ------------------------------
final <- tail(df, 1)
cat("\nFinal (steady-state) values at t =", final$time, "days:\n")
cat(sprintf("  R = %.7f\n", final$R))
cat(sprintf("  B = %.7f\n", final$B))
cat(sprintf("  C = %.7f\n", final$C))
