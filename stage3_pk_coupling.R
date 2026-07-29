## ============================================================
## Stage 3: couple denosumab PK to the Lemaire (2004) bone model
##
## Denosumab is a monoclonal antibody that neutralizes RANKL, mimicking
## the natural decoy function of OPG. Stage 2 drove the bone model's
## existing I_O (exogenous OPG) input as a fixed step input. Stage 3
## replaces that fixed input with a dynamic one: I_O(t) is driven by the
## simulated denosumab concentration from a separate PK model, so the
## drug's own pharmacokinetics (absorption, distribution, elimination,
## target-mediated disposition) now shape the time course of the bone
## response, rather than an arbitrary constant.
##
## Coupling rule (deliberately simple and explicit):
##   I_O(t) = k_scale * C_deno(t)
## where C_deno(t) is the free denosumab concentration and k_scale is a
## single illustrative scaling constant -- NOT calibrated to clinical
## bone-turnover-marker data. See the units/timescale check that preceded
## this script for the reasoning behind k_scale = 100.
##
## The bone sub-model's equations and parameters are UNCHANGED from
## stage1_baseline.R / stage2_denosumab.R. The only modification on the
## bone side is that I_O is now a computed quantity instead of a fixed
## parameter.
##
## Naming note: Project 3's PK model calls its free-drug concentration
## "C", which collides with the bone model's osteoclast state "C". The
## PK variable is renamed C_deno below -- this is a namespace fix only;
## no bone identifiers, equations, or parameters are touched.
##
## Time units: Project 3 is natively in hours; the bone model is in days.
## All PK rate constants are converted h -> day (x24) below so the whole
## system runs on a single day-based clock. Concentrations (nmol/L) and
## volumes (L) are unaffected by this conversion.
## ============================================================

library(rxode2)
library(ggplot2)

## ---- Combined model --------------------------------------------------
mod3 <- rxode2({
  ## ---- PK sub-model: two-compartment QSS-TMDD (Project 3, h -> day) ----
  disc   <- (Ctot - Rtot - Kss)^2 + 4 * Kss * Ctot
  discP  <- max(disc, 0)
  C_deno <- max(0.5 * ((Ctot - Rtot - Kss) + sqrt(discP)), 0)  # free denosumab, nmol/L
  RC     <- Ctot - C_deno

  d/dt(depot) <- -ka * depot
  d/dt(Ctot)  <- (ka * depot) / Vc - (CL / Vc) * C_deno - (Q / Vc) * C_deno +
                 (Q / Vp) * Cp - kint * RC
  d/dt(Cp)    <- (Q / Vc) * C_deno - (Q / Vp) * Cp
  d/dt(Rtot)  <- ksyn - kdeg * (Rtot - RC) - kint * RC

  ## ---- Coupling: denosumab concentration drives I_O ---------------------
  I_O <- k_scale * C_deno

  ## ---- Bone sub-model: Lemaire (2004) -- IDENTICAL to stage1_baseline.R --
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

## ---- Denosumab PK parameters ------------------------------------------
## Population-typical (IIV = 0) estimates, originally in per-hour units.
## Source: denosumab-tmdd-qss project (Table 3 ground truth), ultimately
## from Choi S, Park S, Jung J, Baek S, Lim H-S (2025). Population
## pharmacokinetics/pharmacodynamics analysis confirming biosimilarity of
## SB16 to reference denosumab. Front Pharmacol 16:1631034.
## doi:10.3389/fphar.2025.1631034
pk_ka_h   <- 0.0078   # 1/h    absorption rate constant
pk_Vc     <- 1.58     # L      central volume (apparent, Vc/F)
pk_Vp     <- 6.06     # L      peripheral volume (apparent, Vp/F)
pk_CL_h   <- 0.006    # L/h    clearance (apparent, CL/F)
pk_Q_h    <- 0.20     # L/h    inter-compartmental clearance
pk_kint_h <- 0.022    # 1/h    drug-target complex internalization rate
pk_Kss    <- 1.56     # nmol/L quasi-steady-state constant
pk_ksyn_h <- 0.01     # nmol/L/h  target (RANKL) synthesis rate
pk_R0     <- 15.23    # nmol/L baseline target concentration
pk_kdeg_h <- pk_ksyn_h / pk_R0   # 1/h, derived: ksyn = kdeg * R0 at steady state

pk_pars <- c(
  ka   = pk_ka_h * 24,     # -> 1/day
  Vc   = pk_Vc,            # L (concentration/volume terms unaffected by time unit)
  Vp   = pk_Vp,            # L
  CL   = pk_CL_h * 24,     # -> 1/day
  Q    = pk_Q_h * 24,      # -> 1/day
  kint = pk_kint_h * 24,   # -> 1/day
  Kss  = pk_Kss,           # nmol/L
  ksyn = pk_ksyn_h * 24,   # -> nmol/L/day
  kdeg = pk_kdeg_h * 24    # -> 1/day
)

## ---- Coupling constant --------------------------------------------------
## k_scale = 100: peak typical-patient free concentration (~35-40 nmol/L)
## lands peak I_O around ~3500-4000, inside the moderate-to-visible band
## bracketed by Stage 2's own dose-response scan (I_O = 2000 moderate,
## 16000 strong). Illustrative only -- not fit to any bone endpoint.
k_scale <- 100

## ---- Bone parameters (identical to stage1_baseline.R / stage2_denosumab.R,
## except I_O is no longer a fixed parameter -- it is computed above) ------
bone_pars <- c(
  C_s = 0.005, D_A = 0.7, d_B = 0.7, D_C = 0.0021, D_R = 0.0007, f0 = 0.05,
  I_L = 0, I_P = 0, K = 10,
  k1 = 0.01, k2 = 10, k3 = 0.00058, k4 = 0.017, k5 = 0.02, k6 = 3,
  k_B = 0.189, K_L_P = 3000000, kO = 0.35, K_O_P = 200000, k_P = 86,
  r_L = 1000, S_P = 250
)

all_pars <- c(pk_pars, k_scale = k_scale, bone_pars)

## ---- Dosing: 60 mg SC denosumab at months 0, 6, 12 -----------------------
MW_kDa    <- 147
dose_mg   <- 60
dose_nmol <- (dose_mg * 1e6) / (MW_kDa * 1e3)   # ~408.2 nmol per dose

month_d      <- 30.44                # days/month, matches Project 3's convention
dose_times_d <- c(0, 6, 12) * month_d

## ---- Initial conditions: undosed PK state + baseline bone steady state ---
inits <- c(
  depot = 0, Ctot = 0, Cp = 0, Rtot = pk_R0,
  R = 0.0007734, B = 0.0007282, C = 0.0009127
)

## ---- Simulate: 3 doses + ~7 months of post-last-dose follow-up ----------
times <- seq(0, 700, length.out = 7000)

ev <- eventTable()
for (dt in dose_times_d) ev$add.dosing(dose = dose_nmol, start.time = dt, dosing.to = "depot")
ev$add.sampling(times)

sol <- as.data.frame(rxSolve(mod3, all_pars, ev, inits = inits))

## ---- Plot: concentration and bone response on aligned time axes ---------
long <- do.call(rbind, lapply(c("C_deno", "R", "B", "C"), function(v) {
  data.frame(time = sol$time, variable = v, value = sol[[v]])
}))
long$variable <- factor(long$variable, levels = c("C_deno", "R", "B", "C"))

var_labels <- c(
  C_deno = "Free denosumab (nmol/L)",
  R      = "R: Responding osteoblasts",
  B      = "B: Active osteoblasts",
  C      = "C: Active osteoclasts"
)

p <- ggplot(long, aes(time, value)) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_vline(xintercept = dose_times_d, linetype = "dashed", color = "grey40") +
  facet_wrap(~variable, scales = "free_y", ncol = 1, labeller = as_labeller(var_labels)) +
  labs(
    title = "Stage 3: denosumab PK driving I_O in the Lemaire bone model",
    subtitle = "Dashed lines = 60 mg SC doses (months 0, 6, 12).  I_O(t) = k_scale * C_deno(t),  k_scale = 100",
    x = "Time (days)",
    y = NULL
  ) +
  theme_minimal()

print(p)
ggsave("stage3_pk_coupling.png", p, width = 8, height = 10, dpi = 150)

## ---- Report: peak I_O per dose -------------------------------------------
R0_bone <- inits["R"]; B0_bone <- inits["B"]; C0_bone <- inits["C"]
peak_window_days <- 30

cat("\n=== Peak I_O / free denosumab within", peak_window_days, "days after each dose ===\n")
for (dt in dose_times_d) {
  win <- sol[sol$time >= dt & sol$time <= dt + peak_window_days, ]
  i <- which.max(win$I_O)
  cat(sprintf(
    "Dose at day %6.1f: peak I_O = %8.1f  (C_deno = %6.2f nmol/L)  at day %6.1f (+%.1f d post-dose)\n",
    dt, win$I_O[i], win$C_deno[i], win$time[i], win$time[i] - dt
  ))
}

## ---- Report: trough I_O and bone recovery just before each re-dose -------
cat("\n=== Trough I_O and bone-state recovery just before each subsequent dose ===\n")
for (dt in dose_times_d[-1]) {
  idx <- which.min(abs(sol$time - (dt - 0.1)))
  cat(sprintf(
    "Just before dose at day %6.1f (t=%.1f d): I_O=%7.2f  C_deno=%.4f nmol/L | R=%.7f (%.1f%% of baseline)  B=%.7f (%.1f%% of baseline)  C=%.7f (%.1f%% of baseline)\n",
    dt, sol$time[idx], sol$I_O[idx], sol$C_deno[idx],
    sol$R[idx], 100 * sol$R[idx] / R0_bone,
    sol$B[idx], 100 * sol$B[idx] / B0_bone,
    sol$C[idx], 100 * sol$C[idx] / C0_bone
  ))
}

## ---- Report: full washout after the last dose ----------------------------
idx_end <- nrow(sol)
cat(sprintf(
  "\nEnd of simulation (day %.1f, %.1f days after last dose): I_O=%.3f  C_deno=%.4f nmol/L | R=%.7f (%.1f%%)  B=%.7f (%.1f%%)  C=%.7f (%.1f%%)\n",
  sol$time[idx_end], sol$time[idx_end] - max(dose_times_d), sol$I_O[idx_end], sol$C_deno[idx_end],
  sol$R[idx_end], 100 * sol$R[idx_end] / R0_bone,
  sol$B[idx_end], 100 * sol$B[idx_end] / B0_bone,
  sol$C[idx_end], 100 * sol$C[idx_end] / C0_bone
))
