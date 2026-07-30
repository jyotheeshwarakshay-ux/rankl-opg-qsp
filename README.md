# rankl-opg-qsp

**A virtual-population QSP model of denosumab's effect on bone remodeling, built on the canonical Lemaire (2004) RANK–RANKL–OPG model.**

This repository couples a published bone cell-population model to a denosumab
pharmacokinetic model and uses the combined system to characterize how denosumab
response varies across a virtual patient population. It is the third component of
a connected pharmacometrics programme centered on denosumab (see
[Related work](#related-work)).

---

## Key finding

Under equal 30% inter-individual variability applied to both pharmacokinetic
(PK) and bone-physiology parameters, variability in denosumab response is
driven **~92% by bone physiology and ~8% by drug exposure** — bone-parameter
differences dominate by roughly 3.5×.

The decomposition is verified rather than asserted: an additive-variance check
(Var[bone-only] + Var[PK-only] ≈ Var[full], ratio 1.04) confirms the two sources
act independently, so the split is not an artifact of interaction effects.

**Interpretation.** Within this model, response variability is dominated by
bone-physiology parameters rather than PK exposure — response heterogeneity
originates upstream of pharmacokinetics in the model structure. (This is a
statement about the model's behavior under the stated variability, not a
clinical claim about real patients.)

*(This finding is conditional on the chosen variability structure — see
[Limitations](#limitations).)*

---

## What this is

A four-stage build on the Lemaire et al. (2004) bone-remodeling model
(J Theor Biol 229(3):293–309; BioModels `BIOMD0000000278`), a canonical
three-state (responding osteoblasts *R*, active osteoblasts *B*, active
osteoclasts *C*) cell-population model with algebraic RANK–RANKL–OPG–PTH control.

| Stage | Script | What it establishes |
|-------|--------|---------------------|
| 1 | `stage1_baseline.R` | Reproduces the published untreated steady state; verified both algebraically (derivatives ≈ 0 at published initial conditions) and by simulation. |
| 1b | `stage1_stability_check.R` | Confirms the baseline is a *stable attractor*: +20% perturbations to R/B/C all return to the same fixed point (~14–68 days). |
| 2 | `stage2_denosumab.R` | Denosumab introduced as an OPG-like RANKL blocker via the model's `I_O` input. Graded, stable suppression of osteoclasts; independently reproduces Lemaire's conclusion that anti-resorptive monotherapy suppresses both sides of the remodeling couple. |
| 3 | `stage3_pk_coupling.R` | Couples a denosumab PK model (from the related TMDD study) so that `I_O` is driven by time-varying free drug concentration. Produces dose-driven suppression cycles under a 60 mg q6-month schedule. |
| 4 | `stage4_virtual_population.R`, `stage4_dosed_response.R` | Generates a plausibility-filtered virtual population and decomposes response variability into PK vs bone contributions (the [key finding](#key-finding)). |

---

## The denosumab thread

Denosumab is a monoclonal antibody against RANKL. This project deliberately
reuses the denosumab pharmacokinetics characterized in a companion repository so
that drug disposition and downstream bone response form a single mechanistic
chain:

**denosumab PK (TMDD)  →  free-drug concentration  →  RANKL blockade (`I_O`)  →  bone cell dynamics  →  population response**

This is the point of the programme: the repositories are meant to read as one
connected line of work rather than independent exercises.

---

## Diagnostic highlights

Two verification steps did real work and are documented in the scripts:

- **The plausibility filter caught a genuine bug.** An early population run
  silently initialized every patient from R = B = C = 0 instead of the true
  baseline steady state (an rxode2 name-matching quirk: `c(R = ...)` produced the
  name `"R.R"`, which the solver could not match). The mandatory pre-dosing
  plausibility check rejected 100/100 patients for the same reason, surfacing the
  bug before any downstream number was trusted.

- **"Zero exclusions" was interrogated, not accepted.** At 30% CV the filter
  excluded no patients. Rather than treat this as success, the parameter
  boundaries were probed directly: k_B fails only near ~1% of its typical value
  and D_R only near ~1e-6× typical — i.e. the failure boundaries sit orders of
  magnitude outside any plausible patient range. The zero-exclusion rate reflects
  genuine structural robustness of these parameters, not a lenient filter.

---

## Limitations

Stated plainly, because knowing where a model stops being trustworthy is part of
using it.

1. **Memoryless drug–effect link.** `I_O` is set directly proportional to free
   denosumab concentration (`I_O = k_scale · C_free`), with no off-rate dynamics.
   A consequence is that osteoclast populations recover to ~100% of baseline
   between doses. This is a *property of the coupling choice*, **not** a
   reproduction of the clinical rebound phenomenon — real denosumab produces
   sustained suppression across the dosing interval and a post-discontinuation
   overshoot above baseline, neither of which this model captures. Representing
   those would require hysteresis (slow off-rate) in the RANKL blockade.

2. **Illustrative, not calibrated, scaling.** `k_scale = 100` was chosen so that
   peak `I_O` lands in the moderate-effect band established in Stage 2. It is a
   transparent free parameter, not fit to clinical dose–response data.

3. **Illustrative virtual population.** Parameters are sampled at 30% CV around
   typical values; the resulting population spans a wide physiological range
   (e.g. ~8-fold spread in baseline active-osteoblast equilibrium) and is *not*
   calibrated to observed patient variability data.

4. **The 92/8 split is conditional.** It holds *under equal 30% CV on both PK and
   bone parameters*. A different variability structure (e.g. wider PK, narrower
   bone) would shift the split. The finding is a statement about this variance
   structure, not a claimed universal property of denosumab.

---

## Reproducing

Requires R with `rxode2` and `ggplot2`. See `SETUP.md` for environment notes
(the model was developed on Apple Silicon macOS; `rxode2` requires `cmake` and an
accepted Xcode license to compile).

```r
# from the repository root
source("stage1_baseline.R")
source("stage1_stability_check.R")
source("stage2_denosumab.R")
source("stage3_pk_coupling.R")
source("stage4_virtual_population.R")
source("stage4_dosed_response.R")
```

Each script is self-contained and writes its own plot.

---

## Sources

- **Bone model:** Lemaire V, Tobin DMS, Greller LD, Cho CR, Suva LJ. *Modeling the
  interactions between osteoblast and osteoclast activities in bone remodeling.*
  J Theor Biol. 2004;229(3):293–309. BioModels `BIOMD0000000278`.
  `BIOMD0000000278.ode` in this repository is the human-readable model export,
  retained for provenance.
- **Denosumab PK parameters:** Choi S, Park S, Jung J, Baek S, Lim H-S. *Population
  pharmacokinetics/pharmacodynamics analysis confirming biosimilarity of SB16 to
  reference denosumab.* Front Pharmacol. 2025;16:1631034.
  doi:10.3389/fphar.2025.1631034 (CC BY). The denosumab PK model used here is the
  two-compartment TMDD model with QSS approximation and first-order absorption
  characterized in that biosimilarity analysis. Notably, the source authors report
  that drug-concentration data alone were insufficient to identify the target
  turnover parameters (kdeg, kint) — the same practical-identifiability limitation
  independently examined in the companion `denosumab-tmdd-qss` study — which is why
  the coupled PK parameters here are treated as illustrative rather than definitive.
  The PK block in `stage3_pk_coupling.R` is annotated with this provenance.

## Related work

- **denosumab-tmdd-qss** — the denosumab TMDD / QSS parameter-recovery and
  identifiability study whose PK model feeds Stage 3 here.
- **tmdd-approximations** — TMDD model with QSS and Michaelis–Menten approximations.
- **Theophylline-popPK** — population-PK reproduction.
