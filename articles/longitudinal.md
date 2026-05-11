# Considerations for Longitudinal fMRI Analysis in fmri.pipeline

## Introduction

Longitudinal fMRI data presents unique statistical challenges,
specifically regarding the proper treatment of variance components
across different time points (sessions) versus repeated within-session
observations (runs). This vignette details the current longitudinal
pathways in `fmri.pipeline`, highlights specific statistical and
computational considerations for FSL users, and shows when a hybrid
**FSL L1/L2 + AFNI `3dLMEr` L3** workflow is preferable.

The package now supports two practically important longitudinal routes:

- FSL-native longitudinal analysis, where session-level maps are modeled
  at Level 3 with subject explanatory variables
  (`l3_input_mode = "pooled_sessions_subject_ev"`).
- A hybrid workflow in which FSL produces session-level COPEs at Level 2
  and AFNI `3dLMEr` fits the longitudinal mixed-effects model at Level 3
  (`l3_input_mode = "3dlmer"`).

This distinction matters. AFNI’s current `3dLMEr` documentation
describes it as the more flexible successor to `3dLME`, especially for
complex random-effects structures. In `fmri.pipeline`, however, `3dLMEr`
is used as a **session-level longitudinal modeler** on top of
FSL-derived L2 inputs rather than as a single end-to-end replacement for
all lower-level FSL estimation.

## Statistical Framework: The Hierarchical fMRI Model

Before examining specific implementations, it is useful to formalize the
hierarchical structure underlying multi-run, multi-session fMRI
analyses. A complete longitudinal fMRI dataset features (at least) four
nested levels of variation: time points within runs, runs within
sessions, sessions within subjects, and subjects within the population.
FSL’s architecture, however, provides only three modeling stages,
forcing one of these variance components to be either collapsed or
ignored.

### Level 1: Within-Run Time Series Model

For subject $`i`$, session $`j`$, and run $`r`$, the Level 1 model is a
standard General Linear Model (GLM) applied to the preprocessed BOLD
time series:

``` math
Y_{ijr} = X_{ijr} \beta_{ijr} + \epsilon_{ijr}
```

where:

- $`Y_{ijr}`$ is the $`T \times 1`$ vector of BOLD signal at a given
  voxel ($`T`$ = number of time points),
- $`X_{ijr}`$ is the $`T \times p`$ design matrix of $`p`$ convolved
  regressors,
- $`\beta_{ijr}`$ is the $`p \times 1`$ vector of parameter estimates
  (PEs), and
- $`\epsilon_{ijr} \sim \mathcal{N}(0, \sigma^2_{ijr} V_{ijr})`$ is
  temporally autocorrelated noise, where $`V_{ijr}`$ captures the
  autocorrelation structure (estimated via FILM prewhitening in FSL).

The outputs of Level 1 are the **contrast of parameter estimates**
(COPEs) and their associated **variance of COPEs** (VARCOPEs), for a
given contrast vector $`c`$:

``` math
\hat{\gamma}_{ijr} = c^T \hat{\beta}_{ijr}, 
\qquad
\text{Var}(\hat{\gamma}_{ijr}) = \hat{\sigma}^2_{ijr} \, c^T (X_{ijr}^T \hat{V}_{ijr}^{-1} X_{ijr})^{-1} c
```

The VARCOPE is critical: it carries forward the **measurement
uncertainty** from each run’s time series fit. All higher-level analyses
in FSL condition on these quantities, treating the L1 COPEs as “data”
and the L1 VARCOPEs as known (or estimated) within-run measurement
error.

### Level 2: Combining Runs (Fixed or Mixed Effects)

At Level 2, the L1 COPEs from multiple runs are combined. For a given
subject $`i`$ and session $`j`$ (when using `l2_scope = "id_session"`),
the model is:

``` math
\hat{\gamma}_{ijr} = Z_{ij} \alpha_{ij} + \eta_{ijr} + e_{ijr}
```

where:

- $`\hat{\gamma}_{ijr}`$ is the L1 COPE for run $`r`$ (treated as the
  “observed data” at L2),
- $`Z_{ij}`$ is the L2 design matrix (often simply a column of ones for
  averaging, or a run-specific contrast),
- $`\alpha_{ij}`$ is the L2 parameter of interest (e.g., the
  session-level mean activation),
- $`\eta_{ijr} \sim \mathcal{N}(0, \sigma^2_{\text{runs}|ij})`$ is the
  **between-run** random effect (present only if using mixed effects;
  zero under fixed effects), and
- $`e_{ijr}`$ has **known** variance $`\text{Var}(\hat{\gamma}_{ijr})`$,
  i.e., the L1 VARCOPE propagated upward.

Under FSL’s **fixed effects** (FE) approach—the default at Level 2—the
between-run variance $`\sigma^2_{\text{runs}}`$ is forced to zero, and
the combination reduces to inverse-variance weighted averaging:

``` math
\hat{\alpha}_{ij}^{\text{FE}} 
= \frac{\sum_r w_{ijr} \, \hat{\gamma}_{ijr}}{\sum_r w_{ijr}},
\qquad
w_{ijr} = \frac{1}{\text{Var}(\hat{\gamma}_{ijr})}
```

``` math
\text{Var}(\hat{\alpha}_{ij}^{\text{FE}}) = \frac{1}{\sum_r w_{ijr}}
```

The critical assumption here is that **runs are exchangeable**—each run
provides an equally informative, independent replicate of the same
underlying signal, differing only in measurement noise. The output
VARCOPE from L2 reflects only **propagated L1 uncertainty**, not genuine
between-run variability.

Under FSL’s **mixed effects** (FLAME) approach, the total variance for
each run’s COPE becomes:

``` math
\text{Var}_{\text{total}}(\hat{\gamma}_{ijr}) = \sigma^2_{\text{runs}|ij} + \text{Var}(\hat{\gamma}_{ijr})
```

FLAME estimates $`\sigma^2_{\text{runs}}`$ from the data, and the
resulting L2 VARCOPE is larger (more conservative) because it
incorporates genuine between-run heterogeneity.

### Level 3: Subject-Level or Group-Level Analysis

The L2 outputs feed into Level 3. In a cross-sectional design, L3 is
typically the group level. In a longitudinal design, L3 may serve as
**either** the session-comparison level (within subject) or the
population level, depending on the approach.

Under the standard cross-sectional design:

``` math
\hat{\alpha}_i = W_i \delta + \zeta_i + u_i
```

where:

- $`\hat{\alpha}_i`$ is the L2 COPE for subject $`i`$,
- $`W_i`$ is the group-level design matrix (intercept plus group
  covariates),
- $`\delta`$ is the population-level effect of interest,
- $`\zeta_i \sim \mathcal{N}(0, \sigma^2_{\text{subj}})`$ is the
  **between-subject** random effect, and
- $`u_i`$ has known variance $`\text{Var}(\hat{\alpha}_i)`$, the L2
  VARCOPE.

FLAME estimates $`\sigma^2_{\text{subj}}`$ from the data, and the
group-level inference accounts for both measurement uncertainty
(propagated from L1 through L2) and genuine between-subject variability.

### The Missing Level: Where Does “Session” Belong?

A complete longitudinal model would require **four** distinct variance
components:

| Level | Source | Variance Component |
|----|----|----|
| L1 | Within-run measurement noise | $`\sigma^2_{\epsilon}`$ (from FILM) |
| L2 | Between-run, within-session | $`\sigma^2_{\text{runs}\|\text{session}}`$ |
| L3 | Between-session, within-subject | $`\sigma^2_{\text{sessions}\|\text{subject}}`$ |
| L4 | Between-subject (population) | $`\sigma^2_{\text{subjects}}`$ |

FSL’s three-level architecture cannot directly accommodate all four.
Therefore, at least one variance component must be either:

1.  **Collapsed** (assumed to be zero, via fixed effects), or
2.  **Conflated** with another component (by pooling data from different
    levels).

This is the fundamental tension underlying all longitudinal FSL
approaches introduced below.

------------------------------------------------------------------------

## Choosing a Pathway

The most defensible workflow depends on what scientific question is
primary.

| Use case | Recommended path | Why |
|----|----|----|
| Balanced baseline/follow-up study, categorical session effect, want to stay inside FEAT/FLAME | `l2_scope = "id_session"` + `l3_input_mode = "pooled_sessions_subject_ev"` | Keeps the analysis fully within FSL and properly models subject means at L3 |
| Need separate session-specific group maps rather than an explicit longitudinal slope | `l2_scope = "id_session"` + `l3_input_mode = "per_session"` | Straightforward per-session inference with modular recomputation |
| Unequal spacing between visits, missing follow-up sessions, continuous time covariates, treatment-by-time interactions, or subject-specific random slopes | `l2_scope = "id_session"` + `l3_input_mode = "3dlmer"` | Uses AFNI `3dLMEr` for a more flexible longitudinal mixed-effects model |
| Near-contiguous acquisitions where “sessions” are essentially pseudo-runs and only a grand subject average is needed | `l2_scope = "id"` | Acceptable only when between-session structure is scientifically unimportant |

## The Challenge of Variance in FSL

FSL’s standard higher-level analysis tool, FLAME (FMRIB’s Local Analysis
of Mixed Effects), is designed intrinsically around hierarchical fixed
and mixed-effects variance modeling. For cross-sectional data with
multiple runs per subject, the standard pipeline is: 1. **Level 1
(Run):** General Linear Model (GLM) for time series data. 2. **Level 2
(Subject):** Fixed-effects combination of multiple runs. 3. **Level 3
(Group):** Mixed-effects (FLAME 1/2) combination across subjects.

When an additional level—session or wave—is introduced, researchers must
decide how to fold the extra temporal dimension into FSL’s rigid 3-level
structure. FSL generally struggles with true repeated-measures modeling
(where observations are non-independent due to within-subject
clustering).

To mitigate this, `fmri.pipeline` offers a few workarounds.

------------------------------------------------------------------------

## Existing Approaches in `fmri.pipeline`

### 1. The “Id” Pooling Approach (`l2_scope = "id"`)

In this approach, **all runs across all sessions** for a specific
subject are pooled together into a single, massive Level 2 fixed-effects
model.

#### How it works:

- At Level 1, each run is modeled independently.
- At Level 2, instead of separating data by `session`, the pipeline
  takes every Level 1 cope for a subject and feeds it into a single
  FLAME model (via `fsl_l2_model.R` and `setup_l2_models.R`).
- To test for session effects (e.g., Session 2 \> Session 1), explicit
  contrast weightings must be constructed across all run-level copes for
  the subject.

#### Statistical Implications of Pooling

By pooling all runs across sessions into a single L2 model, the pipeline
implicitly assumes a **single** between-run variance component that
governs all input COPEs:

``` math
\hat{\gamma}_{ir} = Z_i \alpha_i + \eta_{ir}, \qquad \eta_{ir} \sim \mathcal{N}(0, \sigma^2_{\text{pool}})
```

where the index $`r`$ now runs over *all* runs from *all* sessions
(e.g., $`r = 1, \ldots, R \times J`$ for $`R`$ runs per session and
$`J`$ sessions).

Under fixed effects (the typical L2 setting), $`\sigma^2_{\text{pool}}`$
is forced to zero and runs are combined purely by inverse-variance
weighting. Under mixed effects, FLAME would estimate a single
$`\sigma^2_{\text{pool}}`$. In either case, the model treats all runs as
**exchangeable**, irrespective of which session they came from.

In reality, the variance structure is **heterogeneous and
hierarchical**:

``` math
\text{Cov}(\hat{\gamma}_{ir}, \hat{\gamma}_{ir'}) = 
\begin{cases}
\sigma^2_{\text{runs}|\text{session}} + \sigma^2_{\text{sessions}|\text{subject}} & \text{if } r, r' \text{ are in different sessions} \\
\sigma^2_{\text{runs}|\text{session}} & \text{if } r, r' \text{ are in the same session} \\
\sigma^2_{\text{runs}|\text{session}} + \text{Var}(\hat{\gamma}_{ir}) & \text{if } r = r'
\end{cases}
```

The pooled model conflates $`\sigma^2_{\text{runs}|\text{session}}`$ and
$`\sigma^2_{\text{sessions}|\text{subject}}`$ into a single parameter.
When $`\sigma^2_{\text{sessions}} \gg \sigma^2_{\text{runs}}`$ (often
the case in longitudinal studies spanning months or years), the pooled
model will:

1.  **Underestimate standard errors** for session-level contrasts
    (Session 2 − Session 1), because the between-session variance
    inflates the residual variance rather than being modeled as a
    distinct component.
2.  **Overweight** runs from sessions with lower L1 measurement noise,
    implicitly assuming this reflects higher signal-to-noise for the
    session-level effect when it may simply reflect scanner stability
    differences.

#### Pros:

- **Conceptual Simplicity:** All data for a subject is neatly nested
  under one Level 2 `gfeat` directory.
- **Maintains Pipeline Architecture:** Keeps the standard 3-level (Run
  -\> Subject -\> Group) FSL hierarchy intact by collapsing the
  “Session” dimension at Level 2.

#### Cons & Limitations:

- **Violation of Exchangeable Variance:** FSL’s fixed-effects
  combination at L2 assumes the error variance of the input copes is
  conceptually equivalent or exchangeable. Pooling runs from entirely
  different sessions forces FSL to assume that the variability between
  two runs *within* a session is mathematically identical to the
  variance between two runs *across* sessions (e.g., separated by 1
  year). This is rarely true due to state versus trait variances.
- **Efficiency and Computation:** Generating an L2 analysis with a high
  number of inputs (e.g., 4 runs x 3 sessions = 12 inputs) significantly
  bottlenecks FLAME estimation.
- **Fragility:** If a single bad run in Session 1 is detected, the
  entire multi-session L2 model for that subject must be deleted and
  recomputed, costing substantial wall-time on compute clusters.

### 2. The Hierarchical “Session” Approach (`l2_scope = "id_session"`)

In this widely recommended FSL approach, runs are collapsed within
sessions via fixed-effects, followed by a cross-session mixed-effects or
repeated measures setup.

#### How it works:

- **Level 1:** Run modeling.
- **Level 2 (`l2_scope = "id_session"`):** Runs are combined strictly
  within their respective sessions. A subject with 3 sessions will
  generate 3 separate L2 `gfeat` directories.
- **Level 3 (`l3_input_mode = "pooled_sessions_subject_ev"` or
  `"per_session"`):** The L2 session outputs are fed into a Level 3
  analysis. If modeling longitudinal change internally in FSL, you must
  use a subject-mean approach.

#### Variance Decomposition Under This Approach

This approach correctly separates within-session and between-session
variance:

**At L2** (within-session, fixed effects), for subject $`i`$, session
$`j`$:

``` math
\hat{\alpha}_{ij} = \frac{\sum_r w_{ijr} \hat{\gamma}_{ijr}}{\sum_r w_{ijr}}, \qquad \text{Var}(\hat{\alpha}_{ij}) = \frac{1}{\sum_r w_{ijr}}
```

The L2 VARCOPE reflects **only** within-session measurement
precision—runs from different sessions do not contaminate this estimate.

**At L3** (between-session, mixed effects via the subject-mean
approach):

``` math
\hat{\alpha}_{ij} = \mu_i + \tau \, t_{ij} + \xi_{ij} + u_{ij}
```

where:

- $`\mu_i`$ is a subject-specific intercept (modeled via a dedicated EV
  per subject),
- $`\tau`$ is the population-level slope of the session effect (e.g.,
  time, wave),
- $`t_{ij}`$ is the session indicator or time variable for subject
  $`i`$, session $`j`$,
- $`\xi_{ij} \sim \mathcal{N}(0, \sigma^2_{\text{sessions}|\text{subject}})`$
  is the between-session, within-subject random variation, and
- $`u_{ij}`$ has known variance $`\text{Var}(\hat{\alpha}_{ij})`$, the
  L2 VARCOPE.

This formulation cleanly isolates the three estimable variance
components:

| Component | Interpretation | How Estimated |
|----|----|----|
| $`\text{Var}(\hat{\gamma}_{ijr})`$ | L1 measurement noise | FILM + OLS at L1 |
| $`\text{Var}(\hat{\alpha}_{ij})`$ | Propagated L1 noise, averaged over runs | Inverse-variance weighting at L2 |
| $`\sigma^2_{\text{sessions}\|\text{subject}}`$ | True between-session variability | FLAME at L3 |

#### The Subject-Mean EV Design Matrix

A key practical requirement in this approach is that controlling for
within-subject dependence in FSL FLAME requires constructing a
subject-mean explanatory variable (EV). For $`N`$ subjects each
contributing $`J`$ sessions, the L3 design matrix has the form:

``` math
W = \begin{bmatrix} I_N \otimes \mathbf{1}_J & t \end{bmatrix}
```

where $`I_N \otimes \mathbf{1}_J`$ generates a block-diagonal structure
of $`N`$ indicator columns (one per subject), and $`t`$ is the
session/time covariate of interest.

For a study with $`N = 100`$ subjects and $`J = 3`$ sessions, this
produces a $`300 \times 101`$ design matrix. This rapidly becomes
unwieldy and pushes against FEAT’s practical limits.

#### Pros:

- **Statistical Validity:** Properly isolates within-session variance
  from between-session variance.
- **Modularity:** Re-computing one session’s data does not mandate
  re-computing another session.

#### Cons:

- **Subject Means Required:** The FSL GLM requires defining one
  Explanatory Variable (EV) per subject to model out the
  subject-specific intercept (mean). This rapidly balloons the design
  matrix size for large samples.
- **Balanced Designs Needed:** FLAME assumes relatively balanced
  designs. If continuous time covariates (e.g., days since baseline) or
  missing sessions are introduced, FLAME’s subject-mean approach
  struggles to retain degrees of freedom and robust variance estimates.

### 3. The Hybrid AFNI `3dLMEr` Approach (`l2_scope = "id_session"`, `l3_input_mode = "3dlmer"`)

This is the package’s current flexible pathway for longitudinal Level 3
modeling.

#### How it works:

- **Level 1:** Each run is estimated in FSL as usual.
- **Level 2 (`l2_scope = "id_session"`):** Runs are collapsed within
  session, producing one session-level COPE per subject-session
  combination.
- **Level 3 (`l3_input_mode = "3dlmer"`):** `fmri.pipeline` harvests
  those FSL L2 outputs, writes an AFNI `dataTable.txt`, and generates a
  `3dLMEr` script using the requested fixed and random effects.

#### What AFNI is solving here

Let $`\hat{\alpha}_{ij}`$ denote the session-level COPE for subject
$`i`$ at visit $`j`$ after FSL Level 2 averaging across runs. The AFNI
stage then fits a longitudinal mixed model of the form:

``` math
\hat{\alpha}_{ij} = X_{ij}\beta + b_{0i} + b_{1i} t_{ij} + u_{ij}
```

where:

- $`X_{ij}\beta`$ contains the fixed effects (session, time, group,
  covariates, interactions),
- $`b_{0i}`$ is a subject-specific random intercept,
- $`b_{1i}`$ is an optional subject-specific random slope for time or
  session, and
- $`u_{ij}`$ is the residual session-level error.

This is more accurate than pooled FSL L2 for most longitudinal questions
because the repeated session structure is modeled explicitly. At the
same time, it is important to be precise about what is and is not
estimated here: **run-level combination has already happened in FSL**.
In other words, this is a hybrid two-stage workflow, not a single model
that jointly re-estimates run-, session-, and subject-level variance
from raw L1 outputs.

#### Package-specific implementation details

- `3dlmer` is only compatible with `l2_scope = "id_session"`.
- The AFNI backend currently consumes **FSL-produced subject-session
  contrasts** at L3.
- Numeric predictors in the refit L3 model are automatically written
  under AFNI’s `-qVars`.
- Random effects are specified with `random_effects`, then appended to
  the AFNI `-model` formula.
- `fmri.pipeline` uses `Subj` as the subject column in the AFNI data
  table, which matches AFNI’s recommended conventions.

#### Pros:

- **Flexible random-effects structure:** Random intercepts and random
  slopes are easier to express than in FSL’s subject-EV workaround.
- **More natural longitudinal covariates:** Continuous time,
  age-at-scan, treatment duration, and similar predictors fit naturally
  into the model.
- **Missing visits are easier to tolerate:** AFNI mixed models are
  generally much better suited than FSL FEAT to incomplete longitudinal
  schedules.
- **Clearer modeling of treatment-by-time questions:** Interactions such
  as `group * days_since_baseline` are straightforward.

#### Cons:

- **Hybrid, not fully joint:** The model does not explicitly carry L2
  VARCOPEs into AFNI as observation-specific weights.
- **Extra software dependency:** AFNI plus its required R packages must
  be available on the compute environment.
- **Interpretation requires mixed-model discipline:** Centering of
  quantitative variables and selection of random-effects terms matter.

------------------------------------------------------------------------

## Error Propagation Across Levels

A key feature of FSL’s hierarchical framework is that **uncertainty
propagates upward** from level to level. Misspecification at any stage
cascades through all subsequent stages.

### Formal Error Propagation

At each level $`\ell`$, the total variance of a parameter estimate is
the sum of:

1.  **Propagated uncertainty** from level $`\ell - 1`$: the known
    (estimated) measurement variance of the input data.
2.  **New variance** at level $`\ell`$: the between-unit variability
    estimated at the current level.

Concretely, for FLAME at level $`\ell`$ with inputs $`\hat{\theta}_k`$
(the COPEs from level $`\ell - 1`$):

``` math
\hat{\theta}_k = \delta^{(\ell)} + \zeta_k + e_k
```

where $`\zeta_k \sim \mathcal{N}(0, \sigma^2_{(\ell)})`$ and $`e_k`$ has
known variance $`s^2_k`$ (the VARCOPE from level $`\ell-1`$). The
marginal variance of each input is:

``` math
\text{Var}(\hat{\theta}_k) = \sigma^2_{(\ell)} + s^2_k
```

The resulting estimate and its variance are:

``` math
\hat{\delta}^{(\ell)} = \frac{\sum_k \tilde{w}_k \hat{\theta}_k}{\sum_k \tilde{w}_k}, \qquad \tilde{w}_k = \frac{1}{\sigma^2_{(\ell)} + s^2_k}
```

``` math
\text{Var}(\hat{\delta}^{(\ell)}) = \frac{1}{\sum_k \tilde{w}_k}
```

### Consequences of Misspecification

1.  **Underestimated L1 VARCOPEs** (e.g., from inadequate
    autocorrelation modeling): The L1 error $`s^2_k`$ is too small,
    causing FLAME at L2/L3 to **overweight** noisy runs and produce
    **anti-conservative** group-level inference.

2.  **Collapsed variance components** (e.g., the “id” pooling approach):
    When $`\sigma^2_{\text{sessions}}`$ is not explicitly modeled, the
    between-session variance inflates $`\sigma^2_{(\ell)}`$ at the
    pooled level (or is absorbed into the residual). Session-level
    contrasts then have **incorrect standard errors**, since FLAME’s
    single random-effects variance cannot simultaneously capture both
    between-run and between-session variability.

3.  **Compounding bias across levels**: If variance is systematically
    misestimated at L1 (e.g., due to scanner artifact at one session),
    the error propagates multiplicatively:

``` math
\text{Var}(\hat{\delta}_{\text{group}}) = \frac{1}{\sum_i \frac{1}{\hat{\sigma}^2_{\text{subj}} + \text{Var}(\hat{\alpha}_i)}}
```

where $`\text{Var}(\hat{\alpha}_i)`$ itself depends on
$`\text{Var}(\hat{\gamma}_{ijr})`$. A biased L1 VARCOPE shifts the
effective weight for that subject at L3, distorting the group-level
estimate in a manner that is difficult to diagnose post hoc.

### Illustration: Two-Session, Two-Run Design

Consider a simple design with $`J = 2`$ sessions and $`R = 2`$ runs per
session. The complete data vector for subject $`i`$ is
$`\{\hat{\gamma}_{i11}, \hat{\gamma}_{i12}, \hat{\gamma}_{i21}, \hat{\gamma}_{i22}\}`$.

**Under id_session (correct decomposition):**

- L2 produces two session-level COPEs: $`\hat{\alpha}_{i1}`$ and
  $`\hat{\alpha}_{i2}`$, each with VARCOPEs reflecting only
  within-session run variability.
- L3 can estimate $`\sigma^2_{\text{sessions}|i}`$ from the two session
  COPEs, correctly attributing between-session variance.
- A session contrast $`\hat{\alpha}_{i2} - \hat{\alpha}_{i1}`$ has
  variance:

``` math
\text{Var}(\hat{\alpha}_{i2} - \hat{\alpha}_{i1}) = 2 \sigma^2_{\text{sessions}|i} + \text{Var}(\hat{\alpha}_{i1}) + \text{Var}(\hat{\alpha}_{i2})
```

**Under id pooling (conflated variance):**

- L2 pools all four runs. Under fixed effects, the session contrast is:

``` math
\hat{\Delta}_i^{\text{pool}} = \frac{1}{2}\bigl(\hat{\gamma}_{i21} + \hat{\gamma}_{i22}\bigr) - \frac{1}{2}\bigl(\hat{\gamma}_{i11} + \hat{\gamma}_{i12}\bigr)
```

with variance:

``` math
\text{Var}(\hat{\Delta}_i^{\text{pool}}) = \frac{1}{4} \sum_{j,r} c_{jr}^2 \, \text{Var}(\hat{\gamma}_{ijr})
```

where $`c_{jr}`$ are the L2 contrast weights ($`+1/2`$ for Session 2
runs, $`-1/2`$ for Session 1 runs). This expression **ignores**
$`\sigma^2_{\text{sessions}}`$ entirely—it treats the session difference
as solely a function of L1 measurement noise, producing standard errors
that can be dramatically too small.

------------------------------------------------------------------------

## Practical Recipes and Best Practices

Given the rigid limitations of adding a 4th level (Group) natively into
FSL FLAME, and considering FSL’s limited handling of true linear
mixed-effects designs, the most useful question is not whether one
platform is globally “best”, but which route best matches the
longitudinal structure of your study.

### A. The Subject-Mean Approach (`l3_input_mode = "pooled_sessions_subject_ev"`)

This is the statistically preferred FSL-native approach for pooled
session analysis. It separates the session effect from the
subject-specific baseline.

#### How it works:

- **Automatic Injection:** Starting with the current version, if `id` is
  not present in the L3 formula, the pipeline **automatically injects**
  a factor EV for each subject.
- **Formula Transformation:** A formula like `~ session` is internally
  transformed to `~ 0 + id + session`.
- **Statistical Benefit:** This models a separate intercept for every
  subject, allowing FSL to estimate the session effect while controlling
  for subject-level means. It effectively performs the standard
  longitudinal partitioning of variance within the FSL model.

**Why use this over Id Pooling?** It correctly treats sessions as
repeated observations within subjects and models the individual
differences in baseline activation levels.

**Best fit:** Balanced designs with a modest number of sessions,
categorical visit indicators, and a strong preference to keep the entire
analysis inside FEAT/FLAME.

#### Example configuration

``` yaml
l2_models:
  session_summary:
    level: 2
    l2_scope: id_session
    model_formula: ~ 1
    contrasts:
      diagonal: yes

l3_models:
  longitudinal_fsl:
    level: 3
    l3_input_mode: pooled_sessions_subject_ev
    model_formula: ~ session_label
    reference_level:
      session_label: baseline
    contrasts:
      diagonal: yes
      cond_means: session_label
      pairwise_diffs: session_label
      overall_response: yes
    fsl_outlier_deweighting: no
```

**Why PALM?** PALM explicitly supports multi-level exchangeability
blocks (EBs), allowing you to tell the permutation test exactly which
runs belong to which sessions, and which sessions belong to which
subjects.

### B. Route into AFNI’s `3dLMEr`

For more flexible mixed-effects modeling at the session level: 1.
Process Level 1 in FSL. 2. Collapse runs within session at Level 2 using
`l2_scope = "id_session"`. 3. Set the Level 3 model to
`l3_input_mode = "3dlmer"` so that `fmri.pipeline` writes the AFNI
script and `dataTable.txt` automatically.

AFNI’s current documentation recommends `3dLMEr` over `3dLME` when the
random-effects structure becomes more complex. That matches the package
implementation: the supported longitudinal AFNI path in `fmri.pipeline`
is `3dLMEr`.

**Best fit:** Irregular visit timing, continuous time metrics, partially
missing follow-up, or any design where you want subject-specific random
slopes rather than only subject intercept EVs.

#### Example 1: Categorical visit effect with random intercepts

``` yaml
glm_software: [fsl, afni]

l2_models:
  session_summary:
    level: 2
    l2_scope: id_session
    model_formula: ~ 1
    contrasts:
      diagonal: yes

l3_models:
  longitudinal_3dlmer:
    level: 3
    l3_input_mode: 3dlmer
    execution_backend: afni
    producer_backend: fsl
    model_formula: ~ session_label
    random_effects: (1 | Subj)
    reference_level:
      session_label: baseline
    contrasts:
      diagonal: yes
      cond_means: session_label
      pairwise_diffs: session_label
      overall_response: yes
    fsl_outlier_deweighting: no
```

Here, `execution_backend: afni` tells the pipeline that this model is
estimated by AFNI `3dLMEr`, and `producer_backend: fsl` indicates that
the Level 2 inputs it consumes were produced by FSL. These fields can
also be set programmatically on the `gpa` object (see the **Choosing
model backends** section of
[`vignette("design")`](https://uncdependlab.github.io/fmri.pipeline/articles/design.md)).

#### Example 2: Unequally spaced visits with a random slope for time

``` yaml
glm_software: [fsl, afni]

l2_models:
  session_summary:
    level: 2
    l2_scope: id_session
    model_formula: ~ 1
    contrasts:
      diagonal: yes

l3_models:
  longitudinal_time_slope:
    level: 3
    l3_input_mode: 3dlmer
    model_formula: ~ days_since_baseline + group + age_baseline
    random_effects: (1 + days_since_baseline | Subj)
    contrasts:
      diagonal: yes
      overall_response: yes
    fsl_outlier_deweighting: no
```

In this second pattern, `days_since_baseline` and `age_baseline` are
numeric covariates. `fmri.pipeline` will emit them as AFNI `-qVars`
automatically. AFNI also recommends thinking carefully about centering
quantitative variables, so centering time at baseline or at a clinically
meaningful visit is often preferable.

#### Example 3: Treatment-by-time interaction

``` yaml
glm_software: [fsl, afni]

l2_models:
  session_summary:
    level: 2
    l2_scope: id_session
    model_formula: ~ 1
    contrasts:
      diagonal: yes

l3_models:
  treatment_time_interaction:
    level: 3
    l3_input_mode: 3dlmer
    model_formula: ~ days_since_baseline * treatment_group
    random_effects: (1 + days_since_baseline | Subj)
    contrasts:
      diagonal: yes
      overall_response: yes
    fsl_outlier_deweighting: no
```

This is precisely the kind of model that becomes awkward inside FSL’s
subject-EV framework but is natural in `3dLMEr`.

### C. Restrict “Id” Pooling to Specific Use Cases

The `"id"` pooling approach should generally be deprecated as a default
option for long-term longitudinal fMRI studies, as the statistical
tradeoff is severe. It should only be used when: - Sessions are spaced
by merely hours or minutes (essentially making them contiguous
pseudo-runs). - The researcher explicitly wishes to compute a grand
average across all data for a subject with no interest in temporal
dynamics.

## Conclusion

While `fmri.pipeline` can mathematically construct the `"id"` pooling
approach, it sacrifices statistical validity and computational
modularity. In most real longitudinal studies, the package should
prioritize `l2_scope = "id_session"` and then choose between two
downstream paths:

- FSL subject-EV modeling when the design is balanced and relatively
  simple.
- AFNI `3dLMEr` when the longitudinal structure is genuinely
  mixed-effects in nature.

The formal treatment of variance components and error propagation above
makes explicit that these are not mere convenience choices. The main
current lesson is simple: avoid conflating run and session variance, and
use `3dLMEr` when the scientific question demands longitudinal
flexibility that FEAT/FLAME does not provide gracefully.
