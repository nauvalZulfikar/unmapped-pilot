# unmapped-pilot

How many buildings does every footprint source miss? This pilot asks whether
multiple-systems estimation (capture-recapture with several lists) can answer
that in a simulated district where the true count is known. It is the
feasibility check behind the PhD proposal *Counting the Unmapped*.

## Design

60 villages, about 300 buildings each. Four lists try to capture each building:
Microsoft, Google, OpenStreetMap and a property-tax roll. Capture depends on

- a baseline rate per list,
- a shared latent visibility (tree cover, roof material). Microsoft and Google
  load on it heavily, OSM lightly, the roll not at all. This is what makes the
  imagery lists positively dependent,
- building size (heterogeneity common to every list),
- a village effect on OSM (uneven mapping effort),
- matching error. A capture that fails to match its building becomes a
  separate one-list unit, as happens when footprints from two sources are
  joined by geometry.

| Scenario | What changes against S1 |
|---|---|
| S0 | lists independent |
| S1 | shared imagery, OSM village effect (baseline) |
| S2 | + building-size heterogeneity |
| S3 | + OSM much weaker |
| S4 | + 5% of captures fail to match |
| S5 | roll removed, three lists only |

## Estimators

- **naive union**: buildings seen by any list.
- **Chao lower bound**: n + f1²/(2 f2).
- **triple collocation**: per-village counts from the three open lists,
  error variances by the covariance method, inverse-variance merge. It is
  scaled to Microsoft, so it estimates what the lists hold, not what none holds.
  Included because it is the nearest reference-free method used on building data.
- **log-linear**: Poisson model on the 2^k capture-history table, all subsets of
  pairwise interactions, chosen by BIC. The highest-order interaction is zero
  (the identifying assumption). Wald interval on the log scale.
- **latent class**: conditional likelihood with K = 1 or 2 classes chosen by BIC
  (K = 2 only with four lists). Delta-method variance from a finite-difference
  Hessian, Chao-style log interval.

## Results (200 replicates per scenario, seed 20261023)

Bias and RMSE are relative to the true count. Full table in `results.md`.

| Scenario | Seen by any list | Naive | Chao | Log-linear (coverage) | Latent class (coverage) |
|---|---|---|---|---|---|
| S0 independent | 86.9% | -13.1% | +5.1% | +0.0% (85%) | +0.0% (96%) |
| S1 shared imagery | 83.2% | -16.8% | -1.4% | -0.0% (92%) | -0.1% (94%) |
| S2 + size | 79.0% | -21.0% | -7.1% | -0.6% (90%) | -5.1% (1%) |
| S3 + weak OSM | 79.2% | -20.8% | -2.1% | -0.0% (92%) | -0.0% (97%) |
| S4 + 5% matching error | 89.7% | -10.3% | +16.7% | +21.5% (0%) | +20.3% (0%) |
| S5 no roll | 73.1% | -26.9% | -8.3% | -7.9% (0%) | -17.3% (0%) |

Triple collocation sits at -50% in every scenario, because it reports
Microsoft's scale.

What this shows:

1. Dependence between the imagery lists is handled. With four lists, log-linear
   models recover the count within 1% under shared imagery and a weak OSM.
2. **Matching error dominates.** 5% failed matches inflate every estimator by
   about 20% and no interval covers the truth. Carrying matching error through
   the estimate is therefore the core methodological task, not a detail.
3. **The roll matters.** Without an administrative list made independently of
   imagery, both models under-count by 8 to 17%.
4. Latent class with two classes does not absorb continuous size heterogeneity
   (S2). Log-linear interactions did better here.

## Limits

- Wald intervals ignore model selection, so log-linear coverage sits below 95%
  even when the lists are independent (S0, 85%).
- The latent-class Hessian is a finite-difference approximation.
- Simulated data only. The parameters are plausible, not estimated from Bandung.
- R packages (Rcapture, LCMCR) were not used. Both estimators are written in
  Python with numpy and scipy, so results may differ slightly from those packages.

## Run

```
python -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python sim.py --reps 200   # writes results.json and results.md
```

## References

- Chao, A. (1987). Estimating the population size for capture-recapture data with unequal catchability. *Biometrics* 43(4), 783-791. doi:10.2307/2531532
- Fienberg, S. E. (1972). The multiple recapture census for closed populations and incomplete 2^k contingency tables. *Biometrika* 59(3), 591-603. doi:10.2307/2334810
- Manrique-Vallier, D. (2016). Bayesian population size estimation using Dirichlet process mixtures. *Biometrics* 72(4), 1246-1254. doi:10.1111/biom.12502
- Magyari-Sáska, Z. & Haidu, I. (2025). Built-up surface ensemble model for Romania based on OpenStreetMap, Microsoft Building Footprints, and Global Human Settlement Layer data sources using triple collocation analysis. *ISPRS Int. J. Geo-Inf.* 14(11), 420. doi:10.3390/ijgi14110420
