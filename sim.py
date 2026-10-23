"""Pilot: can multiple-systems estimation recover buildings that every footprint list misses?

Villages of buildings with a known count are captured by four dependent lists
(Microsoft, Google, OpenStreetMap, a tax roll). Five estimators are scored on
relative bias, relative RMSE and 95% interval coverage over repeated draws.

    python sim.py            # 200 replicates per scenario, writes results.json + results.md
    python sim.py --reps 20  # quick check
"""
import argparse
import itertools
import json

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit

LISTS = ["microsoft", "google", "osm", "roll"]
Z95 = 1.959964

# Each scenario changes one thing against S1, the realistic baseline.
#   alpha    baseline logit of capture per list
#   imagery  loading of each list on a shared latent visibility (tree cover, roof material);
#            a common loading between lists is what makes them positively dependent
#   size     loading on building size, common to every list (heterogeneity, Mh)
#   osm_vil  sd of a village effect on OSM (mapping effort varies by village)
#   split    probability that a capture fails to match and appears as a separate building
SCENARIOS = {
    "S0 independent": dict(alpha=[0.0, -0.3, -1.0, -0.5], imagery=[0, 0, 0, 0], size=0, osm_vil=0, split=0),
    "S1 shared imagery": dict(alpha=[0.0, -0.3, -1.0, -0.5], imagery=[1.5, 1.5, 0.5, 0], size=0, osm_vil=0.7, split=0),
    "S2 + size heterogeneity": dict(alpha=[0.0, -0.3, -1.0, -0.5], imagery=[1.5, 1.5, 0.5, 0], size=1.0, osm_vil=0.7, split=0),
    "S3 + weak OSM": dict(alpha=[0.0, -0.3, -3.0, -0.5], imagery=[1.5, 1.5, 0.5, 0], size=0, osm_vil=0.7, split=0),
    "S4 + 5% matching error": dict(alpha=[0.0, -0.3, -1.0, -0.5], imagery=[1.5, 1.5, 0.5, 0], size=0, osm_vil=0.7, split=0.05),
    "S5 no roll (3 lists)": dict(alpha=[0.0, -0.3, -1.0, -0.5], imagery=[1.5, 1.5, 0.5, 0], size=0, osm_vil=0.7, split=0, drop_roll=True),
}


def simulate(rng, sc, villages=60, mean_size=300):
    """Return (histories, village ids, true N). histories is units x lists, 0/1."""
    n_v = rng.poisson(mean_size, villages)
    vil = np.repeat(np.arange(villages), n_v)
    N = vil.size
    visibility = rng.normal(size=N)
    size = rng.normal(size=N)
    osm_effect = rng.normal(0, sc["osm_vil"], villages)[vil] if sc["osm_vil"] else 0
    logit = np.array(sc["alpha"])[None, :] + np.outer(visibility, sc["imagery"]) + sc["size"] * size[:, None]
    logit[:, 2] += osm_effect
    caught = rng.random((N, 4)) < expit(logit)
    if sc.get("drop_roll"):
        caught = caught[:, :3]
    if sc["split"]:
        # A failed match removes the capture from the building and adds a one-list unit.
        broken = caught & (rng.random(caught.shape) < sc["split"])
        caught = caught & ~broken
        rows, cols = np.nonzero(broken)
        ghosts = np.zeros((rows.size, caught.shape[1]), bool)
        ghosts[np.arange(rows.size), cols] = True
        caught = np.vstack([caught, ghosts])
        vil = np.concatenate([vil, vil[rows]])
    seen = caught.any(axis=1)
    return caught[seen].astype(int), vil[seen], N


def history_counts(h):
    k = h.shape[1]
    cells = np.array(list(itertools.product([0, 1], repeat=k)))[1:]  # drop the unseen cell
    idx = h @ (1 << np.arange(k)[::-1])
    counts = np.bincount(idx, minlength=2**k)[1:]
    return cells, counts


# ---------- estimators ----------

def naive(h, vil):
    return h.shape[0], None


def chao(h, vil):
    """Chao (1987) lower bound from the frequency of capture counts, with the (k-1)/k factor for k lists."""
    k = h.shape[1]
    f = np.bincount(h.sum(axis=1), minlength=3)
    return h.shape[0] + ((k - 1) / k * f[1] ** 2 / (2 * f[2]) if f[2] else np.nan), None


def poisson_glm(X, y, iters=50):
    """Newton-Raphson Poisson fit. Returns None when a model diverges (a cell driven to zero)."""
    beta = np.linalg.lstsq(X, np.log(y + 0.5), rcond=None)[0]
    with np.errstate(all="ignore"):
        for _ in range(iters):
            mu = np.exp(X @ beta)
            step = np.linalg.solve(X.T @ (mu[:, None] * X), X.T @ (y - mu))
            beta += step
            if not np.isfinite(beta).all() or np.abs(beta).max() > 50:
                return None
            if np.abs(step).max() < 1e-9:
                break
        mu = np.exp(X @ beta)
        cov = np.linalg.inv(X.T @ (mu[:, None] * X))
        ll = np.sum(y * np.log(mu) - mu)
    return (beta, cov, ll) if np.isfinite(cov[0, 0]) else None


def loglinear(h, vil):
    """Fienberg (1972): main effects plus a BIC-chosen set of pairwise interactions.
    Higher-order interactions are fixed at zero, which is what identifies the unseen cell."""
    cells, y = history_counts(h)
    k = cells.shape[1]
    pairs = list(itertools.combinations(range(k), 2))
    best = None
    for r in range(len(pairs) + 1):
        for chosen in itertools.combinations(pairs, r):
            X = np.column_stack([np.ones(len(y)), cells] + [cells[:, a] * cells[:, b] for a, b in chosen])
            if X.shape[1] >= len(y):
                continue
            fit = poisson_glm(X, y)
            if fit is None:
                continue
            beta, cov, ll = fit
            bic = -2 * ll + X.shape[1] * np.log(y.sum())
            if best is None or bic < best[0]:
                best = (bic, beta[0], cov[0, 0])
    _, b0, v0 = best
    n = y.sum()
    m0 = np.exp(b0)
    half = Z95 * np.sqrt(v0)
    return n + m0, (n + m0 * np.exp(-half), n + m0 * np.exp(half))


def _lc_probs(theta, cells, k, classes):
    pi = np.exp(np.r_[0, theta[: classes - 1]])
    pi /= pi.sum()
    p = expit(theta[classes - 1:].reshape(classes, k))
    full = np.array(list(itertools.product([0, 1], repeat=k)))
    lik = np.prod(np.where(full[None], p[:, None, :], 1 - p[:, None, :]), axis=2)  # classes x cells
    return pi @ lik  # probability of each history, unseen first


def _lc_fit(y, cells, k, classes, starts, rng):
    def nll(theta):
        P = _lc_probs(theta, cells, k, classes)
        return -np.sum(y * np.log(P[1:] / (1 - P[0]) + 1e-300))
    fits = [minimize(nll, rng.normal(0, 1, classes - 1 + classes * k), method="BFGS") for _ in range(starts)]
    return min(fits, key=lambda f: f.fun), nll


def _hessian(f, x, eps=1e-4):
    d = x.size
    H = np.empty((d, d))
    E = np.eye(d) * eps
    for i in range(d):
        for j in range(i, d):
            H[i, j] = H[j, i] = (f(x + E[i] + E[j]) - f(x + E[i] - E[j]) - f(x - E[i] + E[j]) + f(x - E[i] - E[j])) / (4 * eps**2)
    return H


def latent_class(h, vil, max_classes=2, starts=6, seed=0):
    """Latent-class MSE on the conditional likelihood, the frequentist analogue of LCMCR.
    One or two classes, chosen by BIC. Two classes need four or more lists: with three,
    they have more parameters than there are observed cells."""
    cells, y = history_counts(h)
    k = cells.shape[1]
    n = y.sum()
    rng = np.random.default_rng(seed)
    best = None
    for classes in range(1, (max_classes if k >= 4 else 1) + 1):
        fit, nll = _lc_fit(y, cells, k, classes, starts, rng)
        bic = 2 * fit.fun + fit.x.size * np.log(n)
        if best is None or bic < best[0]:
            best = (bic, fit, nll, classes)
    _, fit, nll, classes = best
    p0 = lambda t: _lc_probs(t, cells, k, classes)[0]
    P0 = p0(fit.x)
    N_hat = n / (1 - P0)
    m0 = N_hat - n
    # Delta method on the observed information, plus the binomial part of n.
    try:
        cov = np.linalg.inv(_hessian(nll, fit.x))
        eps = 1e-6
        grad = np.array([(p0(fit.x + eps * e) - p0(fit.x - eps * e)) / (2 * eps) for e in np.eye(fit.x.size)])
        var_N = n * P0 / (1 - P0) ** 2 + (n / (1 - P0) ** 2) ** 2 * max(float(grad @ cov @ grad), 0)
    except np.linalg.LinAlgError:
        return N_hat, None
    if m0 <= 0 or not np.isfinite(var_N):
        return N_hat, None
    c = np.exp(Z95 * np.sqrt(np.log(1 + var_N / m0**2)))
    return N_hat, (n + m0 / c, n + m0 * c)


def triple_collocation(h, vil):
    """Per-village counts from the three open lists, merged by triple-collocation weights
    (covariance method, scaled to Microsoft). It weights what the lists hold; it has no
    term for what none of them holds, so its total sits in Microsoft's units."""
    V = vil.max() + 1
    x = np.array([np.bincount(vil[h[:, j] == 1], minlength=V) for j in range(3)], float)
    Q = np.cov(x)
    beta = np.array([1.0, Q[1, 2] / Q[0, 2], Q[1, 2] / Q[0, 1]])
    scaled = (x - x.mean(axis=1, keepdims=True)) / beta[:, None] + x[0].mean()
    err = np.array([Q[0, 0] - Q[0, 1] * Q[0, 2] / Q[1, 2],
                    Q[1, 1] / beta[1] ** 2 - Q[0, 1] * Q[1, 2] / Q[0, 2] / beta[1] ** 2,
                    Q[2, 2] / beta[2] ** 2 - Q[0, 2] * Q[1, 2] / Q[0, 1] / beta[2] ** 2])
    w = 1 / np.clip(err, 1e-9, None)
    return float((w @ scaled / w.sum()).sum()), None


ESTIMATORS = {"naive union": naive, "Chao lower bound": chao, "triple collocation": triple_collocation,
              "log-linear (BIC)": loglinear, "latent class (BIC, K<=2)": latent_class}


def run(reps, seed):
    out = {}
    for s_i, (name, sc) in enumerate(SCENARIOS.items()):
        rows = {e: [] for e in ESTIMATORS}
        for r in range(reps):
            rng = np.random.default_rng([seed, s_i, r])
            h, vil, N = simulate(rng, sc)
            for e, f in ESTIMATORS.items():
                est, ci = f(h, vil)
                rows[e].append((N, est, ci, h.shape[0]))
        out[name] = {}
        for e, res in rows.items():
            N = np.array([a[0] for a in res], float)
            est = np.array([a[1] for a in res], float)
            if np.isnan(est).all():
                continue
            rel = (est - N) / N
            cov = [lo <= n <= hi for n, _, ci, _ in res if ci for lo, hi in [ci]]
            out[name][e] = dict(rel_bias=float(np.nanmean(rel)), rel_rmse=float(np.sqrt(np.nanmean(rel**2))),
                                coverage=float(np.mean(cov)) if cov else None)
        seen = np.mean([a[3] / a[0] for a in rows["naive union"]])
        out[name]["_share_seen"] = float(seen)
        print(f"{name}: done", flush=True)
    return out


def to_markdown(out, reps):
    lines = [f"Replicates per scenario: {reps}. Bias and RMSE relative to the true count.", "",
             "| Scenario | Seen by any list | Estimator | Bias | RMSE | 95% coverage |", "|---|---|---|---|---|---|"]
    for name, res in out.items():
        seen = res.pop("_share_seen")
        for i, (e, m) in enumerate(res.items()):
            cov = "" if m["coverage"] is None else f"{m['coverage']:.0%}"
            lines.append(f"| {name if i == 0 else ''} | {f'{seen:.1%}' if i == 0 else ''} | {e} | "
                         f"{m['rel_bias']:+.1%} | {m['rel_rmse']:.1%} | {cov} |")
        res["_share_seen"] = seen
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=200)
    ap.add_argument("--seed", type=int, default=20261023)
    a = ap.parse_args()
    out = run(a.reps, a.seed)
    json.dump(out, open("results.json", "w"), indent=2)
    open("results.md", "w").write(to_markdown(out, a.reps))
    print(open("results.md").read())
