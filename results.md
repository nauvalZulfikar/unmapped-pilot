Replicates per scenario: 200. Bias and RMSE relative to the true count.

| Scenario | Seen by any list | Estimator | Bias | RMSE | 95% coverage |
|---|---|---|---|---|---|
| S0 independent | 86.9% | naive union | -13.1% | 13.1% |  |
|  |  | Chao lower bound | +0.6% | 0.7% |  |
|  |  | triple collocation | -50.0% | 50.0% |  |
|  |  | log-linear (BIC) | +0.0% | 0.4% | 85% |
|  |  | latent class (BIC, K<=2) | +0.0% | 0.4% | 96% |
| S1 shared imagery | 83.2% | naive union | -16.8% | 16.8% |  |
|  |  | Chao lower bound | -5.2% | 5.3% |  |
|  |  | triple collocation | -50.0% | 50.0% |  |
|  |  | log-linear (BIC) | -0.0% | 0.5% | 92% |
|  |  | latent class (BIC, K<=2) | -0.1% | 0.8% | 94% |
| S2 + size heterogeneity | 79.0% | naive union | -21.0% | 21.0% |  |
|  |  | Chao lower bound | -10.5% | 10.6% |  |
|  |  | triple collocation | -50.0% | 50.0% |  |
|  |  | log-linear (BIC) | -0.6% | 1.2% | 90% |
|  |  | latent class (BIC, K<=2) | -5.1% | 5.2% | 1% |
| S3 + weak OSM | 79.2% | naive union | -20.8% | 20.8% |  |
|  |  | Chao lower bound | -6.8% | 6.8% |  |
|  |  | triple collocation | -50.0% | 50.0% |  |
|  |  | log-linear (BIC) | -0.0% | 0.7% | 92% |
|  |  | latent class (BIC, K<=2) | -0.0% | 1.1% | 97% |
| S4 + 5% matching error | 89.7% | naive union | -10.3% | 10.4% |  |
|  |  | Chao lower bound | +9.9% | 10.0% |  |
|  |  | triple collocation | -50.0% | 50.0% |  |
|  |  | log-linear (BIC) | +21.5% | 21.6% | 0% |
|  |  | latent class (BIC, K<=2) | +20.3% | 20.4% | 0% |
| S5 no roll (3 lists) | 73.1% | naive union | -26.9% | 26.9% |  |
|  |  | Chao lower bound | -14.5% | 14.5% |  |
|  |  | triple collocation | -50.0% | 50.0% |  |
|  |  | log-linear (BIC) | -7.9% | 8.0% | 0% |
|  |  | latent class (BIC, K<=2) | -17.3% | 17.4% | 0% |
