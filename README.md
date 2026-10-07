Get Started with cp-VAR
================
cp-VAR is a MATLAB code package for implementing Bayesian estimation of change-point vector autoregressions (VARs) in my [paper](https://github.com/econdojo/papers/blob/main/pdf/BinUtil.pdf) "Appetite for Treasuries, Debt Cycles, and Fiscal Inflation." The model is a VAR with multiple change points and non-recurrent state transitions, estimated by a Metropolis-within-Gibbs sampler. While the application studies U.S. inflation and fiscal stance, the methodology is general and can be adapted to other data sets. For any issue, suggestion or bug report, please send an email to [econdojo [at] gmail.com](mailto:econdojo@gmail.com).

Package Structure
-----------------------------------
The package consists of the following files in the master folder and a collection of helper routines under the "utils" folder.

1. **mcmc.m** is the main script. It sets the number of draws, burn-in, lag order, number of states, and prior/proposal settings, builds the VAR data, initializes the parameters and state path, runs the sampler, and then performs posterior analysis.

2. **evalPost.m** returns the negative log posterior kernel of the state-dependent population means. It is minimized by **csminwel.m** to center the Metropolis-Hastings proposal.

3. **evalLik.m** runs the filtering recursion of [Chib (1998)](https://doi.org/10.1016/S0304-4076(97)00115-2) for the change-point chain and returns the predictive and updated state probabilities and the log likelihood.

4. **ssEqm.m** imposes the steady-state relation nominal rate = inflation + GDP growth + real rate/surplus-debt ratio.

5. **data.txt** holds the quarterly U.S. data, one column per observable: GDP growth, inflation, interest rate, surplus-debt ratio, and credit spread.

The "utils" folder contains general-purpose routines:

* Optimization: **csminwel.m** (Chris Sims' quasi-Newton optimizer) with **csminit.m**, **bfgsi.m**, and **numgrad.m**.
* Distributions: **mvt_pdf.m** and **mvt_rnd.m** (multivariate Student-t), **prior_pdf.m** and **prior_rnd.m** (priors).
* Linear algebra: **cholmod.m** (modified Cholesky), **logdet.m**, and **sum2one.m**.
* Posterior analysis: **PostStat.m**, **ProbBand.m** (HPD bands), **IneffFactor.m**, and **NeweyWest.m**.
* Plotting: **NBERbc.m**, **NBERbcYY.m**, **tight_subplot.m**, and **progressbar.m**.

Replication
-----
To run the sampler, set the MATLAB directory to the master folder and run **mcmc.m**, which adds the "utils" folder to the search path automatically. The default settings are 11,000 draws with 1,000 burn-in, three states (two change points), and Beta priors on the transition probabilities.

Suggested Citation
-----
Tan, Fei, Appetite for Treasuries, Debt Cycles, and Fiscal Inflation, *Macroeconomic Dynamics*, Volume 30, August 2026, Page 1&ndash;17.
