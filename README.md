# Frequency Estimator

## Problem Description[1]

1. (System) There is a coin of which the head probability oscillates: $$p_{k}(t)
   = \alpha+\beta \cos(2 \pi f_{k}t)$$

1. Frequency $f_{k}$ follows a discrete stochastic process. We assume the Wiener
   process for simplicity.

1. (Observable) For each step, we can toss a coin at time $t_{k}$ we want and get a datum $d
   \in \{0, 1\}$, which represents head and tail, respectively.

1. (Target) Assuming that noise in the Wiener process is small enough respect to $t_{k}$,
   how can we estimate $f_{k}$ efficiently?

## Available Methods and Benchmark Results

### Test condition

1. $\Delta f = \Delta W$, $\Delta W \sim N(f_{0}, 2DT)$

1. 

1. Step 1 and 2 are repeated for $10000$ times and $\textrm{MSE}$ are averaged.

|Method|Result|
|------|------|

## Reference

[1] Shulman, M. D., Harvey, S. P., Nichol, J. M., Bartlett, S. D., Doherty, A.
C., Umansky, V., & Yacoby, A. (2014). Suppressing qubit dephasing using
real-time Hamiltonian estimation. Nature Communications 2014 5:1, 5(1), 1–6.
https://doi.org/10.1038/ncomms6156

[2] Ferrie, C., Granade, C. E., & Cory, D. G. (2012). How to best sample a
periodic probability distribution, or on the accuracy of Hamiltonian finding
strategies. Quantum Information Processing 2012 12:1, 12(1), 611–623.
https://doi.org/10.1007/S11128-012-0407-6
