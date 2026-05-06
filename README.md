# Optimal Factorization for Multistage Decimation

This repository contains a Python implementation of the algorithm proposed in the paper:

**"Optimized Multistage Decimation Based on Optimal Factorization of Decimation Ratio"**
*Vukasin Babic and Djordje Babic*, *IEEE Signal Processing Letters*, 2025.
DOI: 10.1109/LSP.2025.3595512

## Description

The algorithm computes the optimal factorization of a decimation ratio $R$ to minimize computational load, supporting both FIR-only and TMFS + FIR multistage architectures.

## Usage

```python
from optimal_factorization import find_optimal_factorization

delta_p = 0.01
delta_s = 0.001
f_p = 0.9
R = 97
F_in = 1

cost, factors = find_optimal_factorization(delta_p, delta_s, f_p, R, F_in)
print("Computation cost:", cost)
print("FIR factors:", factors)
```

## Poster

ICASSP 2026 poster presentation:

- [ICASSP 2026 Poster (PDF)](icassp2026_poster.pdf)
