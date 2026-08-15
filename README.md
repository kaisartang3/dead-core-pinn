# Physics-Informed Neural Network for Diffusion-Reaction Problems with Dead-Core Formation

This repository contains the implementation used for the manuscript *Physics-Informed Neural Network for Diffusion-Reaction Problems with Dead-Core Formation in Catalyst Slabs*. It includes the structured free-boundary PINN, exact-solution utilities, MATLAB shooting-method benchmarks, and figure-generation scripts.

## Repository Contents

- `deadcore_pinn.ipynb` — Jupyter notebook containing the full PINN implementation and numerical experiments.
- `dead_zone_exact_original.py` — Python script for the exact dead-core solution.
- `requirements.txt` — Python dependencies required to run the PINN notebook.
- `proceeding_comparison.m` — MATLAB script used to generate the analytical comparison figure.
- `matlab/` — MATLAB source files for the transformed shooting-method benchmark, convergence checks, timing study, and shooting figures.
- `figures/` — exported PINN plots and manuscript figures.

## Problem Description

The code solves a nonlinear diffusion-reaction problem with an unknown dead-core boundary. The PINN simultaneously learns:

- the concentration profile;
- the dead-core interface location.

The implementation uses a structured trial solution so that the interface and outer-boundary conditions are satisfied by construction.

## PINN Method

The PINN implementation includes:

- domain transformation to a fixed interval;
- asymptotically informed structured ansatz;
- trainable free-boundary parameterization;
- physics-informed residual loss;
- biased collocation-point sampling;
- Adam optimization followed by L-BFGS refinement.

The same fully connected network architecture is used for all test cases: four hidden layers with 64 neurons per hidden layer, hyperbolic tangent (`tanh`) activations, and a linear output layer.

## Python Requirements

Install the required Python packages with:

```bash
pip install -r requirements.txt
```

Then open and run:

```text
deadcore_pinn.ipynb
```

The notebook reproduces the PINN solutions, interface errors, profile errors, and the main PINN figures used in the manuscript.

## MATLAB Shooting Benchmark

The `matlab/` directory contains the classical shooting-method implementation and supporting scripts.

The main reproducibility script is:

```text
matlab/run_all_shooting_study.m
```

To reproduce the shooting study:

1. Open MATLAB.
2. Set the current folder to the repository's `matlab` directory.
3. Run:

```matlab
run_all_shooting_study
```

The script performs:

- epsilon-convergence checks;
- transformed shooting for the representative cases;
- interface-error calculations;
- execution-time measurements;
- shooting solution and error plots;
- export of the final numerical results.

The timing study uses 20 repeated runs after one untimed warm-up call.

## Reproducibility

The repository contains the Python PINN implementation and the MATLAB shooting-method codes used for the numerical comparisons in the manuscript. The main PINN and shooting calculations can be reproduced by running `deadcore_pinn.ipynb` and `matlab/run_all_shooting_study.m`, respectively.

## Citation

If you use this code, please cite the manuscript:

*Physics-Informed Neural Network for Diffusion-Reaction Problems with Dead-Core Formation in Catalyst Slabs.*
