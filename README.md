# Physics-Informed Neural Network for Diffusion-Reaction Problems with Dead-Core Formation

This repository contains an implementation of a Physics-Informed Neural Network (PINN) for solving nonlinear diffusion-reaction problems with dead-core interface identification.

## Contents
- `deadcore_pinn.ipynb` — Jupyter notebook containing the full PINN implementation
- `requirements.txt` — Python dependencies required to run the notebook
- `dead_zone_exact_original.py` — Python script for the exact dead-zone solution
- `proceeding_comparison.m` — MATLAB script used to generate the comparison figure from the conference proceeding
- `figures/` — exported plots and result figures

## Problem Description
This code solves a nonlinear diffusion-reaction problem with an unknown dead-core boundary. The PINN simultaneously learns:
- the solution profile
- the dead-core interface location

## Methods
The implementation includes:
- domain transformation to a fixed interval
- structured PINN ansatz
- trainable free-boundary parameterization
- physics-informed residual loss
- biased collocation point sampling
- Adam and L-BFGS optimization


## Requirements
Install the required packages using:

```bash
pip install -r requirements.txt
```

## Citation

If you use this code, please cite the manuscript:

*Physics-Informed Neural Network for Diffusion-Reaction Problems with Dead-Zone Formation in Catalyst Slabs.*

