# Physics-Informed Neural Network for Diffusion-Reaction Problems with Dead-Core Formation

This repository contains an implementation of a Physics-Informed Neural Network (PINN) for solving nonlinear diffusion-reaction problems with dead-core interface identification.

## Contents
- `Deadcore PINN.ipynb` — Jupyter notebook containing the full implementation
- `requirements.txt` — Python dependencies required to run the notebook

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

