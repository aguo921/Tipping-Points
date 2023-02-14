# Tipping Points

This project uses numerical simulations to investigate how systems change in response to a slow change in a parameter; in particular, critical transitions, or tipping points, which result in large qualitative changes in the system.

We look at several indicators that could act as early warning signals to upcoming transitions, 
for both mean-field models and spatial models.

## Requirements

This project is in Python 3.

Install the required libraries using
```bash
pip install -r "requirements.txt"
```

Optional: Use a virtual environment.

## Directories

### Mean-field Simulations

The `Mean-Field-Simulations` directory includes simulations of the stochastic differential equation $dx=f(x,y)dt+\sigma dW$, where $x$ is the dependent variable and $y$ is a slow-changing bifurcation parameter.

### Spatial Simulations

The `Spatial-Simulations` directory includes simulations on 2D spatial models.
