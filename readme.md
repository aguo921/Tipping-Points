# Tipping Points

This project uses numerical simulations to investigate how systems change in response to a slow change in a parameter; in particular, critical transitions, or tipping points, which result in large qualitative changes in the system.

We look at several indicators that could act as early warning signals to upcoming transitions, 
for both time-series models and spatial models.

## Authors

This project was undertaken by Angela Guo under the supervision of Associate Professor Graham Donovan in the 2022-2023 Summer Research Programme for the Department of Mathematics in the University of Auckland.

## Requirements

This project is in Python 3.

Install the required libraries using
```bash
pip install -r "requirements.txt"
```

Optional: Use a virtual environment.

## Directories

### Time-Series Simulations

The `Time-Series-Simulations` directory includes simulations of the stochastic differential equation $dx=f(x,y)dt+\sigma dW$, where $x$ is the dependent variable and $y$ is a slow-changing bifurcation parameter.

### Spatial Simulations

The `Spatial-Simulations` directory includes simulations on 2D spatial models.
