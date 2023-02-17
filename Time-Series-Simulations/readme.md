# Time-Series Simulations

Stochastic fast-slow dynamical systems are simulated over a slow change in a bifurcation parameter and temporal indicators are calculated.

The fast-slow dynamical system can be repesented by

$dx=f(x,y)dt+\sigma dW$

$dy=\epsilon dt$

where $x$ is a fast variable and $y$ is a slow (bifurcation) variable. $\sigma$ is the standard deviation of the white noise and $\epsilon$ is the time-scale separation.

## Bifurcations

We simulate the normal forms of three bifurcations.

1. Fold Bifurcation

    The normal form of the fold bifurcation is $f(x,y)=-y-x^2$. Below the bifurcation point ($y<0$), there exists a stable equilibrium ($x=\sqrt{-y}$) and an unstable equilibrium ($x=-\sqrt{-y}$). Above the bifurcation point ($y>0$), there exists no equilibrium.

2. Transcritical Bifurcation

    The normal form of the transcritical bifurcation is $f(x,y)=yx-x^2$. Below the bifurcation point ($y<0$), there exists a stable equlibrium ($x=0$) and an unstable equilibrium ($x=y$). Above the bifurcation point ($y>0$), there exists a stable equilibrium ($x=y$) and an unstable equilibrium ($x=0$). The stabilities of the two equilibria essentially switch at the bifurcation point.

3. Pitchfork Bifurcation

    The normal form of the pitchfork bifurcation is $f(x,y)=yx+x^3$. Below the bifurcation point ($y<0$), there exists two unstable equilibria ($x=\pm\sqrt{-y}$) and one stable equilibrium ($x=0$). Above the bifurcation point ($y>0$), there exists no equilibrium.

## Files

### `simulation_functions.py`

This file defines functions to generate sample paths, calculate the percentage of escaped trajectories, calculate the estimate and analytical variance and calculate the probability density.

### `simulation_plotting.py`

This file defines functions to plot sample paths, percentage of escaped trajectories, distribution of $x$ values at each $y$ value, probability surface plots and variance.

### `simulation.ipynb`

This notebook simulates an ensemble of sample paths for each of the three normal forms of the bifurcations. Sample paths are plotted against the equilibria and on the probability surface plot. Variance and percentage escaped are plotted against the bifurcation parameter $y$.

### `app.py`

This file creates an app interface with user input on the simulation settings. Plots are generated and can be compared between different bifurcations.

Run

```bash
streamlit run Time-Series-Simulations/app.py
```

to open the app.