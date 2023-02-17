# Spatial Simulations

Four spatial models are simulated over a slow change in a bifurcation parameter and spatial indicators are calculated.

## Spatial Models

1. Local Positive Feedback Model

    This model is a discretised stochastic partial differential equation system describing the dynamics of vegetation biomass, $B$, and water level, $w$:

    $\frac{dw_{i,j}}{dt}=R-w_{i,j}-\lambda w_{i,j}B_{i,j}+D(w_{i+1,j}+w_{i-1,j}+w_{i,j+1}+w_{i,j-1}-4w_{i,j})+\sigma_w dW_{i,j}$

    $\frac{dB_{i,j}}{dt}=\rho B_{i,j}\left(w_{i,j}-\frac{B_{i,j}}{B_c}\right)-\mu\frac{B_{i,j}}{B_{i,j}+B_0}+D(B_{i+1,j}+B_{i-1,j}+B_{i,j+1}+B_{i,j-1}-4B_{i,j})+\sigma_B dW_{i,j}$

2. Local Facilitation Model

    This model is a stochastic cellular automaton describing the states of cells: vegetated (+), empty (0) or degraded (-). The grid of cells evolves at each discrete time step through the transition probabilities:

    $w_{[0,+]}=[\delta\rho_+ +(1-\delta)q_{+|0}](b-c\rho_+)$

    $w_{[+,0]}=m$

    $w_{[0,-]}=d$

    $w_{[-,0]}=r+fq_{+|-}$

3. Scale-dependent Feedback Model

    This model is a stochastic partial differential equation system describing the dynamics of plant density, $P$, soil water, $W$ and surface water, $O$:

    $\frac{\partial O}{\partial t}=R(t)-\alpha O\frac{P+W_0 k_2}{P+k_2}+D_O\nabla^2O+\sigma dW$

    $\frac{\partial W}{\partial t}=\alpha O\frac{P+W_0 k_2}{P+k_2}-g_{max}\frac{W}{W+k_1}P-r_W W+D_W\nabla^2W+\sigma dW$

    $\frac{\partial P}{\partial t}=\left(cg_{max}\frac{W}{W+k_1}-d\right)P+D_P\nabla^2P+\sigma dW$

4. Turing Model

    This model is a stochastic partial differential equation system describing the dynamics of an activator, $u$, and an inhibitor, $v$:

    $\frac{\partial u}{\partial t}=u(avu-e)+\nabla^2u+\sigma dW$

    $\frac{\partial v}{\partial t}=v(b-cu^2v)+d\nabla^2v+\sigma dW$

## Spatial Indicators

As a system approaches a critical transition, early warning signals in spatial indicators can arise due to critical slowing down, which occurs as recovery from small perturbations become increasingly slower.

Suppose the state variable $u$ can be represented in a lattice of side length $N$, where $u_{i,j}$ is the value at the location $(i,j)$ and $\bar{u}$ is the spatial mean.

### Spatial Variance

As a system approaches a critical transition, fluctuations around the equilibrium may become stronger due to slower decay.

$$\sigma^2=\frac{1}{N^2}\sum_{i=1}^N\sum_{j=1}^N(u_{i,j}-\bar{u})^2$$

### Spatial Skewness

As a system approaches a critical transition, fluctuations around the equilibrium may become increasingly asymmetric, as the system recovers slower in the direction of the alternative state than the other. Asymmetry can also arise from local flickering events where local units jump between the current and alternative state.

$$\gamma=\frac{1}{N^2}\sum_{i=1}^N\sum_{j=1}^N\frac{(u_{i,j}-\bar{u})^3}{\sigma^3}$$

### Spatial Correlation

As a system approaches a critical transition, local reactions become weaker and diffusion dominates, making neighbouring units more like each other and thus increasingly correlated. Spatial correlation can be measured by nearest-neighbour Moran's I.

$$C=\frac{\sum_{i=1}^N\sum_{j=1}^N (u_{i,j}-\bar{u})(u_{i+1,j}+u_{i-1,j}+u_{i,j+1}+u_{i,j-1}-4\bar{u})}{4\sum_{i=1}^N\sum_{j=1}^N(u_{i,j}-\bar{u})^2}$$

### Power Spectrum

As a system approaches a critical transition, increased memory can cause spectral reddening, where spatial variance becomes increasingly concentrated at lower wavenumbers.

The Discrete Fourier Transform (DFT) decomposes the spatial variation into a summation of sine and cosine waves. The DFT of a spatial state variable is defined for each pair of $x$ and $y$ wavenumbers and is in general a complex number. Hence, we plot the power spectrum, the modulus of the complex matrix, typically up to wavenumbers $p=\frac{N}{2}$ and $q=\frac{N}{2}$ and scaled by the spatial variance $\sigma^2$.

The radial spectrum ($r$-spectrum) can be obtained by the summing the power spectrum at constant distances from the origin. The wavenumber at which the peak occurs is the characteristic spatial frequency of the pattern.

## Function Files

### `spatial_models.py`

This file defines classes for each of the spatial models. Attributes include model parameters and simulation properties. Methods include explicit update equations and equilibria.

This file also includes functions to calculate the Laplace operator approximation, simulate a model for a number of time steps, change a parameter while simulating a model, numerically find the bifurcation point via linearisation and save/load snapshot data.

### `spatial_indicators.py`

This file defines functions to calculate spatial indicators from spatial snapshots, including spatial mean, spatial variance, spatial skewness, spatial correlation and power/radial spectra.

This file also defines functions to bootstrap or coarse grain the snapshot data, which may assist in constructing null models.

Note: Power/radial spectra functions have not been verified.

### `spatial_plotting.py`

This file defines functions to visualise spatial snapshots, plot spatial indicators, plot equilibria and plot maximum eigenvalues.

## Notebooks


### `spatial_simulation.ipynb`

This notebook simulates the spatial models while slowly changing a bifurcation parameter. Spatial snapshots are plotted and animated. Spatial indicators are calculated to see how they change near the bifurcation point. Model information, implementation and observations are documented.

### `local_positive_model.ipynb`

This notebook simulates the local positive feedback model at different levels of noise (both additive and multiplicative) to investigate their effect on the spatial indicators. Simulations are started from different initial conditions to investigate the effect of critical slowing down on recovery to equilibrium. Simulation results are compared to the homogeneous equilibria.

### `local_facilitation_model.ipynb`

This notebook simulates the local facilitation model for different durations at different levels of $b$ to investigate how long the system takes to go to equilibrium from a fully vegetated state as the system approaches the critical transition. Simulation results are compared to the equilibrium from both the mean-field approximation and pair approximation.

### `scale_dependent_model.ipynb`

This notebook simulates the scale-dependent feedback model at different levels of noise (both additive and multiplicative) to investigate the role of noise inducing/breaking patterns. Simulation results are compared to the homogeneous equilibria.

### `turing_model.ipynb`

This notebook simulates the Turing model at different levels of noise (both additive and multiplicative) to investigate the role of noise inducing/breaking patterns. Simulation results are compared to the homogeneous equilibria.

## Simulation Results

The `Results` directory contains the snapshot data simulated in the notebooks. At the end of each notebook, the simulated snapshot data can be saved into the directory. At the start of each notebook, the saved snapshot data can be loaded into the notebook.