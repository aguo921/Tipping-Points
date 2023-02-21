# Spatial Simulations

Four spatial models are simulated over a slow change in a bifurcation parameter and spatial indicators are calculated.

## Spatial Models

### Local Positive Feedback Model

This model is a discretised stochastic partial differential equation system based on the vegetation model by Shnerb et al. (2003) [[1]](https://doi.org/10.1103/PhysRevLett.90.038101) and Guttal and Jayaprakash (2007) [[2]](https://doi.org/10.1016/j.ecolmodel.2006.10.005):

$\frac{dw_{i,j}}{dt}=R-w_{i,j}-\lambda w_{i,j}B_{i,j}+D(w_{i+1,j}+w_{i-1,j}+w_{i,j+1}+w_{i,j-1}-4w_{i,j})+\sigma_w dW_{i,j}$

$\frac{dB_{i,j}}{dt}=\rho B_{i,j}\left(w_{i,j}-\frac{B_{i,j}}{B_c}\right)-\mu\frac{B_{i,j}}{B_{i,j}+B_0}+D(B_{i+1,j}+B_{i-1,j}+B_{i,j+1}+B_{i,j-1}-4B_{i,j})+\sigma_B dW_{i,j}$

Vegetation biomass $B$ grows logistically based on water availability $w$ and is lost due to grazing. Water level $w$ increases from rainfall and decreases due to evaporation/percolation and consumption from plants.

Biomass and water are exchanged between neighbouring cells through diffusion. Hence each cell has a positive effect on its neighbours and a negative effect on itself. The diffusive effect tends to smooth out differences between cells when there is strong enough connectivity between cells.

### Local Facilitation Model

This model is a stochastic cellular automaton (Kefi et al. 2007) [[3]](https://doi.org/10.1016/j.tpb.2006.09.003) describing the states of cells: vegetated (+), empty (0) or degraded (-). The grid of cells evolves at each discrete time step through the transition probabilities:

1. Colonisation: $w_{[0,+]}=(\delta\rho_+ +(1-\delta)q_{+|0})(b-c\rho_+)$

2. Mortality: $w_{[+,0]}=m$

3. Degradation: $w_{[0,-]}=d$

4. Regeneration: $w_{[-,0]}=r+fq_{+|-}$

Colonisation and regeneration probabilities of a cell are positively affected by the presence of vegetation in neighbouring cells.

### Scale-dependent Feedback Model

This model is a stochastic partial differential equation system describing the dynamics of plant density, $P$, soil water, $W$ and surface water, $O$ (Rietkerk et al. 2002) [[4]](https://doi.org/10.1086/342078):

$\frac{\partial O}{\partial t}=R(t)-\alpha O\frac{P+W_0 k_2}{P+k_2}+D_O\nabla^2O+\sigma dW$

$\frac{\partial W}{\partial t}=\alpha O\frac{P+W_0 k_2}{P+k_2}-g_{max}\frac{W}{W+k_1}P-r_W W+D_W\nabla^2W+\sigma dW$

$\frac{\partial P}{\partial t}=\left(cg_{max}\frac{W}{W+k_1}-d\right)P+D_P\nabla^2P+\sigma dW$

Plants grow depending on soil water availability and lost from natural mortality or grazing. Surface water is supplied by rainfall and lost due to soil infiltration and run off. Soil water due to soil infiltration is consumed by plants or lost by run off.

Plant dispersal via seed or vegetative propagation, lateral surface water flow due to pressure differences and lateral subsurface water flow due to capillary forces can be approximated by diffusion terms.

The infiltration rate of water into soil is higher in vegetated areas. This leads to more surface water running off from areas with bare soil and subsequently being infiltrated into more vegetated soil so that water accumulates under vegetation and depletes further away. This can thus lead to the formation of self-organised patterns.

### Turing Model

This model is a stochastic partial differential equation system describing the dynamics of an activator, $u$, and an inhibitor, $v$ (Borgogno et al. 2009) [[5]](https://doi.org/10.1029/2007RG000256):

$\frac{\partial u}{\partial t}=u(avu-e)+\nabla^2u+\sigma dW$

$\frac{\partial v}{\partial t}=v(b-cu^2v)+d\nabla^2v+\sigma dW$

The growth rate of species $u$ increases with increasing values of $u$ and $v$. The species $u$ has a strong negative influence (inhibition) on the growth rate of $v$.

As the inhibitor diffuses much faster in space than the activitor, this leads to short-range positive feedback and long-range negative feedback, which can result in the formation of self-organised patterns.

**Note**: This model is a system of two partial differential equations exhibiting Turing patterns, serving as a simpler version of the scale-dependent feedback model but without ecological context.

## Simulation

### Parameter Values

The model parameters for the local positive feedback model, local facilitation model and scale-dependent feedback model are taken from Dakos et al. (2011) [[6]](https://doi.org/10.1086/659945). The model parameters for the Turing model are taken from Borgogno et al. (2009) [[5]](https://doi.org/10.1029/2007RG000256).

For the local positive feedback model, the parameter $B_c$ was amended from 1 to 10.

### Change in Environmental Parameter

We simulated each model while slowly changing an environmental parameter. For the local positive feedback model and the scale-dependent feedback model, this is the rainfall value $R$. For the local facilitation model, this is the aridity level $b$, where lower values indicate harsher conditions. For the Turing model, this is the parameter $b$, which controls the growth of the inhibitor $v$.

At each control parameter level, we simulated many time steps to allow the system to settle to equilibrium. We took the initial state at each control parameter level to be the final stationary state of the previous control parameter level.

We started the simulation at the fully vegetated state for the local positive feedback model and the local facilitation model. We started the simulation at the homogeneous equilibrium for the scale-dependent feedback model and the Turing model.

### Effect of Noise

We simulated each model while slowly changing an environmental parameter with different levels of noise on the state variables. We also simulated the effect of additive noise and multiplicative, or state-dependent, noise.

We also simulated the scale-dependent feedback model and Turing model near the Turing instability point while slowly increasing the noise to investigate the effect of noise on pattern formation.

### Preventive Measures

We simulated each model while slowly changing an environmental parameter. When the environmental parameter reaches a certain value near the collapse of vegetation, we begin reducing a second parameter representing grazing/mortality to try to prevent collapse.

We simulate different environmental parameter rates with different grazing/mortality reduction rates to see when the preventive measures result in recovery of the system.

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

## Numerical Methods

### Finite Difference Method

The spatial system is discretised into a lattice with periodic boundaries. The finite difference method is used to approximate the Laplace operator.

The 5-point stencil approximation of the Laplace operator

$\nabla^2 u_{i,j}=\frac{1}{\Delta x^2}(u_{i+1,j}+u_{i-1,j}+u_{i,j+1}+u_{i,j-1}-4u_{i,j})$

takes into account diffusion across nearest (top/bottom and left/right) neighbours and is generally numerically stable for sufficiently smooth fields.

The 9-point stencil approximation (Provatas and Elder 2010) [[6]](https://doi.org/10.1002/9783527631520) of the Laplace operator

$\nabla^2 u_{i,j}=\frac{1}{\Delta x^2}[0.5(u_{i+1,j}+u_{i-1,j}+u_{i,j+1}+u_{i,j-1})+0.25(u_{i+1,j+1}+u_{i-1,j+1}+u_{i-1,j+1}+u_{i-1,j-1})-3u_{i,j}]$

takes into account diffusion across nearest and diagonal neighbours. Diagonal neighbours have half weight due to being further from the centre. Due to its more isotropic form, it is more numerically stable and thus more suitable for rapidly varying dynamics. The condition for numerical stability is $\Delta t<\frac{\Delta x^2}{2D}$ where $D$ is the diffusion constant of the variable.

### Euler-Maruyama Method

The Euler-Maruyama method is used to numerically approximate a solution to a stochastic differential equation.

Suppose we have the stochastic differential equation of the form $dx=f(x,t)dt+\sigma dW$. The forward update equation is $x_{n+1}=x_n+f(x,t)\Delta t+\sigma\Delta W$, where $\Delta W\sim N(0,\Delta t)$.

Since we discretise the system of stochastic partial differential equations into a lattice of coupled stochastic ordinary differential equations, we can implement the Euler-Maruyama method on each grid point. 

## Diffusion-Driven Instability

Pattern formation occurs due to diffusion-driven instability in the scale-dependent feedback model and Turing model. Diffusion-driven instability occurs when the spatially homogeneous equilibrium is linearly stable to uniform perturbations but linearly unstable to spatially heterogenous perturbations, which lead to self-organised patterns stabilised by local nonlinearities.

The homogeneous equilibrium is linearly stable to uniform perturbations when the Jacobian matrix of the system without the presence of diffusion evaluated at the homogeneous equilibrium has only negative real parts of the eigenvalues.

The homogeneous equilibrium is linearly unstable to spatially heterogeneous perturbations when the dispersion relation $|I\sigma-J+Dk^2|=0$ has positive real parts of the eigenvalues for any wavenumber $k$. The wavenumber $k_{max}$ at which the maximum positive eigenvalue real part occurs is the dominant mode and thus the characeteristic spatial frequency of the pattern.

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

## References

[[1]](https://doi.org/10.1103/PhysRevLett.90.038101) Shnerb NM, Sarah P, Lavee H, Solomon S (2003) Reactive glass and vegetation patterns. Physical Review Letters 90: 038101.

[[2]](https://doi.org/10.1016/j.ecolmodel.2006.10.005) Guttal V, Jayaprakash C (2007) Impact of noise on bistable ecological systems. Ecological Modelling 201: 420-428.

[[3]](https://doi.org/10.1016/j.tpb.2006.09.003) Kefi S, Rietkerk M, van Baalen M, Loreau M (2007) Local facilitation, bistability and transitions in arid ecosystems. Theoretical Population Biology 71: 367–379.

[[4]](https://doi.org/10.1086/342078) Rietkerk M, Boerlijst MC, van Langevelde F, HilleRisLambers R, van de Koppel J, et al. (2002) Self-organization of vegetation in arid ecosystems. American Naturalist 160: 524-530.

[[5]](https://doi.org/10.1029/2007RG000256) Borgogno, F., P. D’Odorico, F. Laio, and L. Ridolfi (2009), Mathematical models of vegetation pattern formation in ecohydrology, Rev. Geophys., 47, RG1005, doi:10.1029/2007RG000256.

[[6]](https://doi.org/10.1086/659945) Dakos, V., Kéfi, S., Rietkerk, M., van Nes, E. H., & Scheffer, M. (2011). Slowing down in spatially patterned ecosystems at the brink of collapse. The American Naturalist, 177(6). doi:10.1086/659945

[[7]](https://doi.org/10.1002/9783527631520) Provatas, Nikolas; Elder, Ken (2010-10-13). Phase-Field Methods in Materials Science and Engineering. Weinheim, Germany: Wiley-VCH Verlag GmbH & Co. KGaA. p. 219. doi:10.1002/9783527631520. ISBN 978-3-527-63152-0.
