import numpy as np
import pandas as pd
import blosc
import pickle
from torch.autograd.functional import jacobian
from torch import tensor
from spatial_indicators import coarse_graining
from scipy.optimize import fsolve, least_squares


class LocalPositiveFeedbackModel:
    """ The local positive feedback model is based on the vegetation model by Shnerb et 
        al. (2003) [1] and Guttal and Jayaprakash (2007) [2]:

        .. math::
            \frac{dw_{i,j}}{dt}=R-w_{i,j}-\lambda w_{i,j}B_{i,j}+D(w_{i+1,j}+w_{i-1,j}+w_{i,j+1}+w_{i,j-1}-4w_{i,j})+\sigma_w dW_{i,j}
        .. math::
            \frac{dB_{i,j}}{dt}=\rho B_{i,j}\left(w_{i,j}-\frac{B_{i,j}}{B_c}\right)-\mu\frac{B_{i,j}}{B_{i,j}+B_0}+D(B_{i+1,j}+B_{i-1,j}+B_{i,j+1}+B_{i,j-1}-4B_{i,j})+\sigma_B dW_{i,j}

        Space is discretised into a 2D lattice of coupled cells with periodic boundaries.
        
        Vegetation biomass `B` grows logistically based on water availability `w` and 
        is lost due to grazing.

        Biomass and water are exchanged between neighbouring cells. A patch with high biomass 
        will diffuse biomass to its neighbours, leading to a positive effects on its neighbours 
        but a negative effect on the site itself. This results in no patches and a relatively 
        homogeneous state.

        When rainfall falls below a certain value, the cells undergo a shift 
        to a desert state through synchronised fold bifurcations.

        References
        ----------
        .. [1] Shnerb NM, Sarah P, Lavee H, Solomon S (2003) Reactive glass and vegetation patterns. Physical
           Review Letters 90: 038101.

        .. [2] Guttal V, Jayaprakash C (2007) Impact of noise on bistable ecological systems. Ecological Modelling
           201: 420-428.
    """
    def __init__(self, D=0.05, lambda_=0.12, rho=1., Bc=10., mu=2., B0=1., R=2., sigma_w=0.01, sigma_B=0.25, size=100, dt=0.1, noise="additive"):
        """ Construct a local positive feedback model object.

            Parameters
            ----------
            D : float, default=0.05
                Exchange (diffusion) rate, units in day^-1.
            lambda_ : float, default=0.12
                Water consumption rate by vegetation, units in g^-1 day^-1.
            rho : float, default=1
                Maximum vegetation growth rate, units in day^-1.
            Bc : float, default=10
                Vegetation carrying capacity, units in g.
            mu : float, default=2
                Maximum grazing rate, units in day^-1.
            B0 : float, default=1
                Half-saturation constant of vegetation consumption, dimensionless units.
            R : float, default=2
                Mean annual rainfall, units in mm day^-1.
            sigma_w : float, default=0.01
                Standard deviation of white noise on water moisture, dimensionless units.
            sigma_B : float, default=0.25
                Standard deviation of white noise on vegetation biomass, dimensionless units.
            size : int, default=100
                Number of grid points on each side of square 2D grid.
            dt : float, default=0.1
                Time step size, units in days.
            noise : str
                Type of noise (`"additive"` or `"multiplicative"`).
        """
        self.__dict__.update(locals())

    def biomass_reaction(self, B, w):
        """ Calculate the local reaction rate of biomass at each grid point.

            Parameters
            ----------
            B : ndarray
                2D grid of vegetation biomass, units in g.
            w : ndarray
                2D grid of water moisture level, units in mm.

            Returns
            -------
            ndarray
                2D grid of local reaction rate of vegetation biomass, units in g day^-1.
        """
        return self.rho*B*(w - B/self.Bc) - self.mu*B/(B + self.B0)

    def water_reaction(self, B, w):
        """ Calculate the local reaction rate of water moisture level at each grid point.

            Parameters
            ----------
            B : ndarray
                2D grid of vegetation biomass, units in g.
            w : ndarray
                2D grid of water moisture level, units in mm.

            Returns
            -------
            ndarray
                2D grid of local reaction rate of water moisture level, units in mm day^-1.
        """
        return self.R - w - self.lambda_*w*B
    
    def time_step_update(self, B, w):
        """ Update vegetation biomass and water moisture level in one time step at each 
            grid point using the Euler-Maruyama method.

            Parameters
            ----------
            B : ndarray
                2D grid of vegetation biomass, units in g.
            w : ndarray
                2D grid of water moisture level, units in mm.

            Notes
            -----
            The Euler-Maruyama method is used an an approximate numerical solution of a 
            stochastic differential equation by discretising time. This is an explicit numerical scheme.

            .. math::
                w_{i,j}^{n+1}=w_{i,j}^n+(R-w_{i,j}^n-\lambda w_{i,j}^n B_{i,j}+D(w_{i+1,j}^n+w_{i-1,j}^n+w_{i,j+1}^n+w_{i,j-1}^n-4w_{i,j}^n))\Delta t+\sigma_w\Delta W_{i,j}^n
            .. math::
                B_{i,j}^{n+1}=B_{i,j}^n+\left(\rho B_{i,j^n}\left(w_{i,j}^n-\frac{B_{i,j}^n}{B_c}\right)-\mu\frac{B_{i,j}^n}{B_{i,j}^n+B_0}+D(B_{i+1,j}^n+B_{i-1,j}^n+B_{i,j+1}^n+B_{i,j-1}^n-4B_{i,j}^n)\right)\Delta t+\sigma_B\Delta W_{i,j}^n
            .. math::
                W\sim N(0,\Delta t)
        """
        # calculate deterministic change in biomass and water level
        dB = (self.biomass_reaction(B, w) + self.D*laplacian5(B))*self.dt
        dw = (self.water_reaction(B, w) + self.D*laplacian5(w))*self.dt

        # update biomass and water level
        B += dB
        w += dw

        # add white noise to biomass and water level
        add_noise(self, self.sigma_B, B)
        add_noise(self, self.sigma_w, w)

    def initialisation(self):
        """ Initialise grid of vegetation biomass and water moisture level.

            Returns
            -------
            B : ndarray
                2D grid of initial fully vegetated state of vegetation biomass, units in g.
            w : ndarray
                2D grid of initial water moisture level, units in mm.
        """
        B = self.Bc*np.ones((self.size,self.size))
        w = self.R*np.ones((self.size,self.size))

        return B, w

    def equilibria(self):
        """ Find equilibria of system numerically.

            Returns
            -------
            B : float
                Biomass equilibrium, units in g.
            w : float
                Water level equilibrium, units in mm.
        """
        def f(z):
            return [self.biomass_reaction(*z), self.water_reaction(*z)]
        z, _, ier, _ = fsolve(f, [self.Bc, self.R], full_output=True)

        if ier == 1:
            B, w = tuple(z)
        else:
            B, w = (0., 0.)

        return B, w

    def reactions(self, B, w):
        """ Return reaction rates of vegetation biomass and water level.

            Parameters
            ----------
            B : float
                Vegetation biomass, units in g.
            w : float
                Water moisture level, units in mm.

            Returns
            -------
            tuple of floats
                Reaction rates of vegetation biomass and water level.
        """
        return (self.biomass_reaction(B, w), self.water_reaction(B, w))


class LocalFacilitationModel:
    """ The local facilitation model is a stochastic cellular automaton (Kefi et al. 2007) [1]
        with discrete time steps.
    
        The ecosystem is represented by a 2D lattice of cells (with periodic boundaries), each of 
        which can be in one of three states: vegetated (1), empty (0) or degraded (-1).

        Each cell can undergo one of four transformations with a certain probability at each time step:
        1. Colonisation of an empty cell to a vegetated cell
        .. math:: w_{[0,1]}=[\delta\rho_1 +(1-\delta)q_{1|0}](b-c\rho_1)
        2. Mortality of a vegetated cell to an empty cell
        .. math:: w_{[1,0]}=m
        3. Degradation of a empty cell to a degraded cell
        .. math:: w_{[0,-1]}=d
        4. Regeneration of a degraded cell to an empty cell
        .. math:: w_{[-1,0]}=r+fq_{1|-1}

        Colonisation and regeneration probabilities of a cell are positively affected by presence of 
        vegetation in neighbouring cells.

        The facilitation leads to the formation of patches of vegetated grid cells; these patches have 
        size distributions that follow a power law.

        References
        ----------
        .. [1] Kefi S, Rietkerk M, van Baalen M, Loreau M (2007) Local facilitation, bistability and transitions
               in arid ecosystems. Theoretical Population Biology 71: 367–379.
    """
    def __init__(self, m=0.1, f=0.9, delta=0.1, c=0.3, r=0.0001, d=0.2, b=1, size=100):
        """ Construct a local facilitation model object.

            Parameters
            ----------
            m : float, default=0.1
                Mortality probability of a vegetated site.
            f : float, default=0.9
                Local facilitation strength; maximum effect of a neighbouring vegetation site on the 
                regeneration of a degraded site.
            delta : float, default=0.1
                Fraction of seeds globally dispersed.
            c : float, default=0.3
                Intrinsic seed production rate per vegetated site (germination probability) multiplied 
                by the global competitive effect.
            r : float, default=0.0001
                Regeneration probability of a degraded site without vegetated sites in its neighbourhood.
            d : float, default=0.2
                Degradation probability of empty sites.
            b : float, default=1
                Severity of environmental conditions; a lower `b` value reflects a high aridity level.
            size : int, default=100
                Number of grid points on each side of square 2D grid.
        """
        self.__dict__.update(locals())

    def cell_update(self, V, p, q, rho):
        """ Update the state of a cell.

            Parameters
            ----------
            V : int
                Cell in one of three states; vegetated (1), empty (0) or degraded (-1).
            p : float
                Random uniform number between 0 and 1.
            q : float
                Density of vegetated neighbouring sites of the cell.
            rho : float
                Global density of vegetated cells.

            Returns
            -------
            int
                Updated state of cell.
        """
        # cell is vegetated
        if V == 1:
            # moratility
            if p < self.m:
                return 0

        # cell is empty
        elif V == 0:
            # degradation
            if p < self.d:
                return -1

            # colonisation
            elif p < self.d + (self.delta*rho + (1-self.delta)*q)*(self.b-self.c*rho):
                return 1

        # cell is degraded
        else:
            # regeneration
            if p < self.r + self.f*q:
                return 0

        # state of cell remains the same
        return V

    def time_step_update(self, V):
        """ Update state of each cell after one time step.

            Parameters
            ----------
            V : ndarray
                2D array of cell states; vegetated (1), empty (0) and degraded (-1).
        """
        # global vegetation density
        rho = np.count_nonzero(V == 1) / (self.size**2)

        # local vegetation density at each cell
        q = adjacent_neighbour_sum(np.where(V == 1, 1, 0)) / 4

        # grid of uniform random numbers between 0 and 1
        random = np.random.rand(self.size, self.size)

        # update state of each cell
        V[...] = np.vectorize(self.cell_update)(V, random, q, rho)

    def initialisation(self):
        """ Initialise fully vegetated grid of cells.

            Returns
            -------
            ndarray
                2D grid of vegetated (1) cells.
        """
        return np.ones((self.size, self.size)),

    def vegetation_cover(self, snapshots, s=4):
        """ Calculate vegetation cover as a fraction of vegetated cells over each sxs group of cells in each snapshot.

            Parameters
            ----------
            snapshots : Series
                Spatial snapshots at each control parameter level.
            s : int, default=4
                Coarse-graining submatrix size.

            Returns
            -------
            Series
                Vegetation cover for each snapshot.
        """
        cover = []

        for snapshot in snapshots:
            cover.append(coarse_graining(np.where(snapshot == 1, 1, 0), s))

        return pd.Series(cover, index=snapshots.index)

    def vegetation_patchiness(self, snapshots):
        """ Calculate vegetation patchiness in each snapshot.

            Parameters
            ----------
            snapshots : Series
                Spatial snapshots at each control parameter level.

            Returns
            -------
            Series
                Vegetation patchiness for each snapshot.
        """
        patchiness = []

        N = self.size**2  # number of cells in the grid

        for snapshot in snapshots:
            snapshot = np.where(snapshot == 1, 1, 0)

            # calculate the global vegetation density
            global_vegetation = np.sum(snapshot) / N

            # get neighbours
            top = np.roll(snapshot, (0, 1), (0, 1))
            bottom = np.roll(snapshot, (0, -1), (0, 1))
            right = np.roll(snapshot, (1, 0), (0, 1))
            left = np.roll(snapshot, (-1, 0), (0, 1))

            # calculate number of vegetated neighbours if the cell is vegetated
            vegetated_neighbours = np.where(snapshot == 1, top + bottom + right + left, 0)

            # calculate the density of vegetated pairs
            vegetated_pairs = np.sum(vegetated_neighbours) / (4*N)

            # calculate the patchiness
            if global_vegetation != 0:
                patchiness.append(vegetated_pairs / global_vegetation**2)
            else:
                patchiness.append(np.NaN)

        return pd.Series(patchiness, index=snapshots.index)

    def reactions(self, V, D):
        """ Return the reaction rates of the density of vegetated cells and degraded cells of the mean-field approximation.

            Parameters
            ----------
            V : float
                Density of vegetated cells.
            D : float
                Density of degraded cells.

            Returns
            -------
            tuple of floats
                Reaction rates of the density of vegetated cells and degraded cells, respectively.
        """
        return (
            V*(self.b - self.c*V)*(1 - V - D) - self.d*V,
            self.d*(1 - V - D) - (self.r + self.f*V)*D
        )

    def equilibria(self, approximation="mean"):
        """ Find equilibria of system numerically by mean field model or pair approximation model.

            Parameters
            ----------
            approximation : str, default="mean"
                Approximation model to use, `"mean"` for mean field model or `"pair"` for pair approximation model.

            Notes
            -----
            The mean field does not take into account local interactions within the model so is less accurate.

            The implementation of the pair approximation model is not working.
        """
        # TODO: Fix this
        if approximation == "pair":
            def f(z):
                V = z[0]  # density of vegetated cells
                D = z[1]  # density of degraded cells
                VV = z[2]  # density of vegateted-vegetated pairs
                VD = z[3]  # density of vegetated-degraded pairs
                DD = z[4]  # density of degraded-degraded pairs

                E = 1 - V - D  # density of empty cells
                VE = V - VV - DD  # density of vegetated-empty pairs
                ED = D - DD - VD  # density of empty-degraded pairs

                qVD = VD / D  # local density of vegetated neighbours given degraded cell
                qVE = VE / E  # local density of vegetated neighbours given empty cell

                return [
                    # reaction rate of vegetated cells
                    (self.delta*V + (1-self.delta)*qVE)*(self.b-self.c*V)*E - self.m*V,
                    # reaction rate of degraded cells
                    self.d*E - (self.r + self.f*qVD)*D,
                    # reaction rate of vegetated-vegetated pairs
                    2*VE*(self.delta*V + (1-self.delta)/4 + 3/4*(1-self.delta)*qVE)*(self.b-self.c*V) - 2*VV*self.m,
                    # reaction rate of vegetated-degraded pairs
                    self.d*VE + ED*(self.delta + 3/4*(1-self.delta)*qVE)*(self.b-self.c*V) - VD*(self.r + self.f/4 + 3/4*self.f*qVD + self.m),
                    # reaction rate of degraded-degraded pairs
                    2*self.d*ED - 2*DD*(self.r + 3/4*self.f*qVD)
                ]

            z0 = [1, 0, 1, 0, 0]  # initial guess

        elif approximation == "mean":
            def f(z):
                V, D = z  # density of vegetated and degraded cells
                return list(self.reactions(V, D))  # reaction rates of vegetated and degraded cells

            z0 = [1, 0]  # initial guess

        # calculate the variable values with least squares
        res = least_squares(f, z0, bounds=([0]*len(z0), [1]*len(z0)))
        z = res.x

        # return density of vegetated and degraded cells only
        V, D = tuple(z[:2])

        return V, D


class ScaleDependentFeedbackModel:
    """ The scale-dependent feedback model is a system of stochastic partial differential equations 
        describing the dynamics of vegetation biomass `P`, soil water `W` and surface water `O` 
        (Rietkerk et al. 2002) [1].
        
        .. math::
            \frac{\partial O}{\partial t}=R(t)-\alpha O\frac{P+W_0 k_2}{P+k_2}+D_O\nabla^2O+\sigma dW
        .. math::
            \frac{\partial W}{\partial t}=\alpha O\frac{P+W_0 k_2}{P+k_2}-g_{max}\frac{W}{W+k_1}P-r_W W+D_W\nabla^2W+\sigma dW
        .. math::
            \frac{\partial P}{\partial t}=\left(cg_{max}\frac{W}{W+k_1}-d\right)P+D_P\nabla^2P+\sigma dW

        Space is discretised into a 2D lattice of coupled cells with periodic boundaries.

        Plants grow depending on soil water availability and lost from mortality or grazing. Surface water is 
        supplied by rainfall and lost due to soil infiltration and runoff. Soil water due to soil infiltration is taken 
        up by plants or lost by runoff.

        Plants, soil water and surface water all diffuse in 2D space.

        The infiltration rate of water in the soil is higher in areas with vegetation, leading to accumulation of water under 
        vegetation and to its depletion further away (scale-dependent feedback).

        At the Turing instability, the scale-dependent feedback leads to the formation of regular vegetation patterns.

        References
        ----------
        .. [1] Rietkerk M, Boerlijst MC, van Langevelde F, HilleRisLambers R, van de Koppel J, et al. (2002)
               Self-organization of vegetation in arid ecosystems. American Naturalist 160: 524-530.
    """
    def __init__(self, c=10., gmax=0.05, k1=5., d=0.25, alpha=0.2, k2=5., W0=0.2, rW=0.2, DP=0.1, DW=0.1, DO=100, R=2, sigma=0.01, size=100, dt=0.1, dx=5, noise="additive"):
        """ Construct a scale-dependent feedback model object.

            Parameters
            ----------
            c : float, default=10
                Conversion factor for water uptake to plant biomass, units in g mm^-1 m^-2.
            gmax : float, default=0.05
                Maximum specific water uptake, units in mm g^-1 m^-2 day^-1.
            k1 : float, default=5
                Half saturation constant of water uptake by plants, units in mm.
            d : float, default=0.25
                Specific rate of plant density loss due to mortality, units in day^-1.
            alpha : float, default=0.2
                Rate of surface water infiltration, units in day^-1.
            k2 : float, default=5
                Plant density scale determining how surface water infiltration increases with P, units in g m^-2.
            W0 : float, default=0.2
                Minimum surface water infiltration coefficient in the absence of plants, dimensionless units.
            rW : float, default=0.2
                Soil water loss rate due to evaporation and drainage, units in day^-1.
            DP : float, default=0.1
                Plant dispersal diffusion constant, units in m^2 day^-1.
            DW : float, default=0.1
                Soil water diffusion constant, units in m^2 day^-1.
            DO : int, default=100
                Surface water diffusion constant, units in m^2 day^-1.
            R : float, default=2
                Mean annual rainfall, units in mm.
            sigma : float, default=0.01
                Standard deviation of white noise, dimensionless units.
            size : int, default=100
                Number of grid points on each side of square 2D grid.
            dt : float, default=0.1
                Time step size, units in days.
            dx : float, default=5
                Spatial step size, units in m.
            noise : str, default="additive"
                Type of noise (`"additive"` or `"multiplicative"`).
        """
        self.__dict__.update(locals())

    def P_reaction(self, P, W, O):
        """ Calculate the local reaction rate of plant biomass at each grid point.

            Parameters
            ----------
            P : ndarray
                2D grid of plant density, units in g m^-2.
            W : ndarray
                2D grid of soil water level, units in mm.
            O : ndarray
                2D grid of surface water level, units in mm.

            Returns
            -------
            ndarray
                2D grid of local reaction rate of plant density, units in g m^-2 day^-1.
        """
        return self.c*self.gmax*W/(W + self.k1)*P - self.d*P

    def W_reaction(self, P, W, O):
        """ Calculate the local reaction rate of soil water level at each grid point.

            Parameters
            ----------
            P : ndarray
                2D grid of plant density, units in g mm^-2.
            W : ndarray
                2D grid of soil water level, units in mm.
            O : ndarray
                2D grid of surface water level, units in mm.

            Returns
            -------
            ndarray
                2D grid of local reaction rate of soil water level, units in mm day^-1.
        """
        return self.alpha*O*(P + self.k2*self.W0)/(P + self.k2)\
             - self.gmax*W/(W + self.k1)*P - self.rW*W

    def O_reaction(self, P, W, O):
        """ Calculate the local reaction rate of surface water level at each grid point.

            Parameters
            ----------
            P : ndarray
                2D grid of plant density, units in g mm^-2.
            W : ndarray
                2D grid of soil water level, units in mm.
            O : ndarray
                2D grid of surface water level, units in mm.

            Returns
            -------
            ndarray
                2D grid of local reaction rate of surface water level, units in mm day^-1.
        """
        return self.R - self.alpha*O*(P + self.k2*self.W0)/(P + self.k2)

    def time_step_update(self, P, W, O):
        """ Update plant biomass, soil water level and surface water level in one time step 
            at each grid point using the finite-difference method and Euler-Maruyama method.

            Parameters
            ----------
            P : ndarray
                2D grid of plant density, units in g m^-2.
            W : ndarray
                2D grid of soil water level, units in mm.
            O : ndarray
                2D grid of surface water level, units in mm.

            Notes
            -----
            The finite difference method is used to calculate the spatial derivatives (Laplacian) 
            of the variables. The 5-point stencil is used for plant biomass `P` and soil water `W`. 
            Due to the fast dynamics of surface water `O`, the 9-point stencil is used to reduce 
            numerical instability.

            The Euler-Maruyama method is used an an approximate numerical solution of a 
            stochastic differential equation by discretising time. This is an explicit numerical scheme.

            See Also
            --------
            laplacian5 : Laplacian using 5-point stencil (considering adjacent neighbours).
            laplacian9 : Laplacian using 9-point stencil (considering adjacent and diagonal neighbours).
        """
        # calculate deterministic change in variables
        dP = (self.P_reaction(P, W, O) + self.DP*laplacian5(P, self.dx))*self.dt
        dW = (self.W_reaction(P, W, O) + self.DW*laplacian5(W, self.dx))*self.dt
        dO = (self.O_reaction(P, W, O) + self.DO*laplacian9(O, self.dx))*self.dt
        
        # update variables
        P += dP
        W += dW
        O += dO
        
        # add white noise to varaibles
        add_noise(self, self.sigma, P, W, O)

    def equilibria(self):
        """ Find the homogeneous equilibria of the system.

            Returns
            -------
            P : float
                Homogeneous plant equilibrium, units in g mm^-2. If negative set to zero.
            W : float
                Homogeneous soil water equilibrium, units in mm.
            O : float
                Homogeneous surface water equilibrium, units in mm.

            Notes
            -----
            The homogeneous equilibria occur when diffusion and local reaction of each variable is zero.

            .. math:: W=\frac{dk_1}{cg_{max}-d}
            .. math:: P=c\frac{R-r_W W}{d}
            .. math:: O=\frac{R}{\alpha}\frac{P+k_2}{P+W_0 k_2}
        """
        # vegetated equilibrium
        W = self.k1*self.d/(self.c*self.gmax - self.d)
        P = self.c/self.d*(self.R - self.rW*W)
        O = self.R/self.alpha*(P + self.k2)/(P + self.W0*self.k2)

        # plantless equilibrium
        if P <= 0:
            P = 0.
            W = self.R/self.rW
            O = self.R/(self.alpha*self.W0)

        return P, W, O
    
    def initialisation(self):
        """ Initialise grids of plant biomass, soil water level and surface water level as homogeneous equilibria.

            Returns
            -------
            P : ndarray
                2D grid of initial plant biomass, units in g m^-2.
            W : ndarray
                2D grid of initial soil water level, units in mm.
            O : ndarray
                2D grid of initial surface water level, units in mm.
        """
        # find homogeneous equilibria
        P, W, O = self.equilibria()

        # set 2D grid values
        P = P*np.ones((self.size, self.size))
        W = W*np.ones((self.size, self.size))
        O = O*np.ones((self.size, self.size))

        return P, W, O
    
    def diffusion_matrix(self):
        """ Construct a diffusion matrix with the diffusion constants of the variables on the diagonals.

            Returns
            -------
            ndarray
                3x3 diffusion matrix, with plant dispersal diffusion constant on the first diagonal, 
                soil water diffusion constant on the second diagonal and surface water diffusion constant 
                on the thrid diagonal.
        """
        return np.array([[self.DP, 0, 0], [0, self.DW, 0], [0, 0, self.DO]])

    def reactions(self, P, W, O):
        """ Return the reaction rates of plant density, soil water and surface water.

            Parameters
            ----------
            P : float
                Plant density, units in g mm^-2.
            W : float
                Soil water level, units in mm.
            O : float
                Surface water level, units in mm.

            Returns
            -------
            tuple of floats
                Reaction rates of plant density, soil water and surface water.
        """
        return (self.P_reaction(P, W, O), self.W_reaction(P, W, O), self.O_reaction(P, W, O))


class TuringModel:
    """ The Turing model is a system of stochastic version of the partial differential equation system 
        with two species, the activator `u` and the inhibitor `v` (Borgogno et al. 2009) [1].

        .. math::
            \frac{\partial u}{\partial t}=u(avu-e)+\nabla^2u+\sigma dW
        .. math::
            \frac{\partial v}{\partial t}=v(b-cu^2v)+d\nabla^2v+\sigma dW
        
        The growth rate of species `u` increases with increasing values of `u` and `v`.

        The species `u` has a strong negative influence (inhibition) on the growth rate of `v`.

        References
        ----------
        .. [1] Borgogno, F., P. D’Odorico, F. Laio, and L. Ridolfi (2009), Mathematical models of vegetation pattern formation in
               ecohydrology, Rev. Geophys., 47, RG1005, doi:10.1029/2007RG000256.
    """
    def __init__(self, a=22, b=84, c=113.33, d=27.2, e=18, sigma=0.01, size=256, dt=0.1, dx=1, noise="additive"):
        """ Construct a Turing model object.

            Parameters
            ----------
            a : int, default=22
            b : int, default=84
            c : float, default=113.33
            d : float, default=27.2
            e : int, default=18
            sigma : float, default=0.01
                Standard deviation of white noise.
            size : int, default=256
                Number of grid points on each side of square 2D grid.
            dt : float, default=0.1
                Time step size.
            dx : float, default=1
                Spatial step size.
            noise : str
                Type of noise (`"additive"` or `"multiplicative"`).
        """
        self.__dict__.update(locals())

    def u_reaction(self, u, v):
        """ Calculate the local reaction rate of the activator variable u at each grid point.

            Parameters
            ----------
            u : ndarray
                2D grid of activator variable.
            v : ndarray
                2D grid of inhibitor variable.

            Returns
            -------
            ndarray
                2D grid of local reaction rate of activator variable u.
        """
        return u*(self.a*v*u - self.e)

    def v_reaction(self, u, v):
        """ Calculate the local reaction rate of the inhibitor variable v at each grid point.

            Parameters
            ----------
            u : ndarray
                2D grid of activator variable.
            v : ndarray
                2D grid of inhibitor variable.

            Returns
            -------
            ndarray
                2D grid of local reaction rate of inhibitor variable v.
        """
        return v*(self.b - self.c*u**2*v)

    def time_step_update(self, u, v):
        """ Update activator variable u and inhibitor variable v in one time step at each grid point.

            Parameters
            ----------
            u : ndarray
                2D grid of activator variable.
            v : ndarray
                2D grid of inhibitor variable.

            Notes
            -----
            The finite difference method is used to calculate the spatial derivatives (Laplacian) 
            of the variables. The 5-point stencil is used for the activator `u` 
            Due to the fast dynamics of the inhibitor `v`, the 9-point stencil is used to reduce 
            numerical instability.

            The Euler-Maruyama method is used an an approximate numerical solution of a 
            stochastic differential equation by discretising time. This is an explicit numerical scheme.

            See Also
            --------
            laplacian5 : Laplacian using 5-point stencil (considering adjacent neighbours).
            laplacian9 : Laplacian using 9-point stencil (considering adjacent and diagonal neighbours).
        """
        # calculate deterministic change in u and v
        du = (self.u_reaction(u, v) + laplacian5(u, self.dx))*self.dt
        dv = (self.v_reaction(u, v) + self.d*laplacian9(v, self.dx))*self.dt

        # update u and v
        u += du
        v += dv

        # add white noise to u and v
        add_noise(self, self.sigma, u, v)

    def equilibria(self):
        """ Find the homogeneous equilibria of the system.

            Returns
            -------
            u : float
                Homogenous `u` equilibrium.
            v : float
                Homogeneous `v` equilibrium.
        """
        u = self.a*self.b/(self.c*self.e)
        v = self.c*self.e**2/(self.b*self.a**2)

        return u, v

    def initialisation(self):
        """ Initialise grids of u and v as homogeneous equilibria.

            Returns
            -------
            u : ndarray
                2D grid of initial u values.
            v : ndarray
                2D grid of initial v values.
        """
        u, v = self.equilibria()

        u = u*np.ones((self.size, self.size))
        v = v*np.ones((self.size, self.size))

        return u, v

    def diffusion_matrix(self):
        """ Construct a diffusion matrix with the diffusion constants on the diagonals.

            Returns
            -------
            ndarray
                2x2 array, with the diffusion constant for `u` on the first diagonal and the diffusion constant for 
                `v` on the second diagonal.
        """
        return np.array([[1, 0], [0, self.d]])

    def reactions(self, u, v):
        """ Return the reaction rates of `u` and `v`.

            Parameters
            ----------
            u : float
                Activator variable.
            v : float
                Inhibitor variable
            
            Returns
            -------
            tuple of floats
                Reaction rates of `u` and `v`.
        """
        return (self.u_reaction(u, v), self.v_reaction(u, v))


def adjacent_neighbour_sum(array):
    """ Calculate the sum of the four adjacent neighbours to each cell in an array.

        Parameters
        ----------
        array : ndarray
            2D grid of numeric values.
        
        Returns
        -------
        ndarray
            2D grid with sum of four adjacent neighbours in each cell of original array.
    """
    # shift arrays horizontally and vertically
    top = np.roll(array, (0, 1), (0, 1))
    bottom = np.roll(array, (0, -1), (0, 1))
    right = np.roll(array, (1, 0), (0, 1))
    left = np.roll(array, (-1, 0), (0, 1))

    return top + bottom + right + left


def diagonal_neighbour_sum(array):
    """ Calculate the sum of the four diagonal neighbours in each cell.

        Parameters
        ----------
        array : ndarray
            2D grid of numeric values.
        
        Returns
        -------
        ndarray
            2D grid with sum of four diagonal neighbours in each cell of original array.
    """
    # shift arrays diagonally
    top_right = np.roll(array, (1, 1), (0, 1))
    top_left = np.roll(array, (-1, 1), (0, 1))
    bottom_right = np.roll(array, (1, -1), (0, 1))
    bottom_left = np.roll(array, (-1, -1), (0, 1))

    return top_right + top_left + bottom_right + bottom_left


def laplacian5(array, dx=1):
    """ Calculate the Laplacian to each grid point of an array according to the five-point stencil.

        Parameters
        ----------
        array : ndarray
            2D grid of numeric values.
        dx : int or float, default=1
            Spatial step size.

        Returns
        -------
        ndarray
            2D grid of Laplacian at each grid point of original array.

        Notes
        -----
        The Laplacian approximation is obtained by the finite-difference method.

        The 5-point stencil only takes into account the top/bottom and left/right neighbours.

        The 5-point stencil is stable for very smoothly varying fields.

        See Also
        --------
        laplacian9 : Laplacian using 9-point stencil (more numerically stable).
    """
    return (adjacent_neighbour_sum(array) - 4*array) / dx**2


def laplacian9(array, dx=1):
    """ Calculate the Laplacian to each grid point of an array according to the nine-point stencil.
        
        Parameters
        ----------
        array : ndarray
            2D grid of numeric values.
        dx : int or float, default=1.
            Spatial step size.
        
        Returns
        -------
        ndarray
            2D grid of Laplacian at each grid point of original array.

        Notes
        -----
        The 9-point stencil considers the adjacent neighbours and diagonal neighbours. By summing 
        the Taylor series of the neighbours, with the diagonal neighbours each weighted by half, we 
        obtain the 9-point stencil [1].

        The 9-point stencil is a more stable and isotropic form of the Laplacian operator compared to the
        5-point stencil. Thus it is more useful for rapidly varying fields.

        See Also
        --------
        laplacian5 : Laplacian using 5-point stencil.

        References
        ----------
        .. [1] Provatas, Nikolas; Elder, Ken (2010-10-13). Phase-Field Methods in Materials Science and Engineering. 
               Weinheim, Germany: Wiley-VCH Verlag GmbH & Co. KGaA. p. 219. doi:10.1002/9783527631520. ISBN 978-3-527-63152-0.
    """
    return (0.5*adjacent_neighbour_sum(array) + 0.25*diagonal_neighbour_sum(array) - 3*array) / dx**2


def add_noise(model, sigma, *arrays):
    """ Add white noise to arrays of variables in model.

        Parameters
        ----------
        model : obj
            Model object.
        sigma : float
            Standard deviation of noise on variable.
        *arrays : tuple
            Arrays of variables to add noise.
    """
    generate_noise = lambda: sigma*np.sqrt(model.dt)*np.random.normal(size=(model.size, model.size))

    if model.noise == "additive":
        for array in arrays:
            array += generate_noise()
    elif model.noise == "multiplicative":
        for array in arrays:
            array += array*generate_noise()


def simulate_time_steps(model, time_steps, *u):
    """ Simulate a number of time steps on the variables of a model and save the snapshots of the first variable.

        Parameters
        ----------
        model : obj
            Model object.
        time_steps : int
            Number of time steps to simulate.
        *u : tuple
            Arrays of variables to simulate, the first of which to save the snapshots of.
        
        Returns
        -------
        Series
            Snapshots of first variable at each time step.
    """
    snapshots = []

    dt = model.dt if hasattr(model, "dt") else 1  # set time step
    t = dt*np.arange(time_steps)  # range of time values

    # loop through each time step
    for _ in t:
        model.time_step_update(*u)
        snapshots.append(u[0].copy())

    return pd.Series(snapshots, index=t, dtype="object")


def parameter_change(model, p, parameter_name, time_steps=100, warm_up=0, u=None):
    """ Simulate the model a number of time steps at each control parameter level.

        Parameters
        ----------
        model : obj
            Model object.
        p : iterable
            Range of control parameter values.
        parameter_name : str
            Name of parameter, must match attribute in model.
        time_steps : int, default=100
            Number of time steps to simulate at each parameter level.
        warm_up : int, default=0.
            Number of time steps to simulate on initial variable values.
        u : tuple of ndarrays, optional
            Initial variable values.
        return_final : bool, default=False
            Whether to return final variable values.

        Returns
        -------
        snapshots : Series
            Snapshots of first variable at each control parameter level.
        u : tuple of ndarrays
            Variable values at end of final simulation.
        
    """
    # set parameter to first value in parameter range
    setattr(model, parameter_name, p[0])
    
    # initialise variables
    if u is None:
        u = model.initialisation()
    
    snapshots = []  # initialise snapshots list

    simulate_time_steps(model, warm_up, *u)  # warm up

    # iterate through each control parameter value
    for p_step in p:
        setattr(model, parameter_name, p_step)  # set control parameter level
        simulate_time_steps(model, time_steps, *u)
        snapshots.append(u[0].copy())  # store snapshot of first variable (discarding transients)

    return pd.Series(snapshots, index=p, dtype="object")


def bound(x, lower, upper):
    """ Bound a value between a lower and upper bound.

        Parameters
        ----------
        x : float
            Value to be bounded.
        lower : float
            Lower bound.
        upper : float
            Upper bound.

        Returns
        -------
        x : float
            Bounded value.    
    """
    if upper is not None:
        x = min(x, upper)
    if lower is not None:
        x = max(x, lower)
    return x


def preventive_measure(model, time_steps, a0, da, b0, db, a_name, b_name, preventive_point, a_bounds=None, b_bounds=None, u=None, warm_up=100):
    """ Simulate the model while changing a primary parameter level. A secondary parameter level begins changing when the primary parameter passes a certain value.

        Parameters
        ----------
        model : obj
            Model object.
        time_steps : int
            Total number of time steps to simulate.
        a0 : float
            Initial primary parameter value.
        da : float
            Step in primary parameter value per time step.
        b0 : float
            Initial secondary parameter value.
        db : list of floats
            Steps in secondary parameter value per time step.
        preventive_point : float
            Primary parameter value at which secondary parameter begins changing.
        a_bounds : tuple of ints, optional
            Lower and upper bounds of primary parameter.
        b_bounds : tuple of ints, optional
            Lower and upper bounds of secondary parameter.
        u : int, optional
            Initial variable values. If None, warm up the system to the equilibrium.
        warm_up : int, default=100
            Number of time steps to warm up system.

        Returns
        all_snapshots : list of Series
            List of snapshots both before preventive measures and after preventive measures at different rates of the secondary parameter.
    """
    all_snapshots = []

    # initialise variables
    if u is None:
        u = model.initialisation()
        setattr(model, a_name, a0)  # set initial primary parameter
        setattr(model, b_name, b0)  # set initial secondary parameter
    simulate_time_steps(model, warm_up, *u)  # warm up system to equilibrium

    # default parameter bounds
    if a_bounds is None:
        a_bounds = (None, None)
    if b_bounds is None:
        b_bounds = (None, None)
    
    dt = model.dt if hasattr(model, "dt") else 1  # set time step
    t = dt*np.arange(time_steps)  # range of time values

    change_idx = int((preventive_point - a0)/da)  # index at which preventive measures are placed

    snapshots = []
    a = a0  # initial primary parameter

    # simulate system before preventive measures
    for _ in t[:change_idx]:
        setattr(model, a_name, a)  # set primary parameter

        model.time_step_update(*u)
        snapshots.append(u[0].copy())

        a = bound(a + da, *a_bounds)  # update primary parameter
    
    all_snapshots.append(pd.Series(snapshots, index=t[:change_idx]))

    # iterate through each rate of secondary parameter
    for db_step in db:
        new_u = tuple(x.copy() for x in u)  # copy of variable values at preventive point

        snapshots = []
        a = preventive_point  # primary parameter at preventive point
        b = b0  # initial secondary primary
        
        # simulate system after preventive measures
        for _ in t[change_idx:]:
            setattr(model, a_name, a)  # set primary parameter
            setattr(model, b_name, b)  # set secondary parameter

            model.time_step_update(*new_u)
            snapshots.append(new_u[0].copy())

            a = bound(a + da, *a_bounds)  # update primary parameter
            b = bound(b + db_step, *b_bounds)  # update secondary parameter

        all_snapshots.append(pd.Series(snapshots, index=t[change_idx:]))

    return all_snapshots


def jacobian_matrix(model):
    """ Numerically approximate the Jacobian at the equilibrium of a system.

        Parameters
        ----------
        model : obj
            Model object.
        
        Returns
        -------
        tuple of float tensors
            Jacobian approximation.
    """
    U = model.equilibria()
    U_tensor = tuple(tensor(u) for u in U)

    return jacobian(model.reactions, U_tensor)


def max_eigenvalue(model, k=None):
    """ Find the maximum real part of eigenvalue of the Jacobian of a system.

        Parameters
        ----------
        model : obj
            Model object.
        k : iterable, optional
            Range of wavenumber values to find maximum eigenvalues at when diffusion is present.
        
        Returns
        -------
        float or Series
            Maximum real part of eigenvalue of Jacobian matrix. If `k` is set, maximum eigenvalues at each wavenumber.

        Notes
        -----
        When diffusion is not present we find the solution to :math:`|I\sigma-J|=0`. When diffusion is present
        we find the solution to :math:`|I\sigma-J+Dk^2|=0`; this is the dispersal relation.
    """
    J = jacobian_matrix(model)

    if k is not None:
        w_max = []
        for k_step in k:
            w, _ = np.linalg.eig(J - model.diffusion_matrix()*k_step**2)

            w_max.append(max(eigenvalue.real for eigenvalue in w))
        
        return pd.Series(w_max, index=k)
    
    else:
        w, _ = np.linalg.eig(J)

        return max(eigenvalue.real for eigenvalue in w)


def find_bifurcation(model, p0, pf, dp, param_name, tol=1.e-10, k=None):
    """ Find the bifurcation point of a model by checking if the maximum eigenvalue reaches zero over a parameter range.

        Parameters
        ----------
        model : obj
            Model object.
        p0 : float
            Initial parameter value.
        pf : float
            Final parameter value.
        dp : float
            Parameter step size.
        param_name : str
            Name of parameter, must match attribute in model.
        tol : float, default=1.e-10
            Error tolerance.
        k : iterable, optional
            Range of wavenumber values to find maximum eigenvalues.

        Returns
        -------
        p : float
            Parameter at which maximum eigenvalue reaches zero.
    """
    # initial maximum eigenvalue
    setattr(model, param_name, p0)
    w_max = max_eigenvalue(model, k)

    if k is not None:
        w_max = max(w_max.values)

    # iterate over parameter range
    for p in np.arange(p0 + dp, pf, dp):
        # set current maximum eigenvalue
        setattr(model, param_name, p)
        current_w_max = max_eigenvalue(model, k)

        if k is not None:
            current_w_max = max(current_w_max.values)

        # check if maximum eigenvalue is within tolerance of zero
        if current_w_max + tol >= 0:
            return p

        # check if maximum eigenvalue goes back down
        elif current_w_max < w_max:
            # try again with parameter range between current parameter and previous parameter
            return find_bifurcation(model, p - dp, p, dp/10, param_name, tol)
        
        # go to next parameter value
        else:
            w_max = current_w_max

    return None


def save_data(file, data):
    """ Save data into a compressed file.

        Parameters
        ----------
        file : str
            File path.
        data : obj
            Data to save into file.
    """
    pickled_data = pickle.dumps(data)
    compressed_pickle = blosc.compress(pickled_data)

    with open(file, 'wb') as f:
        f.write(compressed_pickle)


def load_data(file):
    """ Load data from a compressed file.

        Parameters
        ----------
        file : str
            File path.
        
        Returns
        -------
        data : obj
            Data from file.
    """
    with open(file, 'rb') as f:
        compressed_pickle = f.read()

    depressed_pickle = blosc.decompress(compressed_pickle)
    data = pickle.loads(depressed_pickle)

    return data