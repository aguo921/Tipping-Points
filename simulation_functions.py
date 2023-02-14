import numpy as np
import pandas as pd
from scipy.integrate import quad


def bound(x, lower, upper):
    """ Check a value is within bounds and return None if not.

        Parameters
        ----------
        x : float
            Input value.
        lower : float
            Lower bound.
        upper : float
            Upper bound.
        
        Returns
        -------
        x : float
            Input value if within bounds, None otherwise.
    """
    if x > lower and x < upper:
        return x
    else:
        return None

def sample_path(f, x0, y0, T, h, sigma, epsilon, bounds=(-np.inf, np.inf)):
    """ Generate a sample path of x values over slow parameter y change using the Euler-Maruyama method.

        Parameters
        ----------
        f : callable
            Deterministic function `f(x,y)`.
        x0 : float
            Initial `x` value (fast variable).
        y0 : float
            Initial `y` value (slow variable).
        T : float
            Run duration.
        h : float
            Time step size.
        sigma : float
            Standard deviation of white noise.
        epsilon : float
            Time scale separation.
        bounds : tuple of floats, default=(-np.inf, np.inf)
            Lower and upper bounds of `x` value.

        Returns
        -------
        Series
            Sample path of x values at each y value level.
    """
    t = np.arange(0, T, h)  # time range
    y = y0 + epsilon*t  # range of y values

    N = len(t) - 1  # number of iterations

    # initiailise x vector
    x = np.zeros(N+1)
    x[0] = x0

    # generate random noise
    dW = np.random.normal(0, np.sqrt(h), N)

    # loop through each iteration
    for i in range(N):
        if x[i] is None:
            x[i+1] = None
        else:
            # update x
            dx = f(x[i],y[i])*h + sigma*dW[i]
            x[i+1] = bound(x[i] + dx, *bounds)

    return pd.Series(x, index=y)


def ensemble_paths(f, x0, y0, T, h, sigma, epsilon, runs=100, bounds=(-np.inf, np.inf)):
    """ Generate an ensemble of sample paths.

        Parameters
        ----------
        f : callable
            Deterministic function `f(x,y)`.
        x0 : float
            Initial `x` value.
        y0 : float
            Initial `y` value.
        T : float
            Run duration.
        h : float
            Time step size.
        sigma : float
            Standard deviation of white noise.
        epsilon : float
            Time scale separation.
        runs : int, default=100
            Number of sample paths to generate.
        bounds : tuple of floats, default=(-np.inf, np.inf)
            Lower and upper bounds of `x` value.

        Returns
        -------
        Dataframe
            Dataframe where columns are each sample path run, row index is parameter y value.
    """
    paths = {}  # initialise dictionary of paths

    # generate paths
    for i in range(runs):
        path = sample_path(f, x0, y0, T, h, sigma, epsilon, bounds)
        paths[i] = path

    paths = pd.DataFrame(paths, index=path.index)

    return paths


def escaped_trajectories(paths, domain, before_bifurcation=True):
    """ Calculate the percentage of escaped trajectories at each parameter value.

        Parameters
        ----------
        paths : Dataframe
            Dataframe of sample paths.
        domain : callable
            Domain within which the path has not escaped given the `y` value.
        before_bifurcation : bool
            Whether to only calculate percentage of escaped trajectories before the 
            bifurcation point at `y = 0`.

        Returns
        -------
        Series
            The percentage of escaped trajectories at each parameter value.
    """
    num_escaped = []

    if before_bifurcation:
        paths = paths[paths.index < 0]

    runs = len(paths.columns)

    for y in paths.index:
        x = paths.loc[y]

        lower, upper = domain(y)
        escaped = np.vectorize(lambda x: np.nan if x < lower or x > upper else x)

        # sum the number of escaped trajectories
        num_escaped.append(sum(np.isnan(escaped(x.values))))
    
    percent_escaped = np.array(num_escaped) / runs * 100
 
    return pd.Series(percent_escaped, index=paths.index)


def distribution(paths, y):
    """ Get the distribution of x values at a given y value.

        Parameters
        ----------
        paths : Dataframe
            Dataframe of sample paths.
        y : float
            Given `y` value.

        Returns
        -------
        ndarray
            1D array of `x` values at the given `y` value.
    """
    return np.array(paths.iloc[np.argmin(np.abs(paths.index - y))])


def probability_density(pdf, domain, y):
    """ Get the probability densities at varying `x` values and a fixed `y` value.

        Parameters
        ----------
        pdf : callable
            Non-normalized probability density function `p(x,y)`.
        domain : callable
            Function to obtain `x` value bounds at a given `y`.
        y : float
            Given `y` value.

        Returns
        -------
        Series
            Normalised probability densities at varying `x` values.
    """
    bounds = domain(y)
    N, _ = quad(lambda x: pdf(x, y), *bounds)
    x = np.linspace(*bounds, 1000)

    return pd.Series(1/N*pdf(x, y), index=x)


def estimate_variance(paths, domain=None):
    """ Get an estimate of the variance from the sample paths for each `y` value.

        Parameters
        ----------
        paths : Dataframe
            Dataframe of sample paths.
        domain : callable, default=None
            Function to obtain `x` bounds at a given `y`.

        Returns
        -------
        Series
            Estimated variance of `x` values at each `y` value.
    """
    variance = []

    paths = paths[paths.index < 0]  # before bifurcation parameter

    for y in paths.index:
        x = np.array(paths.loc[y])

        if domain:
            lower, upper = domain(y)
            x = x[(x > lower) & (x < upper)]  # filter x values within domain

        variance.append(np.var(x))

    return pd.Series(variance, index=paths.index)


def analytical_variance(paths, domain, pdf):
    """ Calculate the analytical variance from the probability density function for each y value.

        Parameters
        ----------
        paths : Dataframe
            Dataframe of sample paths.
        domain : callable
            Function to obtain `x` bounds at a given `y` value.
        pdf : callable
            Non-normalised probabiility density function `p(x,y)`.
        
        Returns
        -------
        Series
            Analytical variance of `x` values at each `y` value.
    """
    variance = []

    paths = paths[paths.index < 0]

    for y in paths.index:
        lower, upper = domain(y)
        N, _ = quad(lambda x: pdf(x, y), lower, upper)  # normalisation constant
        mean, _ = quad(lambda x: x*1/N*pdf(x, y), lower, upper)
        var, _ = quad(lambda x: (x - mean)**2*1/N*pdf(x, y), lower, upper)

        variance.append(var)

    return pd.Series(variance, index=paths.index)


def probability_grid(xlim, ylim, pdf, domain):
    """ Obtain probability densities over (x, y) grid.

        Parameters
        ----------
        xlim : tuple of floats
            `x` value limits.
        ylim : tuple of floats
            `y` value limits.
        pdf : callable
            Probability density function `p(x,y)`.
        domain : callable
            Function to obtain `x` bounds at a given `y` value.

        Returns
        -------
        Y : ndarray
            `y` value at each point in a 2D grid.
        X : ndarray
            `x` value at each point in a 2D grid.
        Z : ndarray
            Probability density at each point in a 2D grid.


    """
    y = np.linspace(*ylim)
    x = np.linspace(*xlim)
    Y, X = np.meshgrid(y, x)
    N = lambda y: quad(lambda x: pdf(x, y), *domain(y))[0]
    Z = np.vectorize(lambda x, y: max(0, 1/N(y) * pdf(x, y)) if N(y) > 0 else 0)(X, Y)

    return Y, X, Z