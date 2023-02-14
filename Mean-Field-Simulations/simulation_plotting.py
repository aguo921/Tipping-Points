from simulation_functions import *
from matplotlib import cm


def plot_sample_path(ax, path):
    """ Plot a sample path on axes.

        Parameters
        ----------
        ax : obj
            Axes object.
        path : Series
            Sample path of `x` values as `y` slowly increases.
    """
    ax.plot(path)
    ax.set_xlabel("y")
    ax.set_ylabel("x")


def plot_ensemble(ax, paths, equilibria=[], labels=True):
    """ Plot ensemble of sample paths and equilibria on axes.

        Parameters
        ----------
        ax : obj
            Axes object.
        paths : Dataframe
            Dataframe of sample paths.
        equilibria : list of dicts
            List of dictionarys of the `y` range, function and stability of each equilibrium.
    """
    for _, x in paths.items():
        ax.plot(x, alpha=0.1, color="black")

    for equilibrium in equilibria:
        y = np.linspace(*equilibrium["range"])
        x = np.vectorize(equilibrium["func"])(y) # [equilibrium["func"](yt) for yt in y]

        ax.plot(y, x, 'r-' if equilibrium["stable"] else 'r--')

    if labels:
        ax.set_xlabel("y")
        ax.set_ylabel("x")


def plot_escaped(ax, percent_escaped, labels=True):
    """ Plot percentage of escaped trajectories at each `y` value on axes.

        Parameters
        ----------
        ax : obj
            Axes object.
        percent_escaped : Series
            Percentage of escaped trajectories at each `y` value.
        labels : bool, optional
            Whether to display axis labels, default is True.
    """
    ax.plot(percent_escaped)
    ax.set_ylim((0, 100))

    if labels:
        ax.set_xlabel("y")
        ax.set_ylabel("% escaped")


def plot_distribution(ax, paths, y, pdf, xlim, domain, labels=True):
    """ Plot histogram of `x` values from all the runs at a given `y` value on axes.

        Parameters
        ----------
        ax : obj
            Axes object.
        paths : Dataframe
            Dataframe of sample paths.
        y : float
            `y` value at which to find `x` values.
        pdf : callable
            Probability density function `p(x,y)`.
        xlim : tuple of floats
            Limits of x-axis.
        domain : callable
            Function to obtain `x` bounds at a given `y` value.
        labels : bool, optional
            Whether to display axis labels, default is True.
    """
    x = distribution(paths, y)
    p = probability_density(pdf, domain, y)

    ax.hist(x, bins=10, density=True)
    ax.plot(p)
    ax.set_xlim(xlim)

    if labels:
        ax.set_xlabel("x(y = {y})")
        ax.set_ylabel("Probability density")


def plot_probability_grid(ax, pxy, paths=None, labels=True, ylim=None):
    """ Plot contour plot of the probability density function `p(x,y)` on (x, y) grid.

        Parameters
        ----------
        ax : obj
            Axes object.
        pxy : tuple of ndarrays
            Tuple containing the x values, y values and probability density values on a 2D grid respectively.
        paths : Dataframe, optional
            Dataframe of sample paths.
        labels : bool, optional
            Whether to display axis labels, default is True.
        yliim : tuple of floats, optional
            Limits of y-axis.

        Returns
        -------
        contour : obj
            Contour object.
    """
    contour = ax.contourf(*pxy, levels=np.arange(0, 6, 0.1), cmap=cm.coolwarm)

    if paths is not None:
        plot_ensemble(ax, paths[paths.index < 0].iloc[:,:30], labels=labels)

    if labels:
        ax.set_xlabel("y")
        ax.set_ylabel("x")

    if ylim is not None:
        ax.set_ylim(*ylim)

    return contour


def plot_variance(ax, estimate_var, analytical_var=None, labels=True):
    """ Plot estimate variance of sample paths and analytical variance (if applicable) on axes.

        Parameters
        ----------
        ax : obj
            Axes object.
        estimate_var : Series
            Estimate variance at each `y` value.
        analytical_var : Series, optional
            Analytical variance at each `y` value.
        labels : bool, optional
            Whether to display axis labels, default is True.
    """
    ax.plot(estimate_var)

    if analytical_var is not None:
        ax.plot(analytical_var)

    if labels:
        ax.set_xlabel("y")
        ax.set_ylabel("Variance")
