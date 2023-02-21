import numpy as np
from matplotlib.animation import ArtistAnimation
import matplotlib.pyplot as plt

from spatial_indicators import *
from spatial_models import *


def plot_snapshots(snapshots, param_values, param_label=None, vmin=None, vmax=None, precision=1):
    """ Plot spatial snapshots at given parameter levels.

        Parameters
        ----------
        snapshots : Series
            Spatial snapshots at different parameter levels.
        param_values : list of floats
            Control parameter values at which to plot snapshots.
        param_label : str, optional
            Parameter name.
        vmin : float, optional
            Minimum spatial value to be mapped.
        vmax : float, optional
            Maximum spatial value to be mapped.
        precision : int, default=1
            Number of decimal points to display of parameter value.
    """
    fig, axs = plt.subplots(ncols = len(param_values))

    for i in range(len(param_values)):
        # get snapshot at parameter value
        snapshot = snapshots.iloc[np.argmin(np.abs(snapshots.index - param_values[i]))]

        # set title of axes
        title = f"{param_values[i]:.{precision}f}"
        if param_label:
            title = f"{param_label} = {title}"

        axs[i].imshow(snapshot, vmin=vmin, vmax=vmax, cmap="gray")
        axs[i].set_title(f"{param_label} = {param_values[i]:.{precision}f}")
        axs[i].axis("off")

    fig.set_size_inches(2*len(param_values), 2)


def plot_distributions(snapshots, param_values, param_label=None, precision=1):
    """ Plot distributions of values in snapshots at given parameter levels.

        Parameters
        ----------
        snapshots : Series
            Spatial snapshots at different parameter levels.
        param_values : list of floats
            Control parameter values at which to plot snapshot distributions.
        param_label : str, optional
            Parameter name.
        precision : int, default=1
            Number of decimal places to display of parameter value.
    """
    fig, axs = plt.subplots(ncols = len(param_values))

    for i in range(len(param_values)):
        # get snapshot at parameter value
        snapshot = snapshots.iloc[np.argmin(np.abs(snapshots.index - param_values[i]))]

        # set title of axes
        title = f"{param_values[i]:.{precision}f}"
        if param_label:
            title = f"{param_label} = {title}"

        axs[i].hist(snapshot.flatten(), density=True)
        axs[i].set_title(f"{param_label} = {param_values[i]:.{precision}f}")

    axs[0].set_ylabel("Density")

    fig.tight_layout()
    fig.set_size_inches(2*len(param_values), 2)


def animate_snapshots(snapshots, vmin=None, vmax=None, step=1):
    """ Animate snapshots as the control parameter changes.

        Parameters
        ----------
        snapshots : Series
            Spatial snapshots at different parameter levels.
        vmin : float, optional
            Minimum spatial value to be mapped.
        vmax : float, optional
            Maximum spatial value to be mapped.
        step : int, default=1
            Number of snapshots between each frame to animate.
    """
    fig, ax = plt.subplots()

    ims = []
    for i in range(0, len(snapshots), step):
        im = ax.imshow(snapshots.iloc[i], vmin=vmin, vmax=vmax, cmap="gray", animated=True)
        ims.append([im])

    plt.close()

    return ArtistAnimation(fig, ims, interval=100, blit=True, repeat=False)


def plot_spatial_indicator(ax, snapshots, indicator, param=None, n=1):
    """ Plot spatial indicator against control parameter level.

        Parameters
        ----------
        ax : obj
            Axes object to plot on.
        snapshots : Series
            Spatial snapshots at each parameter level.
        indicator : str
            Spatial indicator to calculate (`"mean"`, `"variance"`, `"skewness"` or `"correlation"`).
        param : str, optional
            Control parameter name.
        n : int, default=1
            Number of control parameter increments per spatial indicator calculation.

        Returns
        -------
        line : obj
            Line object.
    """
    indicator_functions = {
        "mean": spatial_mean,
        "variance": spatial_variance,
        "skewness": spatial_skewness,
        "correlation": spatial_correlation
    }

    indicator_values = [indicator_functions[indicator](snapshot) for snapshot in snapshots.iloc[::n]]

    ax.set_ylabel(f"Spatial {indicator.capitalize()}")
    if param:
        ax.set_xlabel(param)

    line, = ax.plot(snapshots.iloc[::n].index, indicator_values)

    return line


def plot_spatial_indicator_grid(axs, snapshots, indicators, shape, param, n, legend_loc):
    """ Plot spatial indicators against control parameter level on sets of axes.

        Parameters
        ----------
        axs : obj
            Axes object.
        snapshots : Series
            Spatial snapshots at each parameter level.
        indicators : list of strings
            Spatial indicators to calculate (`"mean"`, `"variance"`, `"skewness"` or `"correlation"`).
        shape : tuple of ints
            Number of rows and columns in grid.
        param : str
            Control parameter name.
        n : int
            Number of control parameter increments per spatial indicator calculation.
        legend_loc : int
            Axes number to place legend on.

        Returns
        -------
        line : obj
            Line object on first set of axes.
    """
    idx = 0
    nrows, ncols = shape

    for i in range(nrows):
        for j in range(ncols):
            # get axis
            if nrows == 1:
                ax = axs[j]
            elif ncols == 1:
                ax = axs[i]
            else:
                ax = axs[i][j]

            # plot line
            if legend_loc == (i, j):
                line = plot_spatial_indicator(ax, snapshots, indicators[idx], param, n)
            elif idx < len(indicators):
                plot_spatial_indicator(ax, snapshots, indicators[idx], param, n)

            idx += 1

    return line


def spatial_indicator_grid(snapshots, param=None, indicators=None, shape=None, levels=None, level_name=None, n=1, legend_loc=None, return_fig=False):
    """ Plot spatial indicators against control parameter level. If there are multiple levels of snapshots, plot multiple 
        lines on each set of axes.

        Parameters
        ----------
        snapshots : Series or list of Series
            Spatial snapshots at each parameter level.
        param : str, optional
            Control parameter name.
        indicators : list of strings, optional
            Spatial indicators to calculate (`"mean"`, `"variance"`, `"skewness"` or `"correlation"`). 
            Default is `["mean", "variance", "skewness", "correlation"]`.
        shape : tuple of ints, optional
            Number of rows and columns in grid. Default is `(2, 2)`.
        levels : list, optional
            Different levels of snapshots series.
        level_name : str, optional
            Level name.
        n : int, default=1
            Number of control parameter increments per spatial indicator calculation.
        legend_loc : tuple of ints, optional
            Location of axes to place legend on.
        return_fig : bool, default=False
            Whether to return the figure object.

        Returns
        -------
        fig : obj
            Figure object if `return_fig` is `True`, otherwise `None`.
    """
    # default indicators
    if indicators is None:
        indicators = ["mean", "variance", "skewness", "correlation"]
    
    # default shape
    if shape is None:
        shape = (2, 2)
    
    nrows, ncols = shape

    # default legend location
    if legend_loc is None:
        legend_loc = (0, 0)

    fig, axs = plt.subplots(nrows=nrows, ncols=ncols)

    # if there is more than one level
    if levels is not None:
        assert len(snapshots) == len(levels), "Number of snapshots series does not match number of levels."

        # plot each level and save the lines
        lines = []
        for snapshots_level in snapshots:
            lines.append(plot_spatial_indicator_grid(axs, snapshots_level, indicators, shape, param, n, legend_loc))
        
        # set up legend
        if level_name:
            labels = [f"{level_name} = {level}" for level in levels]
        else:
            labels = levels

        i, j = legend_loc
        axs[i][j].legend(lines, labels)

    else:
        plot_spatial_indicator_grid(axs, snapshots, indicators, shape, param, n, legend_loc)

    fig.tight_layout()

    if return_fig:
        return fig


def plot_power_spectra(snapshots, param_values, param_label=None, precision=1):
    """ Plot power spectra at given parameter levels.

        Parameters
        ----------
        snapshots : Series
            Spatial snapshots at different parameter levels.
        param_values : list of floats
            Control parameter values at which to plot power spectrums.
        param_label : str, optional
            Parameter name.
        precision : int, default=1
            Number of decimal points to display of parameter value.
    """
    fig, axs = plt.subplots(ncols = len(param_values))

    for i in range(len(param_values)):
        # get snapshot at parameter value
        snapshot = snapshots.iloc[np.argmin(np.abs(snapshots.index - param_values[i]))]

        # set title of axes
        title = f"{param_values[i]:.{precision}f}"
        if param_label:
            title = f"{param_label} = {title}"

        axs[i].imshow(power_spectrum(snapshot), cmap="gray")
        axs[i].set_title(f"{param_label} = {param_values[i]:.{precision}f}")
        axs[i].axis("off")


def plot_radial_spectra(snapshots, dx, param_values, param_label=None):
    """ Plot radial spectra at given parameter levels on single set of axes.

        Parameters
        ----------
        snapshots : Series
            Spatial snapshots at different parameter levels.
        dx : float
            Spatial step size.
        param_values : list of floats
            Control parameter values at which to plot radial spectra.
        param_label : str
            Parameter name.
    """
    fig, ax = plt.subplots()

    legend_labels = []

    for i in range(len(param_values)):
        # get snapshot at parameter value
        snapshot = snapshots.iloc[np.argmin(np.abs(snapshots.index - param_values[i]))]

        # plot the radial spectrum
        ax.plot(radial_spectrum(power_spectrum(snapshot), dx))

        # set legend label
        if param_label:
            legend_labels.append(f"{param_label} = {param_values[i]}")
        else:
            legend_labels.append(f"{param_values[i]}")

    ax.set_xlabel("Wavenumber")
    ax.set_ylabel("Radial Spectrum")
    fig.legend(legend_labels)
            

def plot_equilibria(ax, model, param, param_name, *args):
    """ Plot the homogeneous equilibrium of a model on a set of axes.

        Parameters
        ----------
        ax : obj
            Axes object.
        model : obj
            Model object, must have method `equilibria()`.
        param : iterable
            Range of parameter values to plot equilibrium.
        param_name : str
            Name of parameter, must match attribute in model.
    """
    eq = []
    for p in param:
        setattr(model, param_name, p)
        eq.append(model.equilibria(*args)[0])
    ax.plot(param, eq, '--')


def plot_max_eigenvalues(ax, model, param, param_name, k=None):
    """ Plot maximum real parts of eigenvalues over a parameter range on a set of axes.

        Parameters
        ----------
        ax : obj
            Axes object.
        model : obj
            Model object.
        param : iterable
            Range of parameter values to plot maximum eigenvalues.
        param_name : str
            Name of parameter, must match attribute in model.
        k : iterable, optional
            Range of wavenumbers to plot maximum eigenvalues.
    """
    w_max = []

    if k is not None:
        for p in param:
            setattr(model, param_name, p)
            w_max.append(max_eigenvalue(model, k))
            ax.plot(w_max[-1])

        ax.plot(k, [0] * len(k), '--')

        ax.set_xlabel("k")
        ax.legend([f"{param_name} = {p:.2f}" for p in param])

    else:
        for p in param:
            setattr(model, param_name, p)
            w_max.append(max_eigenvalue(model))

        ax.plot(param, w_max)

        ax.set_xlabel(param_name)
        ax.plot(param, [0] * len(param), '--')

    ax.set_ylabel(r"$Re(\sigma_{max})$")


def plot_time_simulations(model, initial, time_steps):
    """ Plot simulations of a model from different initial conditions.

        Parameter
        ---------
        model : obj
            Model object.
        initial : iterable
            Initial vegetation conditions.
        time_steps : int
            Number of time steps to simulate in each simulation.
    """
    fig, ax = plt.subplots()

    equilibria_values = model.equilibria()

    for init in initial:
        # set up initial values
        initial_values = tuple(value if j > 0 else init for j, value in enumerate(equilibria_values))
        initial_grids = tuple(value*np.ones((model.size, model.size)) for value in initial_values)

        # simulate snapshots
        snapshots = simulate_time_steps(model, time_steps, *initial_grids)

        # calculate and plot mean
        mean = [spatial_mean(snapshot) for snapshot in snapshots]
        ax.plot(snapshots.index, mean, color='k')

    # plot equilibrium
    ax.plot(snapshots.index, [equilibria_values[0]]*len(snapshots.index), linestyle='--', color='r')

    # axis labels
    ax.set_ylabel("Spatial Mean")
    ax.set_xlabel("t")