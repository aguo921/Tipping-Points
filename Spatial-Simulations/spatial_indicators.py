import numpy as np
from scipy.stats import skew
import pandas as pd


def spatial_mean(snapshot):
    """ Calculate mean of values within a 2D spatial snapshot.

        Parameters
        ----------
        snapshot : ndarray
            2D array of numeric values.
        
        Returns
        -------
        float
            Spatial mean.    
    """
    return np.mean(snapshot.flatten())


def spatial_variance(snapshot):
    """ Calculate variance of values within a 2D spatial snapshot.

        Parameters
        ----------
        snapshot : ndarray
            2D array of numeric values.
        
        Returns
        -------
        float
            Spatial variance.
    """
    return np.var(snapshot.flatten())


def spatial_skewness(snapshot):
    """ Calculate skewness of values within a 2D spatial snapshot.

        Parameters
        ----------
        snapshot : ndarray
            2D array of numeric values.

        Returns
        -------
        float
            Spatial skewness.
    """
    return skew(snapshot.flatten())


def spatial_correlation(snapshot):
    """ Calculate autocorrelation of pairs cells within a 2D spatial snapshot using global Moran's coefficient.

        Parameters
        ----------
        snapshot : ndarray
            2D array of numeric values.

        Returns
        -------
        float
            Spatial correlation.
    """
    mean = np.mean(snapshot.flatten())  # spatial mean

    sum_squares = np.sum((snapshot - mean)**2)

    # get nearest neighbours
    top = np.roll(snapshot, (0, 1), (1, 0))
    bottom = np.roll(snapshot, (0, -1), (1, 0))
    left = np.roll(snapshot, (-1, 0), (1, 0))
    right = np.roll(snapshot, (1, 0), (1, 0))

    # calculate Moran's coefficient
    C = np.sum((snapshot - mean)*(top + bottom + left + right - 4*mean))
    if sum_squares == 0:
        return np.NaN
    else:
        return C/(4*sum_squares)


def bootstrapping(snapshot):
    """ Shuffle elements within a 2D spatial snapshot.

        Parameters
        ----------
        snapshot : ndarray
            2D array of numeric values.
        
        Returns
        -------
        ndarray
            Shuffled 2D array.
    """
    flattened_snapshot = snapshot.flatten()
    np.random.shuffle(flattened_snapshot)

    return flattened_snapshot.reshape(snapshot.shape)


def coarse_graining(snapshot, s):
    """ Create a coarse grained version of a 2D spatial snapshot with submatrices of size sxs.

        Parameters
        ----------
        snapshot : ndarray
            2D array of values.
        s : int
            Size of submatrices.
        
        Returns
        -------
        coarse_grained : ndarray
            Smaller coarse grained 2D array of mean values.
    """
    # number of rows and columns within original snapshot
    nrows, ncols = snapshot.shape

    # initialise coarse grained matrix
    coarse_grained = np.zeros((int(np.ceil(nrows / s)), int(np.ceil(ncols / s))))

    # loop through each submatrix of size sxs
    for i in range(0, nrows, s):
        for j in range(0, ncols, s):
            submatrix = snapshot[i:min(i+s,nrows), j:min(j+s,ncols)]

            # add mean of submatrix to coarse grained matrix
            coarse_grained[i//s, j//s] = np.mean(submatrix)

    return coarse_grained


def power_spectrum(snapshot):
    """ Calculate scaled power spectrum of spatial snapshot.

        Parameters
        ----------
        snapshot : ndarray
            2D array of spatial values.

        Returns
        -------
        Scaled power spectra of spatial snapshot.
    """
    m, n = snapshot.shape
    fourier_transform = np.fft.fft2(snapshot, norm="ortho")

    return np.abs(fourier_transform[1:m//2+1,1:n//2+1])**2 / spatial_variance(snapshot)


def wavenumbers(snapshot, dx=1):
    """ Return wavenumbers along each dimension of spatial snapshot.

        Parameters
        ----------
        snapshot : ndarray
            2D array of spatial values.
        dx : float, default=1
            Spatial step size.

        Returns
        -------
        p : ndarray
            1D array of positive wavenumbers in x direction.
        q : ndarray
            1D array of positive wavenumbers in y direction.
    """
    m, n = snapshot.shape

    _, *p = np.fft.rfftfreq(m, dx)
    _, *q = np.fft.rfftfreq(n, dx)

    return p, q


def radial_spectrum(snapshot, dx=1):
    """ Calculate radial spectrum of spatial snapshot.

        Parameters
        ----------
        snapshot : ndarray
            2D array of spatial values.
        dx : float, default=1
            Spatial step size.

        Returns
        -------
        Series
            Sum of power spectra at different distances from the origin.
    """
    # get the power spectrum of the spatial snapshot
    p_spectrum = power_spectrum(snapshot)

    # find wavenumbers
    p, q = wavenumbers(snapshot, dx)
    P, Q = np.meshgrid(p, q, indexing="ij")

    # calculate radial distance for each pair of wavenumbers
    r_distance = np.sqrt(P**2 + Q**2)

    # get range of radial distance
    r, dr = np.linspace(0, max(p), 10, retstep=True)

    r_spectrum = []

    for r_step in r:
        # find pairs of wavenumbers at the current radial distance
        values = p_spectrum[np.where((r_distance >= r_step) & (r_distance < r_step + dr))]

        # sum the power spectra
        r_spectrum.append(np.sum(values))

    return pd.Series(r_spectrum, r)