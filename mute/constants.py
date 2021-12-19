##########################
##########################
###                    ###
###  MUTE              ###
###  William Woodley   ###
###  19 December 2021  ###
###                    ###
##########################
##########################

# Import packages

import os

import numpy as np

# Energies in [MeV]

E_BINS = np.logspace(1.9, 14, 122)[:-30]  # Bin edges
E_WIDTHS = np.diff(E_BINS)  # Bin widths
ENERGIES = np.sqrt(E_BINS[1:] * E_BINS[:-1])  # Bin centers

e_bins = E_BINS
e_widths = E_WIDTHS
energies = ENERGIES

# Slant depths in [km.w.e.] and angles in [degrees]
# These are the defaults used to construct the matrices
# The user can enter their own angles in the function calls, which will interpolate over the grids created by these angles

SLANT_DEPTHS = np.linspace(1, 12, 23)
ANGLES = np.degrees(np.arccos(1 / SLANT_DEPTHS))

slant_depths = SLANT_DEPTHS
angles = ANGLES

# Angles in [degrees] specifically for the calculation of surface flux matrices

ANGLES_FOR_S_FLUXES = np.linspace(0, 90, 20)

# Length values for file comparisons

len_ij = len(ENERGIES) * len(ANGLES_FOR_S_FLUXES)
len_iju = len(ENERGIES) * len(SLANT_DEPTHS) * len(ENERGIES)

# The rest mass of a muon in [MeV]

MU_MASS = 105.6583745

# Define user-set constants

_verbose = 2
_output = True
_directory = os.path.join(os.path.dirname(__file__), "data")
_lab = "Default"
_overburden = "flat"
_vertical_depth = 1
_medium = "rock"
_density = 2.65
_n_muon = 100000

# Setters and getters for global constants

# Verbose


def set_verbose(verbose):

    """
    Set the verbosity level.

    Parameters
    ----------
    verbose : int (default: 2)
        The verbosity level. Options:
        0 = Print no information.
        1 = Print the progress through the propagation loop if propagation muons or the loading of the Monte Carlo underground energies from files.
        2 = Print when calculations start and finish, and where data is read from and written to.
    """

    assert type(verbose) is int

    global _verbose

    _verbose = verbose


def get_verbose():

    """Return the set verbosity level."""

    return _verbose


# Whether or not the results are written to output files in directory


def set_output(output):

    """
    Set whether results are output to files or not. The default is True.

    Most functions also have their own optional output keyword argument to give more control over which output files are created.
    """

    assert type(output) is bool

    global _output

    _output = output


def get_output():

    """Return the set output setting."""

    return _output


# Directory for input and output


def set_directory(directory):

    """Set the working directory for files to be written to and read from. The default is \"data\"."""

    from .__init__ import download_file

    global _directory

    _directory = directory

    if not os.path.isfile(os.path.join(_directory, "data_20211219.txt")):

        download_file(
            "https://github.com/wjwoodley/mute/releases/download/0.1.0/data_20211219.zip",
            _directory,
        )


def get_directory():

    """Return the set working directory for files to be written to and read from."""

    return _directory


# Lab name used in output files


def set_lab(lab):

    """Set the name of the lab. This is used in the output file names as a way to differentiate different data files. The default is \"Default\"."""

    global _lab

    _lab = lab


def get_lab():

    """Return the set name of the lab."""

    return _lab


# Overburden


def set_overburden(overburden):

    """Set the overburden type. The default is  \"flat\" overburden."""

    assert overburden in [
        "flat",
        "mountain",
    ], 'overburden must be set to either "flat" or "mountain".'

    global _overburden

    _overburden = overburden


def get_overburden():

    """Return the set overburden type."""

    return _overburden


# Vertical depth in [km.w.e.]


def set_vertical_depth(vertical_depth):

    """Set the vertical depth, h. The default is 1 km.w.e."""

    global _vertical_depth
    global slant_depths
    global angles

    assert (
        vertical_depth >= 1
    ), "Vertical depths less than 1 km.w.e. are not currently supported."
    assert (
        vertical_depth <= 12
    ), "Vertical depths greater than 12 km.w.e. are not currenlty supported."

    _vertical_depth = vertical_depth

    # Use the set vertical depth to calculate new slant depths and zenith angles

    slant_depths = np.sort(
        np.concatenate(
            ([_vertical_depth], SLANT_DEPTHS[SLANT_DEPTHS > _vertical_depth])
        )
    )
    angles = np.degrees(np.arccos(_vertical_depth / slant_depths))


def get_vertical_depth():

    """Return the set vertical depth."""

    return _vertical_depth


# Medium


def set_medium(medium):

    """Set the medium for the muons to be propagated through. The default is standard rock (\"rock\")."""

    assert medium in [
        "rock",
        "water",
        "ice",
        "air",
    ], 'medium must be set to "rock", "water", "ice", or "air".'

    global _medium

    _medium = medium


def get_medium():

    """Return the set propagation medium."""

    return _medium


# Density in [gcm^-3]


def set_density(density):

    """Set the density of the propagation medium. The default is 2.65 g/cm^3 (the density of standard rock)."""

    assert density > 0, "Media density must be positive."

    global _density

    _density = density


def get_density():

    """Return the set density of the propagation medium."""

    return _density


# Number of muon


def set_n_muon(n_muon):

    """Set the number of muons to be propagated per surface energy-slant depth bin. The default is 10."""

    assert type(n_muon) is int, "Number of muons must be an integer."
    assert n_muon > 0, "Number of muons must be positive."

    global _n_muon

    _n_muon = n_muon


def get_n_muon():

    """Return the set number of muons to be propagated per surface energy-slant depth bin."""

    return _n_muon


# Check if required directories exist


def check_directory(directory, force=False):

    """Check whether a directory required to store output files exists or not. If it does not, ask the user if it should be created."""

    if not os.path.exists(directory):

        # If the directory does not exist, ask the user if they want to create it

        if not force:

            answer = input(
                directory
                + " does not currently exist. Would you like to create it (y/n)?: "
            )

        if force or answer.lower() == "y":

            os.makedirs(directory)

        else:

            print("Directory not created.")

            return

    else:

        return


# Check that the constants are set to correctly before running any calculation functions


def check_constants(force=False):

    """Check that the constants are set correctly before running any calculation functions."""

    # Check that the working directory the user has set exists

    check_directory(get_directory(), force=force)

    # Check that the overburden has been set to one of the available options

    if _overburden is None:

        raise ValueError('Overburden must be set to either "flat" or "mountain".')

    # Check that the vertical depth is set correctly, if needed

    elif _overburden == "flat":

        assert _vertical_depth is not None, "Initial vertical depth must be set."
        assert _vertical_depth > 0, "Initial vertical depth must be positive."

    # Check that the number of muons is set correctly

    if _n_muon <= 0:

        raise ValueError("Number of muons must be a positive integer.")


# Calculate weights


def calc_weights(
    file, units="kmwe", n_bins=[10, 10], max_slant_depth=None, force=False
):

    """
    Calculate the normalised weight matrix used for non-flat overburdens.

    Parameters
    ----------
    file : str
        The full path and name of the .txt file containing the geometry of the mountain. The first column should be the zenith angle in degrees; the second column should be the azimuthal angle in degrees; the third column should be the slant depth in units of units (see below).

    units : str, optional (default: "kmwe")
        The units of the slant depths in file. This must be one of ["m", "km", "mwe", "kmwe"].

    n_bins : array-like, optional (default: [10, 10])
        The number of bins in cos(theta) and slant depth, respectively

    max_slant_depth : float, optional (default: None)
        The maximum slant depth to use when creating the weights. Any data for slant depths above this value will be excluded. If None, all slant depths will be used.

    force : bool
        If True, this will force the creation of a directory to store output if one does not already exist.

    Returns
    -------
    cos_theta : NumPy array
        The centers of the cosine of the zenith angle bins of the weights histogram. The length is equal to n_bins[0].

    slant_depths : NumPy array
        The centers of the slant depth bins of the weights histogram. The length is equal to n_bins[1].

    weights : NumPy ndarray
        A two-dimensional array of shape (len(cos_theta), len(mu_slants)) containing the bin values of the normalised weights histogram.
    """

    # Check values

    assert (
        overburden == "mountain"
    ), "Weights can only be calculated for mountain overburdens."
    assert units_in in ["m", "km", "mwe", "kmwe"], "Units must be m, km, mwe, or kmwe."

    check_constants(force=force)

    # Global variables

    global mountain_zenith
    global mountain_azimuthal
    global mountain_slants

    global weights
    global cos_bin_width
    global slant_bin_width

    global slant_depths
    global angles

    # If the data is in [m] or [km], convert to [m.w.e.] or [km.w.e.]

    density_mult = 1

    if units == "m" or units == "km":
        density_mult = density / 0.997

    # If the data is in [m] or [m.w.e.], convert to [km.w.e.]

    scale = 1

    if units == "m" or units == "mwe":
        scale = 1e-3

    # Load the mountain geometry data from the file

    mountain_zenith_all = np.loadtxt(file)[:, 0]
    mountain_azimuthal_all = np.loadtxt(file)[:, 1]
    mountain_slants_all = np.loadtxt(file)[:, 2] * density_mult * scale

    mountain_zenith = np.unique(mountain_zenith_all)
    mountain_azimuthal = np.unique(mountain_azimuthal_all)
    mountain_slants = np.reshape(
        mountain_slants_all, (len(mountain_zenith), len(mountain_azimuthal))
    )

    # Define bins edges from the data

    if max_slant is None:

        slant_bins = np.sort(np.histogram(mountain_slants_all, bins=n_bins[1])[1])
        cos_bins = np.linspace(
            np.cos(np.radians(max(mountain_zenith_all))), 1, n_bins[0] + 1
        )

    else:

        sel = np.where(np.array(mountain_slants_all) < max_slant_in)

        slant_bins = np.sort(np.histogram(mountain_slants_all[sel], bins=n_bins[1])[1])
        cos_bins = np.linspace(
            np.cos(np.radians(max(mountain_zenith_all[sel]))), 1, n_bins[0] + 1
        )

    # Calculate the weights

    weights = np.histogram2d(
        np.cos(np.radians(mountain_zenith_all)),
        mountain_slants_all,
        bins=[cos_bins, slant_bins],
        density=True,
    )[0]

    # Calculate the bin widths

    slant_bin_width = np.diff(slant_bins)[0]
    cos_bin_width = np.diff(cos_bins)[0]

    # Calculate the bin centers

    slant_depths = (slant_bins[1:] + slant_bins[:-1]) / 2
    angles = np.degrees(np.arccos((cos_bins[1:] + cos_bins[:-1]) / 2))

    return cos_theta, slant_depths, weights


# Clear all variables


def clear():

    """Reset all of the values set or calculated in MUTE to their default values."""

    # Import packages

    import gc

    # Global constants and variables

    global _verbose
    global _output
    global _directory
    global _lab
    global _overburden
    global _vertical_depth
    global _medium
    global _density
    global _n_muon

    _verbose = 2
    _output = True
    _directory = os.path.join(os.path.dirname(__file__), "data")
    _lab = "Default"
    _overburden = "flat"
    _vertical_depth = 1
    _medium = "rock"
    _density = 2.65
    _n_muon = 100000

    # Global weight variables
    # Set the variables to None first, in case they do not already exist

    global mountain_zenith
    global mountain_azimuthal
    global mountain_slants
    global weights
    global cos_bin_width
    global slant_bin_width

    mountain_zenith = None
    mountain_azimuthal = None
    mountain_slants = None
    weights = None
    cos_bin_width = None
    slant_bin_width = None

    del mountain_zenith
    del mountain_azimuthal
    del mountain_slants
    del weights
    del cos_bin_width
    del slant_bin_width

    # Energies, slant depths, zenith angles, and the muon mass

    global E_BINS
    global E_WIDTHS
    global ENERGIES
    global e_bins
    global e_widths
    global energies
    global SLANT_DEPTHS
    global ANGLES
    global slant_depths
    global angles
    global ANGLES_FOR_S_FLUXES
    global len_ij
    global len_iju
    global MU_MASS

    E_BINS = np.logspace(1.9, 14, 122)[:-30]  # [MeV]
    E_WIDTHS = np.diff(E_BINS)  # [MeV]
    ENERGIES = np.sqrt(E_BINS[1:] * E_BINS[:-1])  # [MeV]
    e_bins = E_BINS  # [MeV]
    e_widths = E_WIDTHS  # [MeV]
    energies = ENERGIES  # [MeV]
    SLANT_DEPTHS = np.linspace(1, 12, 23)  # [km.w.e.]
    ANGLES = np.degrees(np.arccos(1 / SLANT_DEPTHS))  # [degrees]
    slant_depths = SLANT_DEPTHS  # [km.w.e.]
    angles = ANGLES  # [degrees]
    ANGLES_FOR_S_FLUXES = np.linspace(0, 90, 20)  # [degrees]
    len_ij = len(ENERGIES) * len(ANGLES_FOR_S_FLUXES)
    len_iju = len(ENERGIES) * len(SLANT_DEPTHS) * len(ENERGIES)
    MU_MASS = 105.6583745  # [MeV]

    gc.collect()
