#########################
#########################
###                   ###
###  MUTE             ###
###  William Woodley  ###
###  15 July 2022     ###
###                   ###
#########################
#########################

# Import packages

import os
from collections import namedtuple

import numpy as np

# Energies in [MeV]

_E_BINS = np.logspace(1.9, 14, 122)[:-30]  # Bin edges
_E_WIDTHS = np.diff(_E_BINS)  # Bin widths
ENERGIES = np.sqrt(_E_BINS[1:] * _E_BINS[:-1])  # Bin centers

# Slant depths in [km.w.e.] and angles in [degrees]
# These are the defaults used to construct the matrices
# The user can enter their own angles in the function calls, which will interpolate over the grids created by these angles

_X_MIN = 0.5
_X_MAX = 14

_SLANT_DEPTHS = np.linspace(_X_MIN, _X_MAX, int(2 * (_X_MAX - _X_MIN) + 1))
_ANGLES = np.degrees(np.arccos(_X_MIN / _SLANT_DEPTHS))

slant_depths = _SLANT_DEPTHS
angles = _ANGLES

# Angles in [degrees] specifically for the calculation of surface flux matrices

ANGLES_FOR_S_FLUXES = np.linspace(0, 90, 20)

# Other constants
# The rest mass of a muon in [MeV]
# Months of the year

_MU_MASS = 105.6583745

MONTHS = [
    "January",
    "February",
    "March",
    "April",
    "May",
    "June",
    "July",
    "August",
    "September",
    "October",
    "November",
    "December",
]
MONTHS_SNAMES = [
    "Jan.",
    "Feb.",
    "Mar.",
    "Apr.",
    "May.",
    "Jun.",
    "Jul.",
    "Aug.",
    "Sep.",
    "Oct.",
    "Nov.",
    "Dec.",
]

# Define default user-set global constants

_verbose = 2
_output = True
_directory = os.path.join(os.path.dirname(__file__), "data")
_lab = "Default"
_overburden = "flat"
_vertical_depth = _X_MIN
_medium = "rock"
_density = 2.65
_n_muon = 100000

# Keep track of which survival probability tensor is loaded

_survival_probability_tensor_configuration = {}
_current_survival_probability_tensor = None

# Keep track of whether or not a mountain profile has been loaded

_mountain_loaded = False

# Use extrapolation for depths lower than 0.5 km.w.e.
# This limitation comes from PROPOSAL, not MUTE

shallow_extrapolation = False

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
        1 = Print the progress through the calculation of surface fluxes if calculating surface fluxes, the propagation loop if propagating muons, or the loading of the Monte Carlo underground energies from files if loading underground energies.
        2 = Print when calculations start and finish, and where data is read from and written to.
    """

    assert isinstance(verbose, int)

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

    assert isinstance(output, bool)

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

    if not os.path.isfile(os.path.join(_directory, "data_20220715.txt")):

        download_file(
            "https://github.com/wjwoodley/mute/releases/download/0.1.0/data_20220715.zip",
            _directory,
        )


def get_directory():

    """Return the set working directory for files to be written to and read from."""

    return _directory


# Lab name used in output files


def set_lab(lab):

    """Set the name of the lab. This is used in the output file names in the \"underground\" data directory as a way to differentiate different data files. The default is \"Default\"."""

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

    """Set the vertical depth, h. The default is 0.5 km.w.e."""

    global _vertical_depth
    global slant_depths
    global angles

    _vertical_depth = vertical_depth

    # Use the set vertical depth to calculate new slant depths and zenith angles
    # Only do this for flat overburdens
    # For mountains, the slant depths and angles will be calculated in load_mountain()

    if vertical_depth < _SLANT_DEPTHS[0] and not shallow_extrapolation:

        raise Exception(
            "The minimum default available slant depth is 0.5 km.w.e. Set constants.shallow_extrapolation to True to enable calculations for depths lower than 0.5 km.w.e. (not recommended)."
        )

    if get_overburden() == "flat":

        slant_depths = np.sort(
            np.concatenate(
                ([_vertical_depth], _SLANT_DEPTHS[_SLANT_DEPTHS > _vertical_depth])
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

    """Set the density of the propagation medium. The default is 2.65 gcm^-3 (the density of standard rock)."""

    assert density > 0, "Media density must be positive."

    global _density

    _density = density


def get_density():

    """Return the set density of the propagation medium."""

    return _density


# Number of muon


def set_n_muon(n_muon):

    """Set the number of muons to be propagated per surface energy-slant depth bin. The default is 100000, as the provided default survival probability tensors in the data directory are given for 100000 muons."""

    assert isinstance(n_muon, int), "Number of muons must be an integer."
    assert n_muon > 0, "Number of muons must be positive."

    global _n_muon

    _n_muon = n_muon


def get_n_muon():

    """Return the set number of muons to be propagated per surface energy-slant depth bin."""

    return _n_muon


# Check if required directories exist


def _check_directory(directory, force=False):

    """This function checks whether a directory required to store output files exists or not. If it does not, ask the user if it should be created."""

    if not os.path.exists(directory):

        # If the directory does not exist, ask the user if they want to create it

        if not force:

            answer = input(
                f"{directory} does not currently exist. Would you like to create it (y/n)?: "
            )

        if force or answer.lower() == "y":

            os.makedirs(directory)

        else:

            print("Directory not created.")

            return

    else:

        return


# Check that the constants are set to correctly before running any calculation functions


def _check_constants(force=False):

    """This function checks that the constants are set correctly before running any calculation functions."""

    # Check that the working directory the user has set exists

    _check_directory(get_directory(), force=force)

    # Check that the overburden has been set to one of the available options

    if _overburden is None:

        raise ValueError('Overburden must be set to either "flat" or "mountain".')

    # Check that the vertical depth is set correctly, if needed

    elif get_overburden() == "flat":

        assert get_vertical_depth() is not None, "Initial vertical depth must be set."
        assert get_vertical_depth() > 0, "Initial vertical depth must be positive."

        assert len(slant_depths) == len(angles) and np.allclose(
            slant_depths, get_vertical_depth() / np.cos(np.radians(angles))
        ), "slant_depths must correspond to angles. Do not assign constants.slant_depths or constants.angles directly. Instead, use constants.set_vertical_depth() in combination with the angles and depths parameters in individual functions. Run constants.clear() to reset the values of slant_depths and / or angles."

    # Check that a mountain profile has been loaded, if needed

    elif get_overburden() == "mountain":

        if not _mountain_loaded:

            raise Exception(
                "The mountain profile must be loaded first by passing a txt file to constants.load_mountain()."
            )

    else:

        raise NotImplementedError(
            'Overburdens of type {0} are not available. The only options are "flat" and "mountain".'.format(
                constants.get_overburden()
            )
        )

    # Check that the number of muons is set correctly

    if _n_muon <= 0:

        raise ValueError("Number of muons must be a positive integer.")


# Load mountain data


def load_mountain(file_name, units="kmwe", density=get_density(), max_slant_depth=14):

    """
    Load mountain data from a profile file.

    The first column should be the zenith angle in [degrees]; the second column should be the azimuthal angle in [degrees]; the third column should be the slant depth in units of units (see below). This function makes available the variables listed under "Sets" below.

    Parameters
    ----------
    file_name : str
        The full path and name of the .txt file containing the profile information of the mountain.

    units : str, optional (default: "kmwe")
        The units of the slant depths in file_name. This must be one of ["m", "km", "mwe", "kmwe"].

    density : float, optional (default: taken from constants.get_density())
        The density in [g cm^-3] of the rock for which the depths are specified in the file.

    max_slant_depth : float, optional (default: 14)
        The maximum slant depth in [km.w.e.] to take from the file. Any data for slant depths above this value will be set to 0 km.w.e. and will be ignored throughout calculations. The default is 14 km.w.e., consistent with the maximum slant depth in constants._SLANT_DEPTHS being 14 km.w.e.

    Sets
    ----
    constants.mountain.zenith : NumPy ndarray
        Unique zenith angles in [degrees] sorted from smallest to largest.

    constants.mountain.azimuthal : NumPy ndarray
        Unique azimuthal angles in [degrees] sorted from smallest to largest.

    constants.mountain.slant_depths : NumPy ndarray
        Unique slant depths in [km.w.e.] in a matrix of shape (len(constants.mountain_zenith), len(constants.mountain_azimuthal)).
    """

    # Check values

    assert (
        get_overburden() == "mountain"
    ), 'The overburden type must be set to "mountain".'
    assert units in [
        "m",
        "km",
        "mwe",
        "kmwe",
    ], 'Units must be "m", "km", "mwe", or "kmwe".'

    # Global variables

    global _mountain_loaded
    global _mountain_zenith_all
    global _mountain_azimuthal_all
    global _mountain_slant_depths_all

    global mountain

    # Create a named tuple to hold immutable results

    Mountain = namedtuple("Mountain", ("zenith", "azimuthal", "slant_depths"))

    # If the data is in [m] or [km], convert to [m.w.e.] or [km.w.e.]

    density_mult = 1

    if units == "m" or units == "km":
        density_mult = density / 0.997

    # If the data is in [m] or [m.w.e.], convert to [km.w.e.]

    scale = 1

    if units == "m" or units == "mwe":
        scale = 1e-3

    # Load the mountain profile data from the file

    _mountain_zenith_all = np.loadtxt(file_name)[:, 0]
    _mountain_azimuthal_all = np.loadtxt(file_name)[:, 1]
    _mountain_slant_depths_all = np.loadtxt(file_name)[:, 2] * density_mult * scale

    # Extract the unique angles and slant depths

    mountain_zenith = np.unique(_mountain_zenith_all)
    mountain_azimuthal = np.unique(_mountain_azimuthal_all)
    mountain_slant_depths = np.reshape(
        np.nan_to_num(_mountain_slant_depths_all),
        (len(mountain_zenith), len(mountain_azimuthal)),
    )

    # Remove the slant depths above the maximum slant depth threshold by setting them to 0

    mountain_slant_depths[mountain_slant_depths > max_slant_depth] = 0

    # Assign the elements of the namedtuple

    mountain = Mountain(mountain_zenith, mountain_azimuthal, mountain_slant_depths)

    # If the mountain profile loads successfully, set _mountain_loaded to True

    _mountain_loaded = True


# Clear all variables


def clear():

    """Reset all of the values set or calculated in MUTE to their default values."""

    # Import packages

    import gc

    # Energies, slant depths, and zenith angles

    global _E_BINS
    global _E_WIDTHS
    global ENERGIES
    global _X_MIN
    global _X_MAX
    global _SLANT_DEPTHS
    global _ANGLES
    global slant_depths
    global angles
    global ANGLES_FOR_S_FLUXES

    _E_BINS = np.logspace(1.9, 14, 122)[:-30]  # [MeV]
    _E_WIDTHS = np.diff(_E_BINS)  # [MeV]
    ENERGIES = np.sqrt(_E_BINS[1:] * _E_BINS[:-1])  # [MeV]
    _X_MIN = 0.5  # [km.w.e.]
    _X_MAX = 14  # [km.w.e.]
    _SLANT_DEPTHS = np.linspace(
        _X_MIN, _X_MAX, int(2 * (_X_MAX - _X_MIN) + 1)
    )  # [km.w.e.]
    _ANGLES = np.degrees(np.arccos(0.5 / _SLANT_DEPTHS))  # [degrees]
    slant_depths = _SLANT_DEPTHS  # [km.w.e.]
    angles = _ANGLES  # [degrees]
    ANGLES_FOR_S_FLUXES = np.linspace(0, 90, 20)  # [degrees]

    # Other constants and variables

    global _MU_MASS
    global MONTHS
    global MONTHS_SNAMES
    global _survival_probability_tensor_configuration
    global _current_survival_probability_tensor
    global _mountain_loaded
    global shallow_extrapolation

    _MU_MASS = 105.6583745  # [MeV]
    MONTHS = [
        "January",
        "February",
        "March",
        "April",
        "May",
        "June",
        "July",
        "August",
        "September",
        "October",
        "November",
        "December",
    ]
    MONTHS_SNAMES = [
        "Jan",
        "Feb",
        "Mar",
        "Apr",
        "May",
        "Jun",
        "Jul",
        "Aug",
        "Sep",
        "Oct",
        "Nov",
        "Dec",
    ]
    _survival_probability_tensor_configuration = {}
    _current_survival_probability_tensor = None
    _mountain_loaded = False
    shallow_extrapolation = False

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
    _vertical_depth = _X_MIN
    _medium = "rock"
    _density = 2.65
    _n_muon = 100000

    # Global mountain variables
    # Set the variables to None first, in case they do not already exist in the namespace

    global _mountain_zenith_all
    global _mountain_azimuthal_all
    global _mountain_slant_depths_all
    global mountain

    _mountain_zenith_all = None
    _mountain_azimuthal_all = None
    _mountain_slant_depths_all = None
    mountain = None

    del _mountain_zenith_all
    del _mountain_azimuthal_all
    del _mountain_slant_depths_all
    del mountain

    gc.collect()
