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

import numpy as np
import scipy.integrate as scii
from tqdm import tqdm

import mute.constants as constants

# Calculate surface fluxes


def calc_s_fluxes(
    primary_model="GSF",
    interaction_model="SIBYLL-2.3c",
    atmosphere="CORSIKA",
    location="USStd",
    month=None,
    output=None,
    file_name="",
    force=False,
    test=False,
):

    """
    Calculate surface fluxes in units of [(cm^2 s sr MeV)^-1] for default surface energy grid and surface flux zenith angles.

    The default surface energy grid is given by constants.ENERGIES, and the default zenith angles are given by constants.ANGLES_FOR_S_FLUXES.

    Parameters
    ----------
    primary_model : str in {"GSF", "HG", "GH", "ZS", "ZSP", "PL27"} or tuple, optional (default: "GSF")
        The primary flux model to use in MCEq. Options:
        GSF  = GlobalSplineFitBeta
        HG   = HillasGaisser2012 (H3a)
        GH   = GaisserHonda
        ZS   = Zatsepin-Sokolskaya (Default)
        ZSP  = Zatsepin-Sokolskaya (PAMELA)
        PL27 = SimplePowerlaw27
        Alternatively, this can be set with a tuple. For example: (pm.GaisserStanevTilav, "3-gen").

    interaction_model : str, optional (default: "SIBYLL-2.3c")
        The hadronic interaction model to use in MCEq. See the Tutorial or MCEq documentation for a list of options.

    atmosphere : {"CORSIKA", "MSIS00"}, optional (default: "CORSIKA")
        The atmospheric model to use in MCEq. For US Standard Atmosphere, use "CORSIKA". For seasonal variations, use "MSIS00".

    location : str, optional (default: "USStd")
        The name of the location for which to calculate the surface fluxes. See the Tutorial or MCEq documentation for a list of options.

    month : str, optional (default: None)
        The month for which to calculate the surface fluxes. For US Standard Atmosphere, use None. For seasonal variations, use the month name.

    output : bool, optional (default: taken from constants.get_output())
        If True, an output file will be created to store the results.

    file_name : str, optional (default: constructed from input parameters)
        The name of the file in which to store the results. If output is False or None, this is ignored.

    force : bool, optional (default: False)
        If True, force the creation of new directories if required.

    test: bool, optional (default: False)
        For use in the file test_s_fluxes.py. If True, this will calculate surface fluxes for only three angles.

    Returns
    -------
    s_fluxes : NumPy ndarray
        A two-dimensional array containing the surface fluxes. The shape will be (91, 20), and the fluxes will be in units of [(cm^2 s sr MeV)^-1].
    """

    # Import packages

    from MCEq.core import MCEqRun
    import mceq_config
    import crflux.models as pm

    # Check values

    constants._check_constants(force=force)

    if output is None:
        output = constants.get_output()

    assert atmosphere in [
        "CORSIKA",
        "MSIS00",
    ], 'atmosphere must be set to either "CORSIKA" or "MSIS00".'

    primary_models = {
        "GSF": (pm.GlobalSplineFitBeta, None),
        "HG": (pm.HillasGaisser2012, "H3a"),
        "GH": (pm.GaisserHonda, None),
        "ZS": (pm.ZatsepinSokolskaya, "default"),
        "ZSP": (pm.ZatsepinSokolskaya, "pamela"),
        "PL27": (pm.SimplePowerlaw27, None),
    }

    if isinstance(primary_model, str):

        assert (
            primary_model in primary_models
        ), "Set primary model not available. See the available options in the Tutorial at {0}.".format(
            "https://github.com/wjwoodley/mute/blob/main/docs/Tutorial.md#changing-the-primary-model"
        )

        primary_model_for_MCEq = primary_models[primary_model]
        pm_sname = primary_model

    elif isinstance(primary_model, tuple) and len(primary_model) == 2:

        primary_model_for_MCEq = primary_model
        pm_sname = primary_model[0](primary_model[1]).sname

    else:

        raise TypeError(
            "Primary model not set correctly. For an explanation, see the Tutorial at {0}.".format(
                "https://github.com/wjwoodley/mute/blob/main/docs/Tutorial.md#changing-the-primary-model"
            )
        )

    # Set the angles for use in MCEq

    angles = constants.ANGLES_FOR_S_FLUXES

    if test:

        angles = [0, 30, 60]

    # Set MCEq up

    mceq_config.enable_default_tracking = False

    if constants.get_verbose() > 1:

        print(
            "Calculating surface fluxes for {0} using {1} and {2}.".format(
                location, pm_sname, interaction_model
            )
        )

        mceq_run = MCEqRun(
            interaction_model=interaction_model,
            primary_model=primary_model_for_MCEq,
            theta_deg=0.0,
            **mceq_config.config,
        )
        mceq_run.set_density_model((atmosphere, (location, month)))

    else:

        from contextlib import redirect_stdout

        with open(os.devnull, "w") as suppress, redirect_stdout(suppress):

            mceq_run = MCEqRun(
                interaction_model=interaction_model,
                primary_model=primary_model_for_MCEq,
                theta_deg=0.0,
                **mceq_config.config,
            )
            mceq_run.set_density_model((atmosphere, (location, month)))

    # Calculate the surface fluxes
    # Zeroth index = Surface energy
    # First index  = Zenith angle

    s_fluxes = np.zeros((len(constants.ENERGIES), len(angles)))

    # Run MCEq
    # Convert the surface fluxes from default [GeV] to [MeV]

    for j in (
        tqdm(range(len(angles))) if constants.get_verbose() >= 1 else range(len(angles))
    ):

        mceq_run.set_theta_deg(angles[j])
        mceq_run.solve()

        s_fluxes[:, j] = (
            1e-3
            * (
                mceq_run.get_solution("total_mu+", mag=0)
                + mceq_run.get_solution("total_mu-", mag=0)
            )[:-30]
        )

    if constants.get_verbose() > 1:
        print("Finished calculating surface fluxes.")

    # Write the results to the file

    if output:

        constants._check_directory(
            os.path.join(constants.get_directory(), "surface"), force=force
        )

        if file_name == "" or not isinstance(file_name, str):

            file_name = os.path.join(
                constants.get_directory(),
                "surface",
                "Surface_Fluxes_{0}_{1}_{2}_{3}.txt".format(
                    location, str(month), interaction_model, pm_sname
                ),
            )

        file_out = open(file_name, "w")

        for i in range(len(constants.ENERGIES)):

            for j in range(len(angles)):

                file_out.write(
                    "{0:1.14f} {1:1.5f} {2:1.14e}\n".format(
                        constants.ENERGIES[i],
                        angles[j],
                        s_fluxes[i, j],
                    )
                )

        file_out.close()

        if constants.get_verbose() > 1:
            print(f"Surface fluxes written to {file_name}.")

    return s_fluxes


# Get surface fluxes


def load_s_fluxes_from_file(
    primary_model="GSF",
    interaction_model="SIBYLL-2.3c",
    atmosphere="CORSIKA",
    location="USStd",
    month=None,
    file_name="",
    force=False,
    test=False,
):

    """
    Retrieve a surface fluxes matrix in units of [(cm^2 s sr MeV)^-1] stored in data/surface based on the input parameters.

    If file_name is not given, this function searches for a file name that matches the set location, month, interaction model, and primary model. If the file does not exist, prompt the user to run calc_s_fluxes().

    Parameters
    ----------
    primary_model : str in {"GSF", "HG", "GH", "ZS", "ZSP", "PL27"} or tuple, optional (default: "GSF")
        The primary flux model to use in MCEq. Options:
        GSF  = GlobalSplineFitBeta
        HG   = HillasGaisser2012 (H3a)
        GH   = GaisserHonda
        ZS   = Zatsepin-Sokolskaya (Default)
        ZSP  = Zatsepin-Sokolskaya (PAMELA)
        PL27 = SimplePowerlaw27
        Alternatively, this can be set with a tuple. For example: (pm.GaisserStanevTilav, "3-gen").

    interaction_model : str, optional (default: "SIBYLL-2.3c")
        The hadronic interaction model to use in MCEq. See the Tutorial or MCEq documentation for a list of options.

    atmosphere : {"CORSIKA", "MSIS00"}, optional (default: "CORSIKA")
        The atmospheric model to use in MCEq. For US Standard Atmosphere, use "CORSIKA". For seasonal variations, use "MSIS00".

    location : str, optional (default: "USStd")
        The name of the location for which to calculate the surface fluxes. See the Tutorial or MCEq documentation for a list of options.

    month : str, optional (default: None)
        The month for which to calculate the surface fluxes. For US Standard Atmosphere, use None. For seasonal variations, use the month name.

    force : bool, optional (default: False)
        If True, force the calculation of a new surface flux matrix and the creation of new directories if required.

    test: bool, optional (default: False)
        For use in the file test_s_fluxes.py. If True, this will calculate surface fluxes for only three angles.

    Returns
    -------
    s_fluxes : NumPy ndarray
        A two-dimensional array containing the surface fluxes. The shape will be (91, 20), and the fluxes will be in units of [(cm^2 s sr MeV)^-1].
    """

    # Check values

    constants._check_constants(force=force)

    if isinstance(primary_model, str):

        pass

    elif isinstance(primary_model, tuple) and len(primary_model) == 2:

        raise TypeError(
            "The primary model must be set with a string. To set it with a tuple, run the surface.calc_s_fluxes() function with output set to True. For an explanation, see the Tutorial at {0}.".format(
                "https://github.com/wjwoodley/mute/blob/main/docs/Tutorial.md#changing-the-primary-model"
            )
        )

    else:

        raise TypeError(
            "Primary model not set correctly. For an explanation, see the Tutorial at {0}.".format(
                "https://github.com/wjwoodley/mute/blob/main/docs/Tutorial.md#changing-the-primary-model"
            )
        )

    # Set the angles

    angles = constants.ANGLES_FOR_S_FLUXES

    if test:

        angles = [0, 30, 60]

    # Define a function to run if there is no surface fluxes file

    def no_file(force):

        # If the file does not exist, ask the user if they want to run MCEq to create it

        if not force:

            answer = input(
                "No surface flux matrix currently exists for these models. Would you like to create one (y/n)?: "
            )

        if force or answer.lower() == "y":

            s_fluxes = calc_s_fluxes(
                primary_model,
                interaction_model,
                atmosphere,
                location,
                month,
                force=force,
                test=test,
            )

            return s_fluxes

        else:

            print("Surface fluxes not calculated.")

            return None

    # Construct a file name based on the user's inputs if one has not been provided

    if file_name == "" or not isinstance(file_name, str):

        file_name = os.path.join(
            constants.get_directory(),
            "surface",
            "Surface_Fluxes_{0}_{1}_{2}_{3}.txt".format(
                location, str(month), interaction_model, primary_model
            ),
        )

    # Check if the file exists

    if os.path.exists(file_name):

        if constants.get_verbose() > 1:
            print(
                "Loading surface fluxes for {0} using {1} and {2}.".format(
                    location, primary_model, interaction_model
                )
            )

        # Check that the file has the correct numbers of energies and zenith angles

        file = open(file_name, "r")
        n_lines = len(file.read().splitlines())
        file.close()

        # If so, read the surface fluxes in from it

        if n_lines == len(constants.ENERGIES) * len(angles):

            s_fluxes = np.reshape(
                np.loadtxt(file_name)[:, 2], (len(constants.ENERGIES), len(angles))
            )

            if constants.get_verbose() > 1:

                print("Loaded surface fluxes.")

            return s_fluxes

        # If the file does not have the correct numbers of energies and zenith angles, run the no_file() function

        else:

            return no_file(force=force)

    # If the file does not exist, run the no_file() function

    else:

        return no_file(force=force)


# Calculate the surface intensities


def calc_s_intensities(
    s_fluxes=None,
    primary_model="GSF",
    interaction_model="SIBYLL-2.3c",
    atmosphere="CORSIKA",
    location="USStd",
    month=None,
    output=None,
    file_name="",
    force=False,
):

    """
    Calculate surface intensities in units of [(cm^2 s sr)^-1] for default surface energy grid and surface flux zenith angles.

    The default surface energy grid is given by constants.ENERGIES, and the default zenith angles are given by constants.ANGLES_FOR_S_FLUXES.

    Parameters
    ----------
    s_fluxes : NumPy ndarray, optional (default: taken from surface.load_s_fluxes_from_file())
        A surface flux matrix of shape (91, 20).

    primary_model : str in {"GSF", "HG", "GH", "ZS", "ZSP", "PL27"} or tuple, optional (default: "GSF")
        The primary flux model to use in MCEq. Options:
        GSF  = GlobalSplineFitBeta
        HG   = HillasGaisser2012 (H3a)
        GH   = GaisserHonda
        ZS   = Zatsepin-Sokolskaya (Default)
        ZSP  = Zatsepin-Sokolskaya (PAMELA)
        PL27 = SimplePowerlaw27
        Alternatively, this can be set with a tuple. For example: (pm.GaisserStanevTilav, "3-gen").

    interaction_model : str, optional (default: "SIBYLL-2.3c")
        The hadronic interaction model to use in MCEq. See the Tutorial or MCEq documentation for a list of options.

    atmosphere : {"CORSIKA", "MSIS00"}, optional (default: "CORSIKA")
        The atmospheric model to use in MCEq. For US Standard Atmosphere, use "CORSIKA". For seasonal variations, use "MSIS00".

    location : str, optional (default: "USStd")
        The name of the location for which to calculate the surface fluxes. See the Tutorial or MCEq documentation for a list of options.

    month : str, optional (default: None)
        The month for which to calculate the surface fluxes. For US Standard Atmosphere, use None. For seasonal variations, use the month name.

    output : bool, optional (default: taken from constants.get_output())
        If True, an output file will be created to store the results.

    file_name : str, optional (default: constructed from input parameters)
        The name of the file in which to store the results. If output is False or None, this is ignored.

    force : bool, optional (default: False)
        If True, force the calculation of a new surface flux matrix and the creation of new directories if required.

    Returns
    -------
    s_intensities : NumPy ndarray
        A one-dimensional array containing the surface intensities. The length will be 20, and the intensities will be in units of [(cm^2 s sr)^-1].
    """

    # Check values

    constants._check_constants()

    if output is None:
        output = constants.get_output()

    # Get the surface flux matrix

    if constants.get_verbose() > 1:
        print("Calculating surface intensities.")

    if s_fluxes is None:

        s_fluxes = load_s_fluxes_from_file(
            primary_model, interaction_model, atmosphere, location, month, force=force
        )

    # Check that the surface flux matrix has been loaded properly

    if s_fluxes is None:

        raise Exception(
            "Surface intensities not calculated. The surface flux matrix was not provided or loaded correctly."
        )

    s_fluxes = np.atleast_2d(s_fluxes)

    # Calculate the surface intensities

    s_intensities = [
        scii.simpson(s_fluxes[:, j], constants.ENERGIES)
        for j in range(len(constants.ANGLES_FOR_S_FLUXES))
    ]

    if constants.get_verbose() > 1:
        print("Finished calculating surface intensities.")

    # Write the results to the file

    if output:

        constants._check_directory(
            os.path.join(constants.get_directory(), "surface"), force=force
        )

        if file_name == "" or not isinstance(file_name, str):

            file_name = os.path.join(
                constants.get_directory(),
                "surface",
                "Surface_Intensities_{0}_{1}_{2}_{3}.txt".format(
                    location, str(month), interaction_model, primary_model
                ),
            )

        file_out = open(file_name, "w")

        for j in range(len(constants.ANGLES_FOR_S_FLUXES)):

            file_out.write(
                "{0:1.5f} {1:1.14e}\n".format(
                    constants.ANGLES_FOR_S_FLUXES[j], s_intensities[j]
                )
            )

        file_out.close()

        if constants.get_verbose() > 1:
            print(f"Surface intensities written to {file_name}.")

    return s_intensities


# Calculate total surface fluxes


def calc_s_tot_flux(
    s_fluxes=None,
    primary_model="GSF",
    interaction_model="SIBYLL-2.3c",
    atmosphere="CORSIKA",
    location="USStd",
    month=None,
    force=False,
):

    """
    Calculate a total surface flux in units of [(cm^2 s)^-1] for default surface energy grid and surface flux zenith angles.

    The default surface energy grid is given by constants.ENERGIES, and the default zenith angles are given by constants.ANGLES_FOR_S_FLUXES.

    Parameters
    ----------
    s_fluxes : NumPy ndarray, optional (default: taken from surface.load_s_fluxes_from_file())
        A surface flux matrix of shape (91, 20).

    primary_model : str in {"GSF", "HG", "GH", "ZS", "ZSP", "PL27"} or tuple, optional (default: "GSF")
        The primary flux model to use in MCEq. Options:
        GSF  = GlobalSplineFitBeta
        HG   = HillasGaisser2012 (H3a)
        GH   = GaisserHonda
        ZS   = Zatsepin-Sokolskaya (Default)
        ZSP  = Zatsepin-Sokolskaya (PAMELA)
        PL27 = SimplePowerlaw27
        Alternatively, this can be set with a tuple. For example: (pm.GaisserStanevTilav, "3-gen").

    interaction_model : str, optional (default: "SIBYLL-2.3c")
        The hadronic interaction model to use in MCEq. See the Tutorial or MCEq documentation for a list of options.

    atmosphere : {"CORSIKA", "MSIS00"}, optional (default: "CORSIKA")
        The atmospheric model to use in MCEq. For US Standard Atmosphere, use "CORSIKA". For seasonal variations, use "MSIS00".

    location : str, optional (default: "USStd")
        The name of the location for which to calculate the surface fluxes. See the Tutorial or MCEq documentation for a list of options.

    month : str, optional (default: None)
        The month for which to calculate the surface fluxes. For US Standard Atmosphere, use None. For seasonal variations, use the month name.

    force : bool, optional (default: False)
        If True, force the calculation of new arrays or matrices and the creation of new directories if required.

    Returns
    -------
    s_tot_flux : float
        The total surface flux in units of [(cm^2 s)^-1].
    """

    # Check values

    constants._check_constants()

    # Calculate the surface intensities

    if constants.get_verbose() > 1:
        print("Calculating total surface flux.")

    s_intensities = calc_s_intensities(
        s_fluxes,
        primary_model,
        interaction_model,
        atmosphere,
        location,
        month,
        force=force,
    )

    if s_intensities is None:

        raise Exception(
            "Total surface flux not calculated. The surface intensities were not calculated properly."
        )

    # Calculate the total surface flux
    # Because constants.ANGLES_FOR_S_FLUXES goes (0..89), cos(angles) goes (1..0)
    # cos(constants.ANGLES_FOR_S_FLUXES) is decreasing, but scii.simpson() wants an increasing array
    # Therefore, integrate backwards, using [::-1] on the integrand and steps
    # Otherwise, the answer will be negative

    s_tot_flux = (
        2
        * np.pi
        * scii.simpson(
            s_intensities[::-1], np.cos(np.radians(constants.ANGLES_FOR_S_FLUXES[::-1]))
        )
    )

    if constants.get_verbose() > 1:
        print("Finished calculating total surface flux.")

    return s_tot_flux
