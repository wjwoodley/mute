#########################
#########################
###                   ###
###  MUTE             ###
###  William Woodley  ###
###  24 May 2025      ###
###                   ###
#########################
#########################

# Import packages

import os

import numpy as np
import scipy.integrate as scii
from tqdm import tqdm
import warnings

import mute.constants as constants

# Calculate surface fluxes


def calc_s_fluxes(
    model="mceq",
    output=None,
    file_name="",
    force=False,
    test=False,
    angles=constants.ANGLES_FOR_S_FLUXES,
    **kwargs,
):
    """
    Calculate surface fluxes in units of [(cm^2 s sr MeV)^-1] for default surface energy grid and surface flux zenith angles.

    The default surface energy grid is given by constants.ENERGIES, and the default zenith angles are given by constants.ANGLES_FOR_S_FLUXES.

    Parameters
    ----------
    model : str in {"daemonflux", "mceq"}, optional (default: "mceq")
        The model to use to calculate surface fluxes. MCeq provides the option to specify primary, interaction, and density model keyword arguments (see Other Parameters), while daemonflux uses "gsf" as the primary model and "ddm" as the interaction model. The default model will be changed to "daemonflux" in v3.1.0. This parameter is case-insensitive.

    output : bool, optional (default: taken from constants.get_output())
        If True, an output file will be created to store the results.

    file_name : str, optional (default: constructed from input parameters)
        The name of the file in which to store the results. If output is False or None, this is ignored.

    force : bool, optional (default: False)
        If True, force the creation of new directories if required.

    test: bool, optional (default: False)
        For use in the file test_s_fluxes.py. If True, this will calculate surface fluxes for only three angles.

    angles : array-like, optional (default: taken from consants.ANGLES_FOR_S_FLUXES)
        An array of zenith angles in [degrees] to calculate the surface fluxes for.

    Other Parameters
    -----------------
    return_error : bool, optional (default: False)
        If True, if model == "daemonflux", return the error on the surface fluxes instead of the surface fluxes.

    primary_model : str in {"gsf", "hg", "h3a", "h4a", "gh", "gst3", "gst4", "zs", "zsp", "pl27"} or tuple, optional (default: "gsf")
        The primary flux model to use in MCEq. This parameter is case-insensitive. Options:
        gsf  = GlobalSplineFitBeta
        hg   = HillasGaisser2012 (H3a)
        h3a  = HillasGaisser2012 (H3a)
        h4a  = HillasGaisser2012 (H4a)
        gh   = GaisserHonda
        gst3 = GaisserStanevTilav (3-gen)
        gst4 = GaisserStanevTilav (4-gen)
        zs   = Zatsepin-Sokolskaya (Default)
        zsp  = Zatsepin-Sokolskaya (PAMELA)
        pl27 = SimplePowerlaw27
        Alternatively, this can be set with a tuple. For example: (pm.GaisserStanevTilav, "3-gen").

    interaction_model : str, optional (default: "sibyll23c")
        The hadronic interaction model to use in MCEq. See the Tutorial or MCEq documentation for a list of options. This parameter is case-insensitive.

    atmosphere : {"corsika", "msis00"}, optional (default: "corsika")
        The atmospheric model to use in MCEq. For US Standard Atmosphere, use "corsika". For seasonal variations, use "msis00". This parameter is case-insensitive.

    location : str or tuple, optional (default: "USStd")
        The name of the location for which to calculate the surface fluxes. See the Tutorial or MCEq documentation for a list of options. Alternatively, this can be set with a tuple of shape (latitude, longitude). This parameter is case-sensitive.

    month : str, optional (default: None)
        The month for which to calculate the surface fluxes. For US Standard Atmosphere, use None. For seasonal variations, use the month name. This parameter is case-sensitive.

    Returns
    -------
    s_fluxes : NumPy ndarray
        A two-dimensional array containing the surface fluxes (or the errors thereof). The default shape will be (91, 20), and the fluxes will be in units of [(cm^2 s sr MeV)^-1].
    """

    # Check values

    constants._check_constants(force=force)

    if output is None:
        output = constants.get_output()

    model = model.lower()

    assert model in [
        "daemonflux",
        "mceq",
    ], "The model must either be 'daemonflux' or 'mceq'."

    # Set the angles

    test_file_name = ""

    if test:

        angles = [0, 30, 60]
        test_file_name = "_pytest"

    # Calculate surface fluxes with DAEMONFLUX

    if model == "daemonflux":

        # Import the Flux class

        import daemonflux

        # Pop the keyword arguments
        # This keyword must be popped instead of gotten, since it should not be passed into the DAEMONFLUX Flux() function

        return_error = kwargs.pop("return_error", False)

        # Calculate surface fluxes and errors

        df_flux = daemonflux.Flux(**kwargs)

        e_grid = 1e-3 * constants.ENERGIES
        e_grid_div = (
            np.reshape(np.repeat(e_grid, len(angles)), (len(e_grid), len(angles)))
        ) ** 3

        if constants.get_verbose() > 1:
            print("Calculating surface fluxes.")

        # Calculate surface fluxes or surface flux errors

        if not return_error:

            s_fluxes = 1e-3 * (
                df_flux.flux(e_grid, angles, "total_muflux") / e_grid_div
            )

        else:

            s_fluxes = 1e-3 * (
                df_flux.error(e_grid, angles, "total_muflux") / e_grid_div
            )

        if constants.get_verbose() > 1:
            print("Finished calculating surface fluxes.")

    # Calculate surface fluxes with MCEq

    elif model == "mceq":

        # Import packages

        import MCEq
        from MCEq.core import MCEqRun
        import crflux.models as pm

        # Get the keyword arguments

        primary_model = kwargs.get("primary_model", "gsf")
        interaction_model = (
            kwargs.get("interaction_model", "sibyll23c")
            .replace("-", "")
            .replace(".", "")
            .lower()
        )
        atmosphere = kwargs.get("atmosphere", "corsika").lower()
        location = kwargs.get("location", "USStd")
        month = kwargs.get("month", None)

        # Check values

        return_error = kwargs.pop("return_error", None)

        if return_error is not None:
            raise ValueError(
                'Errors cannot be calculated with MCEq. To include errors in your calculation, please use daemonflux by setting model = "daemonflux".'
            )

        assert atmosphere in [
            "corsika",
            "msis00",
        ], 'atmosphere must be set to either "corsika" or "msis00".'

        primary_models = {
            "gsf": (pm.GlobalSplineFitBeta, None),
            "hg": (pm.HillasGaisser2012, "H3a"),
            "h3a": (pm.HillasGaisser2012, "H3a"),
            "h4a": (pm.HillasGaisser2012, "H4a"),
            "gh": (pm.GaisserHonda, None),
            "gst3": (pm.GaisserStanevTilav, "3-gen"),
            "gst4": (pm.GaisserStanevTilav, "4-gen"),
            "zs": (pm.ZatsepinSokolskaya, "default"),
            "zsp": (pm.ZatsepinSokolskaya, "pamela"),
            "pl27": (pm.SimplePowerlaw27, None),
        }

        if isinstance(primary_model, str):

            assert (
                primary_model.lower() in primary_models
            ), "Set primary model not available. See the available options in the Tutorial at {0}.".format(
                "https://github.com/wjwoodley/mute/blob/main/docs/Tutorial_Models.md#primary-model"
            )

            if primary_model == "hg":

                if constants.get_verbose() >= 1:

                    warnings.warn(
                        'The "hg" option is deprecated. This option will be removed in v3.1.0. Please use "h3a" or "h4a".',
                        DeprecationWarning,
                        stacklevel=2,
                    )

            primary_model_for_MCEq = primary_models[primary_model.lower()]
            pm_sname = primary_model.lower()

        elif isinstance(primary_model, tuple) and len(primary_model) == 2:

            primary_model_for_MCEq = primary_model
            pm_sname = primary_model[0](primary_model[1]).sname.lower()

        else:

            raise TypeError(
                "Primary model not set correctly. For an explanation, see the Tutorial at {0}.".format(
                    "https://github.com/wjwoodley/mute/blob/main/docs/Tutorial_Models.md#primary-model"
                )
            )

        # Check that the location is set correctly

        location_is_str = isinstance(location, str)
        location_is_tuple = isinstance(location, tuple) and len(location) == 2

        if not (location_is_str or location_is_tuple):

            raise TypeError(
                "Location not set correctly. For an explanation, see the Tutorial at {0}.".format(
                    "https://github.com/wjwoodley/mute/blob/main/docs/Tutorial_Models.md#atmospheric-model"
                )
            )

        if (location_is_str and location != "USStd") or location_is_tuple:

            assert (
                interaction_model != "ddm" and interaction_model != "sibyll23d"
            ), "{0} is not currently supported for varying atmospheres in MCEq. Please use a different interaction model or US Standard Atmosphere.".format(
                interaction_model
            )
            assert (
                atmosphere == "msis00"
            ), 'atmosphere must be set to "msis00" if location is specified.'
            assert (
                month is not None
            ), "month must be specified if location is specified."

        # Set MCEq up

        MCEq.config.enable_default_tracking = False

        if constants.get_verbose() > 1:

            print(
                "Calculating surface fluxes for {0} using {1} and {2}.".format(
                    location, pm_sname, interaction_model
                )
            )
            MCEq.config.debug_level = 1

        else:

            MCEq.config.debug_level = 0

        if isinstance(location, str):

            mceq_run = MCEqRun(
                interaction_model=interaction_model,
                primary_model=primary_model_for_MCEq,
                theta_deg=0.0,
            )

            mceq_run.set_density_model((atmosphere.upper(), (location, month)))

        elif isinstance(location, tuple):

            MCEq.config.density_model = ("MSIS00", ("SouthPole", "January"))

            mceq_run = MCEqRun(
                interaction_model=interaction_model,
                primary_model=primary_model_for_MCEq,
                theta_deg=0.0,
            )

            density_model_MSIS00 = mceq_run.density_model
            density_model_MSIS00.set_location_coord(*location[::-1])
            density_model_MSIS00.set_season(month)

        # Calculate the surface fluxes
        # Zeroth index = Surface energy
        # First index  = Zenith angle

        s_fluxes = np.zeros((len(constants.ENERGIES), len(angles)))

        # Run MCEq
        # Convert the surface fluxes from default [GeV] to [MeV]

        for j in (
            tqdm(range(len(angles)))
            if constants.get_verbose() >= 1
            else range(len(angles))
        ):

            # Solve for the given zenith angle

            mceq_run.set_theta_deg(angles[j])
            mceq_run.solve()

            # Store the results in the s_fluxes matrix

            s_fluxes[:, j] = (
                1e-3
                * (
                    mceq_run.get_solution("total_mu+", mag=0)
                    + mceq_run.get_solution("total_mu-", mag=0)
                )[:-30]
            )

        if constants.get_verbose() > 1:
            print("Finished calculating surface fluxes.")

    else:

        raise NotImplementedError(
            'Model of type {0} is not available. The only options are "daemonflux" and "mceq".'.format(
                model
            )
        )

    # Write the results to the file

    if output:

        constants._check_directory(
            os.path.join(constants.get_directory(), "surface"), force=force
        )

        if file_name == "" or not isinstance(file_name, str):

            if model == "daemonflux":

                error_file = "" if not return_error else "_error"
                file_name = os.path.join(
                    constants.get_directory(),
                    "surface",
                    "surface_fluxes_daemonflux{0}{1}.txt".format(
                        error_file, test_file_name
                    ),
                )

            elif model == "mceq":

                file_name = os.path.join(
                    constants.get_directory(),
                    "surface",
                    "surface_fluxes_{0}_{1}_{2}_{3}{4}.txt".format(
                        location,
                        str(month),
                        interaction_model,
                        pm_sname,
                        test_file_name,
                    ),
                )

            else:

                raise NotImplementedError(
                    'Model of type {0} is not available. The only options are "daemonflux" and "mceq".'.format(
                        model
                    )
                )

        file_out = open(file_name, "w")

        for i in range(len(constants.ENERGIES)):

            for j in range(len(angles)):

                file_out.write(
                    "{0:1.14f} {1:1.5f} {2:1.14e}\n".format(
                        constants.ENERGIES[i], angles[j], s_fluxes[i, j]
                    )
                )

        file_out.close()

        if constants.get_verbose() > 1:
            print(f"Surface fluxes written to {file_name}.")

    return s_fluxes


# Get surface fluxes


def load_s_fluxes_from_file(
    model="mceq", file_name="", force=False, test=False, **kwargs
):
    """
    Retrieve a surface fluxes matrix in units of [(cm^2 s sr MeV)^-1] stored in data/surface based on the input parameters.

    If file_name is not given, this function searches for a file name that matches the set location, month, interaction model, and primary model. If the file does not exist, prompt the user to run calc_s_fluxes().

    Parameters
    ----------
    model : str in {"daemonflux", "mceq"}, optional (default: "mceq")
        The model to use to calculate surface fluxes. MCeq provides the option to specify primary, interaction, and density model keyword arguments (see Other Parameters), while daemonflux uses "gsf" as the primary model and "DDM" as the interaction model. The default model will be changed to "daemonflux" in v3.1.0. This parameter is case-insensitive.

    file_name : str, optional (default: constructed from input parameters)
        The name of the file from which to load the surface fluxes. This must be the full path to the file.

    force : bool, optional (default: False)
        If True, force the calculation of a new surface flux matrix and the creation of new directories if required.

    test: bool, optional (default: False)
        For use in the file test_s_fluxes.py. If True, this will calculate surface fluxes for only three angles.

    Other Parameters
    -----------------
    return_error : bool, optional (default: False)
        If True, if model == "daemonflux", return the error on the surface fluxes instead of the surface fluxes.

    primary_model : str in {"gsf", "hg", "h3a", "h4a", "gh", "gst3", "gst4", "zs", "zsp", "pl27"} or tuple, optional (default: "gsf")
        The primary flux model to use in MCEq. This parameter is case-insensitive. Options:
        gsf  = GlobalSplineFitBeta
        hg   = HillasGaisser2012 (H3a)
        h3a  = HillasGaisser2012 (H3a)
        h4a  = HillasGaisser2012 (H4a)
        gh   = GaisserHonda
        gst3 = GaisserStanevTilav (3-gen)
        gst4 = GaisserStanevTilav (4-gen)
        zs   = Zatsepin-Sokolskaya (Default)
        zsp  = Zatsepin-Sokolskaya (PAMELA)
        pl27 = SimplePowerlaw27
        Alternatively, this can be set with a tuple. For example: (pm.GaisserStanevTilav, "3-gen").

    interaction_model : str, optional (default: "sibyll23c")
        The hadronic interaction model to use in MCEq. See the Tutorial or MCEq documentation for a list of options. This parameter is case-insensitive.

    atmosphere : {"corsika", "msis00"}, optional (default: "corsika")
        The atmospheric model to use in MCEq. For US Standard Atmosphere, use "corsika". For seasonal variations, use "msis00". This parameter is case-insensitive.

    location : str or tuple, optional (default: "USStd")
        The name of the location for which to calculate the surface fluxes. See the Tutorial or MCEq documentation for a list of options. Alternatively, this can be set with a tuple of shape (latitude, longitude). This parameter is case-sensitive.

    month : str, optional (default: None)
        The month for which to calculate the surface fluxes. For US Standard Atmosphere, use None. For seasonal variations, use the month name. This parameter is case-sensitive.

    Returns
    -------
    s_fluxes : NumPy ndarray
        A two-dimensional array containing the surface fluxes. The shape will be (91, 20), and the fluxes will be in units of [(cm^2 s sr MeV)^-1].
    """

    # Check values

    constants._check_constants(force=force)

    model = model.lower()

    if model == "daemonflux":

        # Get the keyword arguments

        return_error = kwargs.get("return_error", False)

    elif model == "mceq":

        # Get the keyword arguments

        primary_model = kwargs.get("primary_model", "gsf")
        interaction_model = (
            kwargs.get("interaction_model", "sibyll23c")
            .replace("-", "")
            .replace(".", "")
            .lower()
        )
        atmosphere = kwargs.get("atmosphere", "corsika").lower()
        location = kwargs.get("location", "USStd")
        month = kwargs.get("month", None)

        # Check values

        if isinstance(primary_model, str):

            pass

        elif isinstance(primary_model, tuple):

            raise TypeError(
                "The primary model must be set with a string. To set it with a tuple, run the surface.calc_s_fluxes() function with output set to True. For an explanation, see the Tutorial at {0}.".format(
                    "https://github.com/wjwoodley/mute/blob/main/docs/Tutorial_Models.md#primary-model"
                )
            )

        else:

            raise TypeError(
                "Primary model not set correctly. For an explanation, see the Tutorial at {0}.".format(
                    "https://github.com/wjwoodley/mute/blob/main/docs/Tutorial_Models.md#primary-model"
                )
            )

    else:

        raise NotImplementedError(
            'Model of type {0} is not available. The only options are "daemonflux" and "mceq".'.format(
                model
            )
        )

    # Set the angles

    angles = constants.ANGLES_FOR_S_FLUXES

    test_file_name = ""

    if test:
        angles = [0, 30, 60]
        test_file_name = "_pytest"

    # Define a function to run if there is no surface fluxes file

    def no_file(force):

        # If the file does not exist, ask the user if they want to run MCEq to create it

        if not force:
            answer = input(
                "No surface flux matrix currently exists for these models. Would you like to create one (y/n)?: "
            )

        if force or answer.lower() == "y":

            s_fluxes = calc_s_fluxes(model=model, force=force, test=test, **kwargs)

            return s_fluxes

        else:

            print("Surface fluxes not calculated.")

            return None

    # Construct a file name based on the user's inputs if one has not been provided

    if file_name == "" or not isinstance(file_name, str):

        if model == "daemonflux":

            error_file = "" if not return_error else "_error"
            file_name = os.path.join(
                constants.get_directory(),
                "surface",
                "surface_fluxes_daemonflux{0}{1}.txt".format(
                    error_file, test_file_name
                ),
            )

        elif model == "mceq":

            file_name = os.path.join(
                constants.get_directory(),
                "surface",
                "surface_fluxes_{0}_{1}_{2}_{3}{4}.txt".format(
                    location,
                    str(month),
                    interaction_model,
                    primary_model.lower(),
                    test_file_name,
                ),
            )

    # Check if the file exists

    if os.path.exists(file_name):

        if constants.get_verbose() > 1:

            if model == "daemonflux":

                print("Loading surface fluxes for daemonflux.")

            elif model == "mceq":

                print(
                    "Loading surface fluxes for {0} using {1} and {2}.".format(
                        location, primary_model.lower(), interaction_model
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
    output=None,
    file_name="",
    force=False,
    **kwargs,
):
    """
    Calculate surface intensities in units of [(cm^2 s sr)^-1] for default surface energy grid and surface flux zenith angles.

    The default surface energy grid is given by constants.ENERGIES, and the default zenith angles are given by constants.ANGLES_FOR_S_FLUXES.

    Parameters
    ----------
    output : bool, optional (default: taken from constants.get_output())
        If True, an output file will be created to store the results.

    file_name : str, optional (default: constructed from input parameters)
        The name of the file in which to store the results. If output is False or None, this is ignored.

    force : bool, optional (default: False)
        If True, force the calculation of a new surface flux matrix and the creation of new directories if required.

    Other Parameters
    ----------------
    s_fluxes : NumPy ndarray, optional (default: taken from surface.load_s_fluxes_from_file())
        A surface flux matrix of shape (91, 20).

    model : str in {"daemonflux", "mceq"}, optional (default: "mceq")
        The model to use to calculate surface fluxes. MCeq provides the option to specify primary, interaction, and density model keyword arguments, while DAEMONFLUX uses "gsf" as the primary model and "ddm" as the interaction model. The default model will be changed to "daemonflux" in v3.1.0. This parameter is case-insensitive.

    return_error : bool, optional (default: False)
        If True, if model == "daemonflux", return the error on the surface fluxes instead of the surface fluxes.

    primary_model : str in {"gsf", "hg", "h3a", "h4a", "gh", "gst3", "gst4", "zs", "zsp", "pl27"} or tuple, optional (default: "gsf")
        The primary flux model to use in MCEq. This parameter is case-insensitive. Options:
        gsf  = GlobalSplineFitBeta
        hg   = HillasGaisser2012 (H3a)
        h3a  = HillasGaisser2012 (H3a)
        h4a  = HillasGaisser2012 (H4a)
        gh   = GaisserHonda
        gst3 = GaisserStanevTilav (3-gen)
        gst4 = GaisserStanevTilav (4-gen)
        zs   = Zatsepin-Sokolskaya (Default)
        zsp  = Zatsepin-Sokolskaya (PAMELA)
        pl27 = SimplePowerlaw27
        Alternatively, this can be set with a tuple. For example: (pm.GaisserStanevTilav, "3-gen").

    interaction_model : str, optional (default: "sibyll23c")
        The hadronic interaction model to use in MCEq. See the Tutorial or MCEq documentation for a list of options. This parameter is case-insensitive.

    atmosphere : {"corsika", "msis00"}, optional (default: "corsika")
        The atmospheric model to use in MCEq. For US Standard Atmosphere, use "corsika". For seasonal variations, use "msis00". This parameter is case-insensitive.

    location : str or tuple, optional (default: "USStd")
        The name of the location for which to calculate the surface fluxes. See the Tutorial or MCEq documentation for a list of options. Alternatively, this can be set with a tuple of shape (latitude, longitude). This parameter is case-sensitive.

    month : str, optional (default: None)
        The month for which to calculate the surface fluxes. For US Standard Atmosphere, use None. For seasonal variations, use the month name. This parameter is case-sensitive.

    Returns
    -------
    s_intensities : NumPy ndarray
        A one-dimensional array containing the surface intensities. The length will be 20, and the intensities will be in units of [(cm^2 s sr)^-1].
    """

    # Check values

    constants._check_constants()

    if output is None:
        output = constants.get_output()

    # Get the keyword arguments

    s_fluxes = kwargs.get("s_fluxes", None)
    model = kwargs.get("model", "mceq").lower()

    # Set the model to None if an s_fluxes matrix is provided by the user
    # This will prevent model information from being used in the construction of file_name

    if s_fluxes is not None:
        model = None

    # Get the keyword arguments
    # Construct output file names for use if output is True

    if output:
        constants._check_directory(
            os.path.join(constants.get_directory(), "surface"), force=force
        )

    if model == "daemonflux":

        return_error = kwargs.get("return_error", False)

        if output:

            if file_name == "" or not isinstance(file_name, str):

                error_file = "" if not return_error else "_error"
                file_name = os.path.join(
                    constants.get_directory(),
                    "surface",
                    "surface_intensities_daemonflux{0}.txt".format(error_file),
                )

    elif model == "mceq":

        primary_model = kwargs.get("primary_model", "gsf").lower()
        interaction_model = (
            kwargs.get("interaction_model", "sibyll23c")
            .replace("-", "")
            .replace(".", "")
            .lower()
        )
        location = kwargs.get("location", "USStd")
        month = kwargs.get("month", None)

        if output:

            if file_name == "" or not isinstance(file_name, str):

                file_name = os.path.join(
                    constants.get_directory(),
                    "surface",
                    "surface_intensities_{0}_{1}_{2}_{3}.txt".format(
                        location, str(month), interaction_model, primary_model
                    ),
                )

    elif model is None:

        if output:

            if file_name == "" or not isinstance(file_name, str):

                file_name = os.path.join(
                    constants.get_directory(), "surface", "surface_intensities.txt"
                )

    else:

        NotImplementedError(
            'Model of type {0} is not available. The only options are "daemonflux" and "mceq".'.format(
                model
            )
        )

    # Get the surface flux matrix

    if constants.get_verbose() > 1:
        print("Calculating surface intensities.")

    if s_fluxes is None:
        s_fluxes = load_s_fluxes_from_file(force=force, **kwargs)

    # Check that the surface flux matrix has been loaded properly

    if s_fluxes is None:
        raise Exception(
            "Surface intensities not calculated. The surface flux matrix was not provided or loaded correctly."
        )

    s_fluxes = np.atleast_2d(s_fluxes)

    # Calculate the surface intensities

    s_intensities = scii.simpson(s_fluxes, x=constants.ENERGIES, axis=0)

    if constants.get_verbose() > 1:
        print("Finished calculating surface intensities.")

    # Write the results to the file

    if output:

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


# Calculate surface energy spectra


def calc_s_e_spect(output=None, file_name="", force=False, **kwargs):
    """
    Calculate a surface energy spectrum in units of [(cm^2 s MeV)^-1] for default surface energy grid.

    The default surface energy grid is given by constants.ENERGIES.

    Parameters
    ----------
    output : bool, optional (default: taken from constants.get_output())
        If True, an output file will be created to store the results.

    file_name : str, optional (default: constructed from input parameters)
        The name of the file in which to store the results. If output is False or None, this is ignored.

    force : bool, optional (default: False)
        If True, force the calculation of a new surface flux matrix and the creation of new directories if required.

    Other Parameters
    ----------------
    s_fluxes : NumPy ndarray, optional (default: taken from surface.load_s_fluxes_from_file())
        A surface flux matrix of shape (91, 20).

    model : str in {"daemonflux", "mceq"}, optional (default: "mceq")
        The model to use to calculate surface fluxes. MCeq provides the option to specify primary, interaction, and density model keyword arguments, while DAEMONFLUX uses "gsf" as the primary model and "ddm" as the interaction model. The default model will be changed to "daemonflux" in v3.1.0. This parameter is case-insensitive.

    return_error : bool, optional (default: False)
        If True, if model == "daemonflux", return the error on the surface fluxes instead of the surface fluxes.

    primary_model : str in {"gsf", "hg", "h3a", "h4a", "gh", "gst3", "gst4", "zs", "zsp", "pl27"} or tuple, optional (default: "gsf")
        The primary flux model to use in MCEq. This parameter is case-insensitive. Options:
        gsf  = GlobalSplineFitBeta
        hg   = HillasGaisser2012 (H3a)
        h3a  = HillasGaisser2012 (H3a)
        h4a  = HillasGaisser2012 (H4a)
        gh   = GaisserHonda
        gst3 = GaisserStanevTilav (3-gen)
        gst4 = GaisserStanevTilav (4-gen)
        zs   = Zatsepin-Sokolskaya (Default)
        zsp  = Zatsepin-Sokolskaya (PAMELA)
        pl27 = SimplePowerlaw27
        Alternatively, this can be set with a tuple. For example: (pm.GaisserStanevTilav, "3-gen").

    interaction_model : str, optional (default: "sibyll23c")
        The hadronic interaction model to use in MCEq. See the Tutorial or MCEq documentation for a list of options. This parameter is case-insensitive.

    atmosphere : {"corsika", "msis00"}, optional (default: "corsika")
        The atmospheric model to use in MCEq. For US Standard Atmosphere, use "corsika". For seasonal variations, use "msis00". This parameter is case-insensitive.

    location : str or tuple, optional (default: "USStd")
        The name of the location for which to calculate the surface fluxes. See the Tutorial or MCEq documentation for a list of options. Alternatively, this can be set with a tuple of shape (latitude, longitude). This parameter is case-sensitive.

    month : str, optional (default: None)
        The month for which to calculate the surface fluxes. For US Standard Atmosphere, use None. For seasonal variations, use the month name. This parameter is case-sensitive.

    Returns
    -------
    s_e_spect : NumPy ndarray
        A one-dimensional array containing the surface energy spectrum. The length will be 91, and the energy spectrum values will be in units of [(cm^2 s GeV)^-1].
    """

    # Check values

    constants._check_constants()

    if output is None:
        output = constants.get_output()

    # Get the keyword arguments

    s_fluxes = kwargs.get("s_fluxes", None)
    model = kwargs.get("model", "mceq").lower()

    # Set the model to None if an s_fluxes matrix is provided by the user
    # This will prevent model information from being used in the construction of file_name

    if s_fluxes is not None:
        model = None

    # Get the keyword arguments
    # Construct output file names for use if output is True

    if output:
        constants._check_directory(
            os.path.join(constants.get_directory(), "surface"), force=force
        )

    if model == "daemonflux":

        return_error = kwargs.get("return_error", False)

        if output:

            if file_name == "" or not isinstance(file_name, str):

                error_file = "" if not return_error else "_error"
                file_name = os.path.join(
                    constants.get_directory(),
                    "surface",
                    "surface_energy_spectrum_daemonflux{0}.txt".format(error_file),
                )

    elif model == "mceq":

        primary_model = kwargs.get("primary_model", "gsf").lower()
        interaction_model = (
            kwargs.get("interaction_model", "sibyll23c")
            .replace("-", "")
            .replace(".", "")
            .lower()
        )
        location = kwargs.get("location", "USStd")
        month = kwargs.get("month", None)

        if output:

            if file_name == "" or not isinstance(file_name, str):

                file_name = os.path.join(
                    constants.get_directory(),
                    "surface",
                    "surface_energy_spectrum_{0}_{1}_{2}_{3}.txt".format(
                        location, str(month), interaction_model, primary_model
                    ),
                )

    elif model is None:

        if output:

            if file_name == "" or not isinstance(file_name, str):

                file_name = os.path.join(
                    constants.get_directory(), "surface", "surface_energy_spectrum.txt"
                )

    else:

        NotImplementedError(
            'Model of type {0} is not available. The only options are "daemonflux" and "mceq".'.format(
                model
            )
        )

    # Get the surface flux matrix

    if constants.get_verbose() > 1:
        print("Calculating surface energy spectrum.")

    if s_fluxes is None:

        s_fluxes = load_s_fluxes_from_file(force=force, **kwargs)

    # Check that the surface flux matrix has been loaded properly

    if s_fluxes is None:
        raise Exception(
            "Surface intensities not calculated. The surface flux matrix was not provided or loaded correctly."
        )

    s_fluxes = np.atleast_2d(s_fluxes)

    # Calculate the surface energy spectrum

    s_e_spect = (
        2
        * np.pi
        * abs(
            scii.simpson(
                s_fluxes, x=np.cos(np.radians(constants.ANGLES_FOR_S_FLUXES)), axis=1
            )
        )
    )

    if constants.get_verbose() > 1:
        print("Finished calculating surface energy spectrum.")

    # Write the results to the file

    if output:

        file_out = open(file_name, "w")

        for i in range(len(constants.ENERGIES)):

            file_out.write(
                "{0:1.14f} {1:1.14e}\n".format(constants.ENERGIES[i], s_e_spect[i])
            )

        file_out.close()

        if constants.get_verbose() > 1:
            print(f"Surface energy spectrum written to {file_name}.")

    return s_e_spect


# Calculate total surface fluxes


def calc_s_tot_flux(
    s_fluxes=None,
    model="mceq",
    force=False,
    **kwargs,
):
    """
    Calculate a total surface flux in units of [(cm^2 s)^-1] for default surface energy grid and surface flux zenith angles.

    The default surface energy grid is given by constants.ENERGIES, and the default zenith angles are given by constants.ANGLES_FOR_S_FLUXES.

    Parameters
    ----------
    s_fluxes : NumPy ndarray, optional (default: taken from surface.load_s_fluxes_from_file())
        A surface flux matrix of shape (91, 20).

    model : str in {"daemonflux", "mceq"}, optional (default: "mceq")
        The model to use to calculate surface fluxes. MCeq provides the option to specify primary, interaction, and density model keyword arguments (see Other Parameters), while daemonflux uses "gsf" as the primary model and "DDM" as the interaction model. The default model will be changed to "daemonflux" in v3.1.0. This parameter is case-insensitive.

    force : bool, optional (default: False)
        If True, force the calculation of new arrays or matrices and the creation of new directories if required.

    Other Parameters
    ----------------
    return_error : bool, optional (default: False)
        If True, if model == "daemonflux", return the error on the surface fluxes instead of the surface fluxes.

    primary_model : str in {"gsf", "hg", "h3a", "h4a", "gh", "gst3", "gst4", "zs", "zsp", "pl27"} or tuple, optional (default: "gsf")
        The primary flux model to use in MCEq. This parameter is case-insensitive. Options:
        gsf  = GlobalSplineFitBeta
        hg   = HillasGaisser2012 (H3a)
        h3a  = HillasGaisser2012 (H3a)
        h4a  = HillasGaisser2012 (H4a)
        gh   = GaisserHonda
        gst3 = GaisserStanevTilav (3-gen)
        gst4 = GaisserStanevTilav (4-gen)
        zs   = Zatsepin-Sokolskaya (Default)
        zsp  = Zatsepin-Sokolskaya (PAMELA)
        pl27 = SimplePowerlaw27
        Alternatively, this can be set with a tuple. For example: (pm.GaisserStanevTilav, "3-gen").

    interaction_model : str, optional (default: "sibyll23c")
        The hadronic interaction model to use in MCEq. See the Tutorial or MCEq documentation for a list of options. This parameter is case-insensitive.

    atmosphere : {"corsika", "msis00"}, optional (default: "corsika")
        The atmospheric model to use in MCEq. For US Standard Atmosphere, use "corsika". For seasonal variations, use "msis00". This parameter is case-insensitive.

    location : str or tuple, optional (default: "USStd")
        The name of the location for which to calculate the surface fluxes. See the Tutorial or MCEq documentation for a list of options. Alternatively, this can be set with a tuple of shape (latitude, longitude). This parameter is case-sensitive.

    month : str, optional (default: None)
        The month for which to calculate the surface fluxes. For US Standard Atmosphere, use None. For seasonal variations, use the month name. This parameter is case-sensitive.

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
        s_fluxes=s_fluxes, model=model, force=force, **kwargs
    )

    if s_intensities is None:
        raise Exception(
            "Total surface flux not calculated. The surface intensities were not calculated properly."
        )

    # Calculate the total surface flux

    s_tot_flux = float(
        2
        * np.pi
        * abs(
            scii.simpson(
                s_intensities, x=np.cos(np.radians(constants.ANGLES_FOR_S_FLUXES))
            )
        )
    )

    if constants.get_verbose() > 1:
        print("Finished calculating total surface flux.")

    return float(s_tot_flux)
