##########################
##########################
###                    ###
###  MUTE              ###
###  William Woodley   ###
###  10 November 2021  ###
###                    ###
##########################
##########################

# Import packages

import mute.constants as constants

import numpy as np
from tqdm import tqdm

# Calculate surface fluxes


def calc_s_fluxes(
    location="USStd",
    month=None,
    interaction_model="SIBYLL-2.3c",
    primary_model="GSF",
    atmosphere="CORSIKA",
    output=None,
    force=False,
):

    """
    Calculate surface fluxes for default surface energy grid and zenith angles.

    The default surface energy grid are given by constants.ENERGIES, and the default zenith angles are given by constants.ANGLES_FOR_S_FLUXES.

    Parameters
    ----------
    location : str, optional (default: "USStd")
        The name of the location for which to calculate the surface fluxes. See the Tutorial or MCEq documentation for a list of options.

    month : str, optional (default: None)
        The month for which to calculate the surface fluxes. For US Standard Atmosphere, use None. For seasonal variations, use the month name.

    interaction_model: str, optional (default: "SIBYLL-2.3c")
        The hadronic interaction model to use in MCEq. See the Tutorial or MCEq documentation for a list of options.

    primary_model : {"GSF", "HG", "GH"}, optional (default: "GSF")
        The primary flux model to use in MCEq. Options:
        GSF = GlobalSplineFitBeta
        HG  = HillasGaisser2012
        GH  = GaisserHonda

    atmosphere : {"CORSIKA", "MSIS00"}, optional (default: "CORSIKA")
        The atmospheric model. For US Standard Atmosphere, use "CORSIKA". For seasonal variations, use "MSIS00".

    output : bool, optional (default: taken from constants.get_output())
        If True, an output file will be created to store the results.

    force : bool
        If True, this will force the creation of a surface_fluxes directory if one does not already exist.

    Returns
    -------
    s_fluxes : NumPy ndarray
        A two-dimensional array containing the surface fluxes.
    """

    # Import MCEq

    from MCEq.core import MCEqRun
    import mceq_config
    import crflux.models as pm

    # Check values

    constants.check_constants(force=force)

    if output is None:
        output = constants.get_output()

    primary_models = {
        "GSF": (pm.GlobalSplineFitBeta, None),
        "HG": (pm.HillasGaisser2012, "H3a"),
        "GH": (pm.GaisserHonda, None),
    }

    assert primary_model in primary_models
    assert atmosphere in ["CORSIKA", "MSIS00"]

    # Set MCEq up

    if constants.get_verbose() > 1:

        print(
            "Calculating surface fluxes for "
            + location
            + " using "
            + interaction_model
            + " and "
            + primary_model
            + "."
        )

        mceq_run = MCEqRun(
            interaction_model=interaction_model,
            primary_model=primary_models[primary_model],
            theta_deg=0.0,
            **mceq_config.config
        )

        mceq_run.set_density_model((atmosphere, (location, month)))

    else:

        from contextlib import redirect_stdout

        with open(os.devnull, "w") as suppress, redirect_stdout(suppress):

            mceq_run = MCEqRun(
                interaction_model=interaction_model,
                primary_model=primary_models[primary_model],
                theta_deg=0.0,
                **mceq_config.config
            )

            mceq_run.set_density_model((atmosphere, (location, month)))

    # Run MCEq
    # Convert the surface fluxes from [GeV] to [MeV] to use in PROPOSAL

    s_fluxes = np.zeros((len(constants.ENERGIES), len(constants.ANGLES_FOR_S_FLUXES)))

    # Run MCEq

    for j in (
        tqdm(range(len(constants.ANGLES_FOR_S_FLUXES)))
        if constants.get_verbose() >= 1
        else range(len(constants.ANGLES_FOR_S_FLUXES))
    ):

        mceq_run.set_theta_deg(constants.ANGLES_FOR_S_FLUXES[j])
        mceq_run.solve()

        s_fluxes[:, j] = (
            mceq_run.get_solution("total_mu+", mag=0)
            + mceq_run.get_solution("total_mu-", mag=0)
        ) * 1e-3

    if constants.get_verbose() > 1:
        print("Finished calculating surface fluxes.")

    # Write the results to a file

    if output:

        constants.check_directory(
            constants.get_directory() + "/surface_fluxes", force=force
        )

        if month is None:

            file_name = (
                constants.get_directory()
                + "/surface_fluxes/Surface_Fluxes_"
                + location
                + "_"
                + interaction_model
                + "_"
                + primary_model
                + ".txt"
            )

        else:

            file_name = (
                constants.get_directory()
                + "/surface_fluxes/Surface_Fluxes_"
                + location
                + "_"
                + month
                + "_"
                + interaction_model
                + "_"
                + primary_model
                + ".txt"
            )

        file_out = open(file_name, "w")

        for i in range(len(constants.ENERGIES)):

            for j in range(len(constants.ANGLES_FOR_S_FLUXES)):

                file_out.write(
                    "{0:1.14f} {1:1.5f} {2:1.14e}\n".format(
                        constants.ENERGIES[i],
                        constants.ANGLES_FOR_S_FLUXES[j],
                        s_fluxes[i, j],
                    )
                )

        file_out.close()

        if constants.get_verbose() > 1:
            print("Surface fluxes written to " + file_name + ".")

    return s_fluxes


# Get surface fluxes


def load_s_fluxes_from_file(
    location="USStd",
    month=None,
    interaction_model="SIBYLL-2.3c",
    primary_model="GSF",
    atmosphere="CORSIKA",
    force=False,
):

    """
    Retrieve a surface fluxes matrix stored in data/surface_fluxes based on the input parameters.

    The function searches for a file name that matches the set location, month, interaction model, and primary model. If the file does not exist, prompt the user to run calc_s_fluxes().

    Parameters
    ----------
    location : str, optional (default: "USStd")
        The name of the location for which to calculate the surface fluxes. See the Tutorial or MCEq documentation for a list of options.

    month : str, optional (default: None)
        The month for which to calculate the surface fluxes. For US Standard Atmosphere, use None. For seasonal variations, use the month name.

    interaction_model: str, optional (default: "SIBYLL-2.3c")
        The hadronic interaction model to use in MCEq. See the Tutorial or MCEq documentation for a list of options.

    primary_model : {"GSF", "HG", "GH"}, optional (default: "GSF")
        The primary flux model to use in MCEq. Options:
        GSF = GlobalSplineFitBeta
        HG  = HillasGaisser2012
        GH  = GaisserHonda

    atmosphere : {"CORSIKA", "MSIS00"}, optional (default: "CORSIKA")
        The atmospheric model. For US Standard Atmosphere, use "CORSIKA". For seasonal variations, use "MSIS00".

    force : bool
        If True, force the calculation of a new surface fluxes matrix if required.

    Returns
    -------
    s_fluxes : NumPy ndarray
        A two-dimensional array containing the surface fluxes.
    """

    # Check values

    constants.check_constants(force=force)

    # Define a function to run if there is no surface fluxes file

    def no_file(force):

        # If the file does not exist, ask the user if they want to run MCEq to create it

        if not force:

            answer = input(
                "No surface flux matrix currently exists for these models. Would you like to create one (y/n)?: "
            )

        if force or answer.lower() == "y":

            s_fluxes_full = calc_s_fluxes(
                location,
                month,
                interaction_model,
                primary_model,
                atmosphere,
                force=force,
            )

            return s_fluxes_full

        else:

            print("Surface fluxes not calculated.")

            return None

    # Construct the file name based on the user's inputs

    if month is None:

        file_name = (
            constants.get_directory()
            + "/surface_fluxes/Surface_Fluxes_"
            + location
            + "_"
            + interaction_model
            + "_"
            + primary_model
            + ".txt"
        )

    else:

        file_name = (
            constants.get_directory()
            + "/surface_fluxes/Surface_Fluxes_"
            + location
            + "_"
            + month
            + "_"
            + interaction_model
            + "_"
            + primary_model
            + ".txt"
        )

    # Check if the file exists

    if os.path.exists(file_name):

        if constants.get_verbose() > 1:
            print(
                "Loading surface fluxes for "
                + location
                + " using "
                + interaction_model
                + " and "
                + primary_model
                + "."
            )

        # If the file exists, read the fluxes in from it

        file = open(file_name, "r")
        n_lines = len(file.read().splitlines())

        file.close()

        # Check that the file has the correct number of energies and zenith angles

        if n_lines == constants.len_ij:

            s_fluxes = np.reshape(
                np.loadtxt(file_name)[:, 2],
                (len(constants.ENERGIES), len(constants.ANGLES_FOR_S_FLUXES)),
            )

            if constants.get_verbose() > 1:
                print("Loaded surface fluxes.")

            return s_fluxes

        else:

            return no_file(force=force)

    else:

        return no_file(force=force)


# Read the energies and angles from a surface fluxes file


def print_s_fluxes_grids(file_name):

    """Return the surface energy grid and zenith angles in a surface fluxes file."""

    file_contents = np.loadtxt(
        constants.get_directory() + "/surface_fluxes/" + file_name
    )

    file_s_energies = np.unique(file_contents[:, 0])
    file_angles = np.unique(file_contents[:, 1])

    print("This file has " + str(len(file_s_energies)) + " surface energies:")
    print(file_s_energies)
    print("This file has " + str(len(file_angles)) + " zenith angles:")
    print(file_angles)
