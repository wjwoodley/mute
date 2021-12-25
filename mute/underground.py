##########################
##########################
###                    ###
###  MUTE              ###
###  William Woodley   ###
###  24 December 2021  ###
###                    ###
##########################
##########################

# Import packages

import os

import numpy as np
import scipy.integrate as scii

import mute.constants as constants
import mute.surface as surface
import mute.propagation as propagation

# Calculate underground fluxes


def calc_u_fluxes(
    angles=None,
    location="USStd",
    month=None,
    interaction_model="SIBYLL-2.3c",
    primary_model="GSF",
    atmosphere="CORSIKA",
    output=None,
    force=False,
):

    """
    Calculate underground fluxes using surface flux and survival probability matrices stored in /data based on the input parameters and the set global parameters.

    The function searches for the proper surface flux and survival probability files based on the set parameters using the mts.load_s_fluxes_from_file() and mtp.load_survival_probability_tensor_from_file() functions. If a file does not exist, the user is prompted for it to be calculated.

    Parameters
    ----------
    angles : array-like, optional (default: taken from consants.angles)
        An array of zenith angles in degrees to calculate the underground fluxes for.

    location : str, optional (default: "USStd")
        The name of the location for which to calculate the surface fluxes. See the Tutorial or MCEq documentation for a list of options.

    month : str, optional (default: None)
        The month for which to calculate the surface fluxes. For US Standard Atmosphere, use None. For seasonal variations, use the month name.

    interaction_model: str, optional (default: "SIBYLL-2.3c")
        The hadronic interaction model to use in MCEq. See the Tutorial or MCEq documentation for a list of options.

    primary_model : str in {"GSF", "HG", "GH", "ZS", "ZSP"} or tuple, optional (default: "GSF")
        The primary flux model to use in MCEq. Options:
        GSF = GlobalSplineFitBeta
        HG  = HillasGaisser2012 (H3a)
        GH  = GaisserHonda
        ZS  = Zatsepin-Sokolskaya (Default)
        ZSP = Zatsepin-Sokolskaya (PAMELA)
        Alternatively, this can be set with a tuple. For example: (pm.GaisserStanevTilav, "3-gen")

    atmosphere : {"CORSIKA", "MSIS00"}, optional (default: "CORSIKA")
        The atmospheric model. For US Standard Atmosphere, use "CORSIKA". For seasonal variations, use "MSIS00".

    output : bool, optional (default: taken from constants.get_output())
        If True, an output file will be created to store the results.

    force : bool
        If True, force the calculation of a new surface flux matrix and survival probability tensor if required.

    Returns
    -------
    u_fluxes : tuple of NumPy ndarray
        A two-element tuple. The first element is a two-dimensional array containing the underground fluxes. The second element is a two-dimensional array containing the underground fluxes calculated with an angle of zero degrees for use in calculation of the true vertical intensities.
    """

    # Import packages

    import scipy.interpolate as scii

    # Check values

    if output is None:
        output = constants.get_output()

    # Get the surface flux and survival prbability matrices
    # Getter functions do not take output in as an argument
    # If the calculation functions have to be run through the getter functions, the global output setting is used

    if constants.get_verbose() > 1:
        print("Calculating underground fluxes.")

    s_fluxes = surface.load_s_fluxes_from_file(
        location, month, interaction_model, primary_model, atmosphere, force=force
    )
    survival = propagation.load_survival_probability_tensor_from_file(force=force)

    if (s_fluxes is None or survival is None) and constants.get_verbose() > 1:

        print("Underground fluxes not calculated.")

        return

    # Interpolate the matrices

    interp_at_angles = np.linspace(
        np.min(constants.angles), np.max(constants.angles), 300
    )
    interp_s_fluxes = scii.interp1d(
        constants.ANGLES_FOR_S_FLUXES, s_fluxes, axis=1, kind="cubic"
    )(interp_at_angles)
    interp_s_fluxes = scii.RectBivariateSpline(
        constants.ENERGIES, interp_at_angles, interp_s_fluxes
    )(constants.ENERGIES, constants.angles)

    interp_survival = scii.interp1d(
        constants.SLANT_DEPTHS, survival, axis=1, kind="cubic"
    )(constants.slant_depths)

    # Reshape the matrices

    interp_s_fluxes = np.nan_to_num(
        np.reshape(interp_s_fluxes, (len(constants.ENERGIES), len(constants.angles)))
    )
    interp_survival = np.nan_to_num(
        np.reshape(
            interp_survival,
            (
                len(constants.ENERGIES),
                len(constants.slant_depths),
                len(constants.ENERGIES),
            ),
        )
    )

    # Calculate the underground fluxes for a flat overburden

    if constants.get_overburden() == "flat":

        # Initialise the underground flux matrices
        # First index  = Surface energy
        # Second index = Zenith angle
        # Third index  = Underground energy

        u_fluxes = np.zeros(
            (len(constants.ENERGIES), len(constants.angles), len(constants.ENERGIES))
        )
        u_fluxes_tr = np.zeros(
            (len(constants.ENERGIES), len(constants.angles), len(constants.ENERGIES))
        )

        # Calculate the underground fluxes

        for i in range(len(constants.ENERGIES)):

            for j in range(len(constants.angles)):

                for u in range(len(constants.ENERGIES)):

                    u_fluxes[i, j, u] = (
                        interp_survival[i, j, u]
                        * interp_s_fluxes[i, j]
                        * (constants.E_WIDTHS[i] / constants.E_WIDTHS[u])
                    )
                    u_fluxes_tr[i, j, u] = (
                        interp_survival[i, j, u]
                        * interp_s_fluxes[i, 0]
                        * (constants.E_WIDTHS[i] / constants.E_WIDTHS[u])
                    )

        # Sum over the surface energy grid axis
        # This reduces the number of axes down to two
        # First index = Zenith angle
        # Second index  = Underground energy

        u_fluxes = np.sum(u_fluxes, axis=0)
        u_fluxes_tr = np.sum(u_fluxes_tr, axis=0)

        # Set the angles to default values
        # If user has input singular angle, turn into an array

        if angles is None:
            angles = constants.angles

        angles = np.atleast_1d(angles)

        # If the user has not specified angles, return the whole matrix with no interpolation
        # Check to make sure both arrays have the same number of values, and that those values are equal at each index

        if len(angles) == len(constants.ANGLES) and np.allclose(
            angles, constants.ANGLES
        ):

            if constants.get_verbose() > 1:
                print("Finished calculating underground fluxes.")

            # Write the results to the file

            if output:

                constants.check_directory(
                    os.path.join(constants.get_directory(), "underground"), force=force
                )

                file_name = os.path.join(
                    constants.get_directory(),
                    "underground",
                    "{0}_Underground_Fluxes.txt".format(constants.get_lab()),
                )
                file_out = open(file_name, "w")

                for j in range(len(constants.ANGLES)):

                    for u in range(len(constants.ENERGIES)):

                        file_out.write(
                            "{0:1.5f} {1:1.14f} {2:1.14e}\n".format(
                                angles[j], constants.ENERGIES[u], u_fluxes[j, u]
                            )
                        )

                file_out.close()

                if constants.get_verbose() > 1:
                    print("Underground fluxes written to " + file_name + ".")

            return u_fluxes, u_fluxes_tr

        else:

            # Interpolate at the angles the user has requested

            interp_at_angles = np.linspace(
                np.min(constants.angles), np.max(constants.angles), 300
            )
            interp_u_fluxes = scii.interp1d(
                constants.angles, u_fluxes, axis=0, kind="cubic"
            )(interp_at_angles)
            interp_u_fluxes = scii.RectBivariateSpline(
                interp_at_angles, constants.ENERGIES, interp_u_fluxes
            )(angles, constants.ENERGIES)
            interp_u_fluxes_tr = scii.interp1d(
                constants.angles, u_fluxes_tr, axis=0, kind="cubic"
            )(interp_at_angles)
            interp_u_fluxes_tr = scii.RectBivariateSpline(
                interp_at_angles, constants.ENERGIES, interp_u_fluxes_tr
            )(angles, constants.ENERGIES)

            # Reshape into a matrix of zeroth dimension len(angles)

            interp_u_fluxes = np.nan_to_num(
                np.reshape(interp_u_fluxes, (len(angles), len(constants.ENERGIES)))
            )
            interp_u_fluxes_tr = np.nan_to_num(
                np.reshape(interp_u_fluxes_tr, (len(angles), len(constants.ENERGIES)))
            )

            # Set the global variables

            u_fluxes = interp_u_fluxes
            u_fluxes_tr = interp_u_fluxes_tr

            if constants.get_verbose() > 1:
                print("Finished calculating underground fluxes.")

            # Write the results to the file

            if output:

                constants.check_directory(
                    os.path.join(constants.get_directory(), "underground"), force=force
                )

                file_name = os.path.join(
                    constants.get_directory(),
                    "underground",
                    "{0}_Underground_Fluxes.txt".format(constants.get_lab()),
                )

                file_out = open(file_name, "w")

                for j in range(len(angles)):

                    for u in range(len(constants.ENERGIES)):

                        file_out.write(
                            "{0:1.5f} {1:1.14f} {2:1.14e}\n".format(
                                angles[j], constants.ENERGIES[u], u_fluxes[j, u]
                            )
                        )

                file_out.close()

                if constants.get_verbose() > 1:
                    print("Underground fluxes written to " + file_name + ".")

            return u_fluxes, u_fluxes_tr

    # Calculate the underground fluxes for a non-flat overburden

    else:

        raise NotImplementedError("Non-flat overburdens are not yet implemented.")


# Read the angles and energies from an underground fluxes file


def print_u_fluxes_grids(file_name):

    """Return the zenith angles and underground energies in an underground fluxes file."""

    file_contents = np.loadtxt(constants.get_directory() + "/underground/" + file_name)

    file_angles = np.unique(file_contents[:, 0])
    file_u_energies = np.unique(file_contents[:, 1])

    print("This file has " + str(len(file_angles)) + " zenith angles:")
    print(file_angles)
    print("This file has " + str(len(file_u_energies)) + " underground energies:")
    print(file_u_energies)

    return


# Calculate bins for use in Geant4


def print_geant4_bins(histtype):

    """Print underground energy or zenith angle histogram bins for use in a Geant4 macro with gps."""

    u_fluxes_sum = np.sum(u_fluxes)

    if histtype == "theta":

        u_fluxes_theta_sum = np.sum(u_fluxes, axis=1) / u_fluxes_sum

        print("/gps/hist/type theta")
        print("/gps/hist/point {0:1.5e} {1:1.5e}".format(constants.ANGLES[0], 0))

        for j in range(len(constants.ANGLES)):

            print(
                "/gps/hist/point {0:1.5e} {1:1.5e}".format(
                    constants.ANGLES[j], u_fluxes_theta_sum[j]
                )
            )

    elif histtype == "energy":

        u_fluxes_energy_sum = np.sum(u_fluxes, axis=0) / u_fluxes_sum

        print("/gps/hist/type energy")
        print("/gps/hist/point {0:1.5e} {1:1.5e}".format(constants.E_BINS[0], 0))

        for u in range(len(constants.ENERGIES)):

            print(
                "/gps/hist/point {0:1.5e} {1:1.5e}".format(
                    constants.E_BINS[u + 1], u_fluxes_energy_sum[u]
                )
            )

    else:

        raise ValueError('histtype must be "energy" or "theta".')


# Calculate the underground intensities


def calc_u_intensities(
    angles=None,
    location="USStd",
    month=None,
    interaction_model="SIBYLL-2.3c",
    primary_model="GSF",
    atmosphere="CORSIKA",
    output=None,
    force=False,
):

    """
    Calculate underground intensities.

    Parameters
    ----------
    angles : array-like, optional (default: taken from consants.angles)
        An array of zenith angles in degrees to calculate the underground intensities for.

    location : str, optional (default: "USStd")
        The name of the location for which to calculate the surface fluxes. See the Tutorial or MCEq documentation for a list of options.

    month : str, optional (default: None)
        The month for which to calculate the surface fluxes. For US Standard Atmosphere, use None. For seasonal variations, use the month name.

    interaction_model: str, optional (default: "SIBYLL-2.3c")
        The hadronic interaction model to use in MCEq. See the Tutorial or MCEq documentation for a list of options.

    primary_model : str in {"GSF", "HG", "GH", "ZS", "ZSP"} or tuple, optional (default: "GSF")
        The primary flux model to use in MCEq. Options:
        GSF = GlobalSplineFitBeta
        HG  = HillasGaisser2012 (H3a)
        GH  = GaisserHonda
        ZS  = Zatsepin-Sokolskaya (Default)
        ZSP = Zatsepin-Sokolskaya (PAMELA)
        Alternatively, this can be set with a tuple. For example: (pm.GaisserStanevTilav, "3-gen")

    atmosphere : {"CORSIKA", "MSIS00"}, optional (default: "CORSIKA")
        The atmospheric model. For US Standard Atmosphere, use "CORSIKA". For seasonal variations, use "MSIS00".

    output : bool, optional (default: taken from constants.get_output())
        If True, an output file will be created to store the results.

    force : bool
        If True, force the calculation of new matrices if required.

    Returns
    -------
    u_intensities : NumPy array
        An array containing the underground intensities.
    """

    # Check values

    if output is None:
        output = constants.get_output()

    # Calculate the underground intensities for a flat overburden

    if constants.get_overburden() == "flat":

        # Set the angles to default values
        # If the user has input a singular angle, turn it into an array

        if angles is None:
            angles = constants.angles

        angles = np.atleast_1d(angles)

        # Calculate the underground fluxes

        u_fluxes = calc_u_fluxes(
            angles,
            location,
            month,
            interaction_model,
            primary_model,
            atmosphere,
            output,
            force=force,
        )

        if u_fluxes is None and constants.get_verbose() > 1:

            print("Underground intensities not calculated.")

            return

        u_fluxes = u_fluxes[0]

        # Initialise the underground intensities array

        u_intensities = np.zeros(len(angles))

        # Calculate the underground intensities

        if constants.get_verbose() > 1:
            print("Calculating underground intensities.")

        for j in range(len(angles)):

            u_intensities[j] = scii.simpson(u_fluxes[j, :], constants.ENERGIES)

        if constants.get_verbose() > 1:
            print("Finished calculating underground intensities.")

        # Write the results to a file

        if output:

            constants.check_directory(
                os.path.join(constants.get_directory(), "underground"), force=force
            )

            file_name = os.path.join(
                constants.get_directory(),
                "underground",
                "{0}_Underground_Intensities.txt".format(constants.get_lab()),
            )

            file_out = open(file_name, "w")

            for j in range(len(angles)):

                file_out.write(
                    "{0:1.5f} {1:1.14e}\n".format(angles[j], u_intensities[j])
                )

            if constants.get_verbose() > 1:
                print("Underground intensities written to " + file_name + ".")

        return u_intensities

    # Calculate the underground intensities for a non-flat overburden

    else:

        raise NotImplementedError("Non-flat overburdens are not yet implemented.")


# Read the angles from an underground intensities file


def print_u_intensities_grid(file_name):

    """Return the zenith angles in an underground intensities file."""

    file_contents = np.loadtxt(constants.get_directory() + "/underground/" + file_name)
    file_angles = np.unique(file_contents[:, 0])

    print("This file has " + str(len(file_angles)) + " zenith angles:")
    print(file_angles)

    return


# Calculate true vertical underground intensities


def calc_u_intensities_tr(
    depths=None,
    location="USStd",
    month=None,
    interaction_model="SIBYLL-2.3c",
    primary_model="GSF",
    atmosphere="CORSIKA",
    output=None,
    force=False,
):

    """
    Calculate true vertical underground intensities.

    Parameters
    ----------
    depths : array-like, optional (default: taken from consants.slant_depths)
        An array of slant depths in [km.w.e.] to calculate the true vertical underground intensities for.

    location : str, optional (default: "USStd")
        The name of the location for which to calculate the surface fluxes. See the Tutorial or MCEq documentation for a list of options.

    month : str, optional (default: None)
        The month for which to calculate the surface fluxes. For US Standard Atmosphere, use None. For seasonal variations, use the month name.

    interaction_model: str, optional (default: "SIBYLL-2.3c")
        The hadronic interaction model to use in MCEq. See the Tutorial or MCEq documentation for a list of options.

    primary_model : str in {"GSF", "HG", "GH", "ZS", "ZSP"} or tuple, optional (default: "GSF")
        The primary flux model to use in MCEq. Options:
        GSF = GlobalSplineFitBeta
        HG  = HillasGaisser2012 (H3a)
        GH  = GaisserHonda
        ZS  = Zatsepin-Sokolskaya (Default)
        ZSP = Zatsepin-Sokolskaya (PAMELA)
        Alternatively, this can be set with a tuple. For example: (pm.GaisserStanevTilav, "3-gen")

    atmosphere : {"CORSIKA", "MSIS00"}, optional (default: "CORSIKA")
        The atmospheric model. For US Standard Atmosphere, use "CORSIKA". For seasonal variations, use "MSIS00".

    output : bool, optional (default: taken from constants.get_output())
        If True, an output file will be created to store the results.

    force : bool
        If True, force the calculation of new matrices if required.

    Returns
    -------
    u_intensities_tr : NumPy array
        An array containing the true vertical underground intensities.
    """

    # Check values

    if output is None:
        output = constants.get_output()
    if depths is None:
        depths = constants.slant_depths

    depths = np.atleast_1d(depths)

    assert all(
        depths >= np.min(constants.slant_depths)
    ), "Cannot calculate for depths lower than the set vertical depth."
    assert all(depths <= 12), "Cannot calculate for depths greater than 12 km.w.e."

    # Calculate the true vertical underground intensities for a flat overburden

    if constants.get_overburden() == "flat":

        # Calculate the underground fluxes

        u_fluxes = calc_u_fluxes(
            np.degrees(np.arccos(constants.get_vertical_depth() / depths)),
            location,
            month,
            interaction_model,
            primary_model,
            atmosphere,
            output,
            force=force,
        )

        if u_fluxes is None and constants.get_verbose() > 1:

            print("True vertical underground intensities not calculated.")

            return

        u_fluxes_tr = u_fluxes[1]

        # Initialise the underground intensity arrays

        u_intensities_tr = np.zeros(len(depths))

        # Calculate the true vertical underground intensities

        if constants.get_verbose() > 1:
            print("Calculating true vertical underground intensities.")

        for x in range(len(depths)):

            u_intensities_tr[x] = scii.simpson(u_fluxes_tr[x, :], constants.ENERGIES)

        if constants.get_verbose() > 1:
            print("Finished calculating true vertical underground intensities.")

        # Write the results to a file

        if output:

            constants.check_directory(
                os.path.join(constants.get_directory(), "underground"), force=force
            )

            file_name = os.path.join(
                constants.get_directory(),
                "underground",
                "{0}_Underground_Intensities_TR.txt".format(constants.get_lab()),
            )

            file_out = open(file_name, "w")

            for x in range(len(depths)):

                file_out.write(
                    "{0:1.5f} {1:1.14e}\n".format(depths[x], u_intensities_tr[x])
                )

            if constants.get_verbose() > 1:

                print(
                    "True vertical underground intensities written to "
                    + file_name
                    + "."
                )

        return u_intensities_tr

    # Calculate the true vertical underground intensities for a non-flat overburden

    else:

        raise NotImplementedError("Non-flat overburdens are not yet implemented.")


# Read the depths from a true vertical underground intensities file


def print_u_intensities_tr_grid(file_name):

    """Return the slant depths in a true vertical underground intensities file."""

    file_contents = np.loadtxt(constants.get_directory() + "/underground/" + file_name)
    file_depths = np.unique(file_contents[:, 0])

    print("This file has " + str(len(file_depths)) + " slant depths:")
    print(file_depths)

    return


# Calculate vertical-equivalent underground intensities


def calc_u_intensities_eq(
    angles=None,
    location="USStd",
    month=None,
    interaction_model="SIBYLL-2.3c",
    primary_model="GSF",
    atmosphere="CORSIKA",
    output=None,
    force=False,
):

    """
    Calculate vertical-equivalent underground intensities.

    Parameters
    ----------
    angles : array-like, optional (default: taken from consants.angles)
        An array of zenith angles in degrees to calculate the vertical-equivalent underground intensities for.

    location : str, optional (default: "USStd")
        The name of the location for which to calculate the surface fluxes. See the Tutorial or MCEq documentation for a list of options.

    month : str, optional (default: None)
        The month for which to calculate the surface fluxes. For US Standard Atmosphere, use None. For seasonal variations, use the month name.

    interaction_model: str, optional (default: "SIBYLL-2.3c")
        The hadronic interaction model to use in MCEq. See the Tutorial or MCEq documentation for a list of options.

    primary_model : str in {"GSF", "HG", "GH", "ZS", "ZSP"} or tuple, optional (default: "GSF")
        The primary flux model to use in MCEq. Options:
        GSF = GlobalSplineFitBeta
        HG  = HillasGaisser2012 (H3a)
        GH  = GaisserHonda
        ZS  = Zatsepin-Sokolskaya (Default)
        ZSP = Zatsepin-Sokolskaya (PAMELA)
        Alternatively, this can be set with a tuple. For example: (pm.GaisserStanevTilav, "3-gen")

    atmosphere : {"CORSIKA", "MSIS00"}, optional (default: "CORSIKA")
        The atmospheric model. For US Standard Atmosphere, use "CORSIKA". For seasonal variations, use "MSIS00".

    output : bool, optional (default: taken from constants.get_output())
        If True, an output file will be created to store the results.

    force : bool
        If True, force the calculation of new matrices if required.

    Returns
    -------
    u_intensities_eq : NumPy array
        An array containing the vertical-equivalent underground intensities.
    """

    # Check values

    if output is None:
        output = constants.get_output()

    # Calculate the vertical-equivalent underground intensities for a flat overburden

    if constants.get_overburden() == "flat":

        # Set the angles to the default values
        # If the user has input a singular angle, turn into an array

        if angles is None:
            angles = constants.angles

        angles = np.atleast_1d(angles)

        # Calculate the underground fluxes

        u_fluxes = calc_u_fluxes(
            angles,
            location,
            month,
            interaction_model,
            primary_model,
            atmosphere,
            output,
            force=force,
        )

        if u_fluxes is None and constants.get_verbose() > 1:

            print("Vertical-equivalent underground intensities not calculated.")

            return

        u_fluxes = u_fluxes[0]

        # Initialise the vertical-equivalent underground intensity array

        u_intensities_eq = np.zeros(len(angles))

        # Calculate the vertical-equivalent underground intensities

        if constants.get_verbose() > 1:

            print("Calculating vertical-equivalent underground intensities.")

        for j in range(len(angles)):

            u_intensities_eq[j] = scii.simpson(
                u_fluxes[j, :], constants.ENERGIES
            ) * np.cos(np.radians(angles[j]))

        if constants.get_verbose() > 1:

            print("Finished calculating vertical-equivalent underground intensities.")

        # Write the results to a file

        if output:

            constants.check_directory(
                os.path.join(constants.get_directory(), "underground"), force=force
            )

            file_name = os.path.join(
                constants.get_directory(),
                "underground",
                "{0}_Underground_Intensities_EQ.txt".format(constants.get_lab()),
            )

            file_out = open(file_name, "w")

            for j in range(len(angles)):

                file_out.write(
                    "{0:1.5f} {1:1.14e}\n".format(angles[j], u_intensities_eq[j])
                )

            if constants.get_verbose() > 1:

                print(
                    "Vertical-equivalent underground intensities written to "
                    + file_name
                    + "."
                )

        return u_intensities_eq

    # Calculate the vertical-equivalent underground intensities for a non-flat overburden

    else:

        raise NotImplementedError("Non-flat overburdens are not yet implemented.")


# Read the depths from a vertical-equivalent underground intensities file


def print_u_intensities_eq_grid(file_name):

    """Return the zenith angles in a vertical-equivalent underground intensities file."""

    file_contents = np.loadtxt(constants.get_directory() + "/underground/" + file_name)
    file_angles = np.unique(file_contents[:, 0])

    print("This file has " + str(len(file_angles)) + " zenith angles:")
    print(file_angles)

    return


# Calculate the total fluxes


def calc_u_tot_flux(
    angles=None,
    location="USStd",
    month=None,
    interaction_model="SIBYLL-2.3c",
    primary_model="GSF",
    atmosphere="CORSIKA",
    force=False,
):

    """
    Calculate the total underground flux by calling calc_u_intensities().

    Parameters
    ----------
    angles : array-like, optional (default: taken from consants.angles)
        An array of zenith angles in degrees to calculate the underground intensities for.

    location : str, optional (default: "USStd")
        The name of the location for which to calculate the surface fluxes. See the Tutorial or MCEq documentation for a list of options.

    month : str, optional (default: None)
        The month for which to calculate the surface fluxes. For US Standard Atmosphere, use None. For seasonal variations, use the month name.

    interaction_model: str, optional (default: "SIBYLL-2.3c")
        The hadronic interaction model to use in MCEq. See the Tutorial or MCEq documentation for a list of options.

    primary_model : str in {"GSF", "HG", "GH", "ZS", "ZSP"} or tuple, optional (default: "GSF")
        The primary flux model to use in MCEq. Options:
        GSF = GlobalSplineFitBeta
        HG  = HillasGaisser2012 (H3a)
        GH  = GaisserHonda
        ZS  = Zatsepin-Sokolskaya (Default)
        ZSP = Zatsepin-Sokolskaya (PAMELA)
        Alternatively, this can be set with a tuple. For example: (pm.GaisserStanevTilav, "3-gen")

    atmosphere : {"CORSIKA", "MSIS00"}, optional (default: "CORSIKA")
        The atmospheric model. For US Standard Atmosphere, use "CORSIKA". For seasonal variations, use "MSIS00".

    force : bool
        If True, force the calculation of new matrices if required.

    Returns
    -------
    u_tot_flux : float
        The total underground flux in units of [cm^(-2) s^(-1)].
    """

    # Calculate the total underground flux for a flat overburden

    if constants.get_overburden() == "flat":

        # Set the angles to default values

        if angles is None:

            angles = constants.angles

        # Calculate the underground intensities

        u_intensities = calc_u_intensities(
            angles,
            location,
            month,
            interaction_model,
            primary_model,
            atmosphere,
            force=force,
        )

        if u_intensities is None and constants.get_verbose() > 1:

            print("Total underground fluxes not calculated.")

            return

        # Calculate the total underground flux
        # Because angles goes (0..89), cos(angles) goes (1..0)
        # cos(angles) is decreasing, but scii.simpson() wants an increasing array
        # Therefore, integrate backwards, using [::-1] on the integrand and steps
        # Otherwise, the answer will be negative

        u_tot_flux = (
            2
            * np.pi
            * scii.simpson(u_intensities[::-1], np.cos(np.radians(angles[::-1])))
        )

    # Calculate the total underground flux for a non-flat overburden

    else:

        raise NotImplementedError("Non-flat overburdens are not yet implemented.")

    return u_tot_flux
