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
import scipy.interpolate as sciint

import mute.constants as constants
import mute.surface as surface
import mute.propagation as propagation

# Perform the convolution for underground fluxes


def _do_convolution(s_fluxes, survival_probability_tensor):

    """This function performs the convolution to calculate underground fluxes given a surface flux matrix and survival probability tensor."""

    # Check values

    constants._check_constants()

    s_fluxes = np.atleast_2d(s_fluxes)

    assert s_fluxes.shape[0] == len(
        constants.ENERGIES
    ), f"The surface flux matrix must use the default MUTE energy grid, which has length {len(constants.ENERGIES)}, not {s_fluxes.shape[0]}."

    # In order for the dot product in the convolution to work, move the axes
    # Move axis 0 (surface energies) to position 2
    # Original: (i, x, u)
    # Moved:    (x, u, i)

    survival_probability_tensor = np.moveaxis(survival_probability_tensor, 0, 2)

    # Construct a matrix of energy bin widths

    widths = np.repeat(constants._E_WIDTHS, s_fluxes.shape[1]).reshape(
        (len(constants.ENERGIES), s_fluxes.shape[1])
    )

    # Calculate the underground flux tensor

    u_fluxes_from_convolution = (
        survival_probability_tensor.dot((s_fluxes.T * constants._E_WIDTHS).T) / widths
    )

    return u_fluxes_from_convolution


# Calculate underground fluxes


def calc_u_fluxes(
    s_fluxes=None,
    survival_probability_tensor=None,
    full_tensor=False,
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
    Calculate underground fluxes in units of [(cm^2 s sr MeV)^-1] for default surface energy grid and zenith angles.

    The default surface energy grid is given by constants.ENERGIES, and the default zenith angles are given by constants._ANGLES.

    Parameters
    ----------
    s_fluxes : NumPy ndarray, optional (default: taken from surface.load_s_fluxes_from_file())
        A surface flux matrix of shape (91, 20).

    survival_probability_tensor : NumPy ndarray, optional (default: taken from propagation.load_survival_probability_tensor_from_file())
        A survival probability tensor of shape (91, 28, 91).

    full_tensor : bool, optional (default: False)
        If True, the full tensor of shape (28, 91, 20) will be returned. Otherwise, if the overburden type is flat, a two-dimensional matrix of shape (91, 28) will be returned.

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

    force : bool
        If True, force the calculation of new matrices and the creation of new directories if required.

    Returns
    -------
    u_fluxes : NumPy ndarray
        A multi-dimensional array containing the underground fluxes. For flat overburdens, the shape will be (91, 28) for the default zenith angles, and the fluxes will be in units of [(cm^2 s sr MeV)^-1]. For mountain overburdens (or if full_tensor is True), the shape will be (28, 91, 20), and the fluxes will be in units of [(cm^2 s sr MeV^2)^-1].
    """

    # Check values

    constants._check_constants()

    if output is None:
        output = constants.get_output()

    # Get the surface flux matrix and survival probability tensor
    # If the calculation functions have to be run through the getter functions, the global output setting is used

    if constants.get_verbose() > 1:
        print("Calculating underground fluxes.")

    if s_fluxes is None:

        s_fluxes = surface.load_s_fluxes_from_file(
            primary_model, interaction_model, atmosphere, location, month, force=force
        )

    if survival_probability_tensor is None:

        survival_probability_tensor = (
            propagation.load_survival_probability_tensor_from_file(force=force)
        )

    # Check that the surface flux matrix and survival probability tensor have been loaded properly

    if s_fluxes is None or survival_probability_tensor is None:

        raise Exception(
            "Underground fluxes not calculated. The surface flux matrix or survival probability tensor was not provided or loaded correctly."
        )

    assert (
        s_fluxes.shape[0] == survival_probability_tensor.shape[2]
    ), f"The surface flux matrix and survival probability tensor must both use the default MUTE energy grid, which has length {len(constants.ENERGIES)}, not {s_fluxes.shape[0]} or {survival_probability_tensor.shape[2]}."

    # Perform the convolution

    u_fluxes_from_convolution = _do_convolution(s_fluxes, survival_probability_tensor)

    # Return an interpolated matrix for flat overburdens

    if constants.get_overburden() == "flat" and not full_tensor:

        # Interpolate the u_fluxes tensor to a constant grid
        # This allows the np.diagonal() function to be used to calculate the matrix of relevant fluxes from the full tensor

        interp_u_fluxes_from_convolution = sciint.interp1d(
            constants.ANGLES_FOR_S_FLUXES,
            u_fluxes_from_convolution,
            axis=2,
            kind="linear",
            fill_value=0,
        )(constants._ANGLES)

        # Extract the (u_energies, angles) diagonal from the convolution tensor

        u_fluxes = np.diagonal(interp_u_fluxes_from_convolution, axis1=0, axis2=2)

        # Slice the angles corresponding to constants.angles

        u_fluxes = u_fluxes[:, (len(constants._ANGLES) - len(constants.angles)) :]

        if constants.get_verbose() > 1:
            print("Finished calculating underground fluxes.")

        # Write the results to the file

        if output:

            constants._check_directory(
                os.path.join(constants.get_directory(), "underground"), force=force
            )

            if file_name == "" or not isinstance(file_name, str):

                file_name = os.path.join(
                    constants.get_directory(),
                    "underground",
                    "{0}_Underground_Fluxes.txt".format(constants.get_lab()),
                )

            file_out = open(file_name, "w")

            for u in range(len(constants.ENERGIES)):

                for j in range(len(constants.angles)):

                    file_out.write(
                        "{0:1.14f} {1:1.5f} {2:1.14e}\n".format(
                            constants.ENERGIES[u], constants.angles[j], u_fluxes[u, j]
                        )
                    )

            file_out.close()

            if constants.get_verbose() > 1:
                print(f"Underground fluxes written to {file_name}.")

        return u_fluxes

    # Return the whole tensor for mountains

    elif constants.get_overburden() == "mountain" or full_tensor:

        u_fluxes = np.squeeze(u_fluxes_from_convolution)

        if constants.get_verbose() > 1:
            print("Finished calculating underground fluxes.")

        # Write the results to the file

        if output:

            constants._check_directory(
                os.path.join(constants.get_directory(), "underground"), force=force
            )

            if file_name == "" or not isinstance(file_name, str):

                file_name = os.path.join(
                    constants.get_directory(),
                    "underground",
                    "{0}_Underground_Fluxes.txt".format(constants.get_lab()),
                )

            file_out = open(file_name, "w")

            for x in range(len(constants._SLANT_DEPTHS)):

                for u in range(len(constants.ENERGIES)):

                    for j in range(len(constants.ANGLES_FOR_S_FLUXES)):

                        file_out.write(
                            "{0:1.5f} {1:1.14f} {2:1.5f} {3:1.14e}\n".format(
                                constants._SLANT_DEPTHS[x],
                                constants.ENERGIES[u],
                                constants.ANGLES_FOR_S_FLUXES[j],
                                u_fluxes[x, u, j],
                            )
                        )

            file_out.close()

            if constants.get_verbose() > 1:
                print(f"Underground fluxes written to {file_name}.")

        return u_fluxes_from_convolution

    else:

        raise NotImplementedError(
            'Overburdens of type {0} are not available. The only options are "flat" and "mountain".'.format(
                constants.get_overburden()
            )
        )


# Create the interpolator to calculate double-differential underground intensities


def _create_interpolator(u_intensities):

    """This function creates the interpolator object to use when computing double-differential underground intensities."""

    global _interpolator

    _interpolator = sciint.interp2d(
        constants._SLANT_DEPTHS,
        constants.ANGLES_FOR_S_FLUXES,
        np.log(u_intensities.T),
        kind="linear",
        fill_value=-np.inf,
    )

    return _interpolator


# Calculate the underground intensities


def calc_u_intensities(method, output=None, file_name="", force=False, **kwargs):

    """
    Calculate underground intensities in units of [(cm^2 s sr)^-1].

    Parameters
    ----------
    method : str
        The type of underground intensities to calculate. Options:
        sd = Single-differential underground intensities for flat overburdens
        eq = Vertical-equivalent underground intensities for flat overburdens
        tr = True vertical underground intensities for flat overburdens
        dd = Double-differential underground intensities for mountain overburdens

    output : bool, optional (default: taken from constants.get_output())
        If True, an output file will be created to store the results.

    file_name : str, optional (default: constructed from input parameters)
        The name of the file in which to store the results. If output is False or None, this is ignored.

    force : bool
        If True, force the calculation of new matrices and the creation of new directories if required.

    Other Parameters
    -----------------
    u_fluxes : NumPy ndarray, optional (default: calcaulated from default surface fluxes and survival probabilities)
        An underground flux tensor of shape (28, 91, 20).

    s_fluxes : NumPy ndarray, optional (default: taken from surface.load_s_fluxes_from_file())
        A surface flux matrix of shape (91, 20).

    survival_probability_tensor : NumPy ndarray, optional (default: taken from propagation.load_survival_probability_tensor_from_file())
        A survival probability tensor of shape (91, 28, 91).

    angles : array-like, optional (default: taken from consants.angles)
        An array of zenith angles in [degrees] to calculate the underground intensities for. These are only used if method is "sd" or "eq".

    depths : array-like, optional (default: taken from consants.slant_depths)
        An array of slant depths in [km.w.e.] to calculate the underground intensities for. These are only used if method is "tr".

    E_th : float, optional (default: 0)
        An energy threshold in [MeV].

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

    Returns
    -------
    u_intensities : NumPy ndarray
        An array containing the underground intensities. For flat overburdens (if method is "sd", "eq", or "tr"), the length will be that of angles, and the intensities will be in units of [(cm^2 s sr)^-1]. For mountain overburdens (if method is "dd"), the shape will be (len(constants.mountain.zenith), len(constants.mountain.azimuthal)), and the intensities will be in units of [(cm^2 s sr km.w.e.)^-1].
    """

    # Check values

    constants._check_constants()

    method = method.lower()

    assert method in [
        "sd",
        "eq",
        "tr",
        "dd",
    ], '"{0}" is not a valid method. Use "sd", "tr", "eq", or "dd". See the Tutorial at {1} for an explanation.'.format(
        method,
        "https://github.com/wjwoodley/mute/blob/main/docs/Tutorial.md#calculating-underground-intensities",
    )

    if output is None:
        output = constants.get_output()

    # Pop the keyword arguments

    u_fluxes = kwargs.pop("u_fluxes", None)
    s_fluxes = kwargs.pop("s_fluxes", None)
    survival_probability_tensor = kwargs.pop("survival_probability_tensor", None)
    angles = kwargs.pop("angles", constants.angles)
    depths = kwargs.pop("depths", constants.slant_depths)
    E_th = kwargs.pop("E_th", 0)
    primary_model = kwargs.pop("primary_model", "GSF")
    interaction_model = kwargs.pop("interaction_model", "SIBYLL-2.3c")
    atmosphere = kwargs.pop("atmosphere", "CORSIKA")
    location = kwargs.pop("location", "USStd")
    month = kwargs.pop("month", None)

    if kwargs:

        from warnings import warn as warning

        warning("{0} arguments not used.".format(kwargs))

    # Ensure the arguments are in the proper form
    # Turn the energy threshold value in [MeV] into an index in order to slice the u_fluxes matrix and the energy grid properly when integrating

    angles = np.atleast_1d(angles)
    depths = np.atleast_1d(depths)
    E_th_i = np.argmin(np.abs(constants.ENERGIES - E_th))

    if E_th_i != 0 and constants.get_verbose() > 1:

        print(
            "Setting energy threshold to {0:1.3f} MeV (the closest value in the energy grid to {1} MeV).".format(
                constants.ENERGIES[E_th_i], E_th
            )
        )

    # Calculate the underground fluxes
    # Get the surface flux matrix and survival probability tensor
    # If the calculation functions have to be run through the getter functions, the global output setting is used

    if u_fluxes is None:

        if constants.get_verbose() > 1:
            print("Calculating underground fluxes.")

        if s_fluxes is None:

            s_fluxes = surface.load_s_fluxes_from_file(
                primary_model,
                interaction_model,
                atmosphere,
                location,
                month,
                force=force,
            )

        if survival_probability_tensor is None:

            survival_probability_tensor = (
                propagation.load_survival_probability_tensor_from_file(force=force)
            )

        # Check that the surface flux matrix and survival probability tensor have been loaded properly

        if s_fluxes is None or survival_probability_tensor is None:

            raise Exception(
                "Underground intensities not calculated. The surface flux matrix or survival probability tensor was not provided or loaded correctly."
            )

        assert (
            s_fluxes.shape[0] == survival_probability_tensor.shape[2]
        ), f"The surface flux matrix and survival probability tensor must both use the default MUTE energy grid, which has length {len(constants.ENERGIES)}, not {s_fluxes.shape[0]} or {survival_probability_tensor.shape[2]}."

        # Perform the convolution

        u_fluxes = _do_convolution(s_fluxes, survival_probability_tensor)

        if constants.get_verbose() > 1:
            print("Finished calculating underground fluxes.")

    # Check that the underground flux tensor has the correct shape

    assert u_fluxes.shape == (
        28,
        91,
        20,
    ), f"The underground flux tensor does not have the correct shape. The shape must be (28, 91, 20), not {u_fluxes.shape}. Set full_tensor to True in the underground.calc_u_fluxes() function in order to return a tensor of shape (28, 91, 20). Alternatively, try passing a surface flux matrix and / or a survival probability tensor in as parameters."

    # Calculate the intensities according to the method type
    # sd = Single-differential intensities
    # tr = True vertical intensities
    # eq = Vertical-equivalent intensities
    # dd = Double-differential intensities (for mountains)

    # Calculate the single-differential standard underground intensities or vertical-equivalent underground intensities

    if method == "sd" or method == "eq":

        return _calc_u_intensities_sd(
            method,
            u_fluxes,
            angles,
            E_th_i,
            primary_model,
            interaction_model,
            atmosphere,
            location,
            month,
            output,
            file_name,
            force,
        )

    # Calculate the true vertical underground intensities

    elif method == "tr":

        return _calc_u_intensities_tr(
            u_fluxes,
            depths,
            E_th_i,
            primary_model,
            interaction_model,
            atmosphere,
            location,
            month,
            output,
            file_name,
            force,
        )

    # Calculate the double-differential underground intensities

    elif method == "dd":

        return _calc_u_intensities_dd(
            u_fluxes,
            E_th_i,
            primary_model,
            interaction_model,
            atmosphere,
            location,
            month,
            output,
            file_name,
            force,
        )

    else:

        raise ValueError(
            '"{0}" is not a valid method. Use "sd", "tr", "eq", or "dd". See the Tutorial at {1} for an explanation.'.format(
                method,
                "https://github.com/wjwoodley/mute/blob/main/docs/Tutorial.md#calculating-underground-intensities",
            )
        )


# Calculate single-differential underground intensities


def _calc_u_intensities_sd(
    method,
    u_fluxes,
    angles,
    E_th_i,
    primary_model,
    interaction_model,
    atmosphere,
    location,
    month,
    output,
    file_name,
    force,
):

    """This function calculates single-differential or vertical-equivalent underground intensities for flat overburdens."""

    # Check values

    constants._check_constants()

    if constants.get_overburden() != "flat":

        raise ValueError(
            "Single-differential underground intensities can only be calculated for flat overburdens. See the Tutorial at {0} for an explanation.".format(
                "https://github.com/wjwoodley/mute/blob/main/docs/Tutorial.md#calculating-underground-intensities"
            )
        )

    assert method in ["sd", "eq"]

    if u_fluxes is None:

        raise Exception(
            "Underground intensities not calculated. The underground flux matrix was not calculated correctly."
        )

    # Calculate the underground intensities
    # Zeroth index = Slant depth
    # First index  = Zenith angle (from constants.ANGLES_FOR_S_FLUXES)

    u_intensities = np.zeros(
        (len(constants._SLANT_DEPTHS), len(constants.ANGLES_FOR_S_FLUXES))
    )

    if constants.get_verbose() > 1:
        print("Calculating underground intensities.")

    # Whether provided or calculated, the u_fluxes matrix has to be of the shape (28, 91, 20)
    # Otherwise the cos(theta) approximation results in incorrect intensities

    for x in range(len(constants._SLANT_DEPTHS)):

        for j in range(len(constants.ANGLES_FOR_S_FLUXES)):

            u_intensities[x, j] = scii.simpson(
                u_fluxes[x, E_th_i:, j], constants.ENERGIES[E_th_i:]
            )

    # Interpolate intensities matrix to consistent grid given by set vertical depth
    # This allows the np.diagonal() function to be used to calculate the array of relevant intensities from the full matrix
    # Define a meshgrid for grid depths and angles
    # Define a meshgrid for adjusted depths and angles

    gd, ga = np.meshgrid(constants._SLANT_DEPTHS, constants.ANGLES_FOR_S_FLUXES)
    ad, aa = np.meshgrid(constants.slant_depths, constants.angles)

    u_intensities = np.nan_to_num(
        np.exp(
            sciint.griddata(
                (gd.flatten(), ga.flatten()),
                np.log(u_intensities.T.flatten()),
                (ad, aa),
            )
        )
    )
    u_intensities = np.diagonal(u_intensities)

    # Interpolate intensities to user-input angles
    # If the user has not specified angles, return the whole matrix with no interpolation
    # Check to make sure both arrays have the same number of values, and that those values are equal at each index

    if not (
        len(angles) == len(constants.angles) and np.allclose(angles, constants.angles)
    ):

        u_intensities = np.interp(angles, constants.angles, u_intensities)

    if method == "eq":

        u_intensities = np.cos(np.radians(angles)) * u_intensities

    if constants.get_verbose() > 1:
        print("Finished calculating underground intensities.")

    # Write the results to the file

    if output:

        constants._check_directory(
            os.path.join(constants.get_directory(), "underground"), force=force
        )

        if file_name == "" or not isinstance(file_name, str):

            file_name = os.path.join(
                constants.get_directory(),
                "underground",
                "{0}_Underground_Intensities_{1}.txt".format(
                    constants.get_lab(), method
                ),
            )

        file_out = open(file_name, "w")

        for j in range(len(angles)):

            file_out.write("{0:1.5f} {1:1.14e}\n".format(angles[j], u_intensities[j]))

        file_out.close()

        if constants.get_verbose() > 1:
            print(f"Underground intensities written to {file_name}.")

    return u_intensities


# Calculate true vertical underground intensities


def _calc_u_intensities_tr(
    u_fluxes,
    depths,
    E_th_i,
    primary_model,
    interaction_model,
    atmosphere,
    location,
    month,
    output,
    file_name,
    force,
):

    """This function calculates true vertical underground intensities for flat overburdens."""

    # Check values

    constants._check_constants()

    if constants.get_overburden() != "flat":

        raise ValueError(
            "True vertical underground intensities can only be calculated for flat overburdens. See the Tutorial at {0} for an explanation.".format(
                "https://github.com/wjwoodley/mute/blob/main/docs/Tutorial.md#calculating-underground-intensities"
            )
        )

    if u_fluxes is None:

        raise Exception(
            "Underground intensities not calculated. The underground flux matrix was not calculated correctly."
        )

    # Calculate the underground intensities

    if constants.get_verbose() > 1:
        print("Calculating true vertical underground intensities.")

    u_fluxes_tr = u_fluxes[:, :, 0]
    u_intensities_tr = [
        scii.simpson(u_fluxes_tr[x, E_th_i:], constants.ENERGIES[E_th_i:])
        for x in range(len(constants._SLANT_DEPTHS))
    ]

    if constants.get_verbose() > 1:
        print("Finished calculating true vertical underground intensities.")

    # Interpolate to the slant depths given by the set vertical depth or the user-defined slant depths

    u_intensities_tr = np.exp(
        np.interp(depths, constants._SLANT_DEPTHS, np.log(u_intensities_tr))
    )

    # Write the results to the file

    if output:

        constants._check_directory(
            os.path.join(constants.get_directory(), "underground"), force=force
        )

        if file_name == "" or not isinstance(file_name, str):

            file_name = os.path.join(
                constants.get_directory(),
                "underground",
                "{0}_Underground_Intensities_tr.txt".format(constants.get_lab()),
            )

        file_out = open(file_name, "w")

        for x in range(len(depths)):

            file_out.write(
                "{0:1.2f} {1:1.14e}\n".format(depths[x], u_intensities_tr[x])
            )

        file_out.close()

        if constants.get_verbose() > 1:

            print(f"True vertical underground intensities written to {file_name}.")

    return u_intensities_tr


# Calculate double-differential underground intensities


def _calc_u_intensities_dd(
    u_fluxes,
    E_th_i,
    primary_model,
    interaction_model,
    atmosphere,
    location,
    month,
    output,
    file_name,
    force,
):

    """This function calculates double-differential underground intensities for mountain overburdens."""

    # Check values

    constants._check_constants()

    if constants.get_overburden() != "mountain":

        raise ValueError(
            "Double-differential underground intensities can only be calculated for non-flat overburdens. See the Tutorial at {0} for an explanation.".format(
                "https://github.com/wjwoodley/mute/blob/main/docs/Tutorial.md#calculating-underground-intensities"
            )
        )

    if u_fluxes is None:

        raise Exception(
            "Underground intensities not calculated. The underground flux matrix was not calculated correctly."
        )

    # Calculate the underground intensities

    u_intensities = np.zeros(
        (len(constants._SLANT_DEPTHS), len(constants.ANGLES_FOR_S_FLUXES))
    )

    if constants.get_verbose() > 1:
        print("Calculating underground intensities.")

    for x in range(len(constants._SLANT_DEPTHS)):

        for j in range(len(constants.ANGLES_FOR_S_FLUXES)):

            u_intensities[x, j] = scii.simpson(
                u_fluxes[x, E_th_i:, j], constants.ENERGIES[E_th_i:]
            )

    if constants.get_verbose() > 1:
        print("Finished calculating underground intensities.")

    # Create an interpolator object

    try:

        _interpolator

    except NameError:

        _interpolator = _create_interpolator(u_intensities)

    # Interpolate the intensities to the slant depths and zenith angles in the mountain profile file
    # Stack the output into an array for reshaping below

    interpolation_result = np.hstack(
        [
            np.exp(_interpolator(depth, angle))
            for depth, angle in zip(
                constants._mountain_slant_depths_all, constants._mountain_zenith_all
            )
        ]
    )

    # Calculate the intensities at the zenith and azimuthal angles specified in the mountain profile file

    interp_u_intensities = np.nan_to_num(
        np.reshape(
            interpolation_result,
            (len(constants.mountain.zenith), len(constants.mountain.azimuthal)),
        )
    )

    # Write the results to the file

    if output:

        constants._check_directory(
            os.path.join(constants.get_directory(), "underground"), force=force
        )

        if file_name == "" or not isinstance(file_name, str):

            file_name = os.path.join(
                constants.get_directory(),
                "underground",
                "{0}_Underground_Intensities_dd.txt".format(constants.get_lab()),
            )

        file_out = open(file_name, "w")

        for j in range(len(constants.mountain.zenith)):

            for az in range(len(constants.mountain.azimuthal)):

                file_out.write(
                    "{0:1.5f} {1:1.5f} {2:1.14e}\n".format(
                        constants.mountain.zenith[j],
                        constants.mountain.azimuthal[az],
                        interp_u_intensities[j, az],
                    )
                )

        file_out.close()

        if constants.get_verbose() > 1:
            print(f"Underground intensities written to {file_name}.")

    return interp_u_intensities


# Calculate total underground fluxes


def calc_u_tot_flux(force=False, **kwargs):

    """
    Calculate a total underground flux in units of [(cm^2 s)^-1].

    Parameters
    ----------
    force : bool
        If True, force the calculation of new matrices and the creation of new directories if required.

    Other Parameters
    -----------------
    u_fluxes : NumPy ndarray, optional (default: calcaulated from default surface fluxes and survival probabilities)
        An underground flux tensor of shape (28, 91, 20).

    s_fluxes : NumPy ndarray, optional (default: taken from surface.load_s_fluxes_from_file())
        A surface flux matrix of shape (91, 20).

    survival_probability_tensor : NumPy ndarray, optional (default: taken from propagation.load_survival_probability_tensor_from_file())
        A survival probability tensor of shape (91, 28, 91).

    angles : array-like, optional (default: taken from consants.angles)
        An array of zenith angles in [degrees] to calculate the underground intensities for. These are only used if method is "sd" or "eq".

    E_th : float, optional (default: 0)
        An energy threshold in [MeV].

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

    Returns
    -------
    u_tot_flux : float
        The total underground flux in units of [(cm^2 s)^-1].
    """

    # Check values

    constants._check_constants()

    # Calculate the total underground flux for a flat overburdens

    if constants.get_overburden() == "flat":

        # Get the angles keyword argument
        # For use in the integral over the angles

        angles = kwargs.get("angles", constants.angles)
        angles = np.atleast_1d(angles)

        # Calculate the underground intensities

        u_intensities = calc_u_intensities(method="sd", force=force, **kwargs)

        if u_intensities is None:

            raise Exception(
                "Total underground flux not calculated. The underground intensities were not calculated properly."
            )

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

        return u_tot_flux

    # Calculate the total underground flux for mountains

    elif constants.get_overburden() == "mountain":

        # Calculate the underground intensities

        u_intensities = calc_u_intensities(method="dd", force=force, **kwargs)

        if u_intensities is None:

            raise Exception(
                "Total underground flux not calculated. The underground intensities were not calculated properly."
            )

        # Calculate the total underground flux
        # Integrate over all angles as defined by the mountain profile file

        u_tot_flux = scii.simpson(
            [
                scii.simpson(
                    u_intensities[j, :], np.radians(constants.mountain.azimuthal)
                )
                for j in range(len(constants.mountain.zenith))
            ][::-1],
            np.cos(np.radians(constants.mountain.zenith[::-1])),
        )

        return u_tot_flux

    else:

        raise NotImplementedError(
            'Overburdens of type {0} are not available. The only options are "flat" and "mountain".'.format(
                constants.get_overburden()
            )
        )
