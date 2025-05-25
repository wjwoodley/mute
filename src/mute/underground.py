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
import scipy.interpolate as sciint
import scipy.optimize as scio

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
    output=None,
    file_name="",
    force=False,
    **kwargs,
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

    output : bool, optional (default: taken from constants.get_output())
        If True, an output file will be created to store the results.

    file_name : str, optional (default: constructed from input parameters)
        The name of the file in which to store the results. If output is False or None, this is ignored.

    force : bool
        If True, force the calculation of new matrices and the creation of new directories if required.

    Other Parameters
    ----------------
    model : str in {"daemonflux", "mceq"}, optional (default: "mceq")
        The model to use to calculate surface fluxes. MCeq provides the option to specify primary, interaction, and density model keyword arguments, while daemonflux uses "gsf" as the primary model and "ddm" as the interaction model. The default model will be changed to "daemonflux" in v3.1.0. This parameter is case-insensitive.

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
    u_fluxes : NumPy ndarray
        A multi-dimensional array containing the underground fluxes. For flat overburdens, the shape will be (91, 28) for the default zenith angles, and the fluxes will be in units of [(cm^2 s sr MeV)^-1]. For mountain overburdens (or if full_tensor is True), the shape will be (28, 91, 20), and the fluxes will be in units of [(cm^2 s sr MeV^2)^-1].
    """

    # Check values

    constants._check_constants()

    if output is None:
        output = constants.get_output()

    # Get the surface flux matrix and survival probability tensor
    # If the calculation functions have to be run through the load functions, the global output setting is used

    if constants.get_verbose() > 1:
        print("Calculating underground fluxes.")

    if s_fluxes is None:
        s_fluxes = surface.load_s_fluxes_from_file(force=force, **kwargs)

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

        # Interpolate to the angles corresponding to the slant depths

        u_fluxes = sciint.interp1d(
            constants.ANGLES_FOR_S_FLUXES, u_fluxes_from_convolution
        )(constants.angles)

        # Extract the (u_energies, angles) diagonal from the convolution tensor

        u_fluxes = np.diagonal(u_fluxes, axis1=0, axis2=2)

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
                    "{0}_underground_fluxes.txt".format(constants.get_lab()),
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
                    "{0}_underground_fluxes.txt".format(constants.get_lab()),
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

    # bounds_error must be set to False because constants._mountain_slant_depths_all may contain values less than np.min(constants._SLANT_DEPTHS). This occurs when values above 14.0 km.w.e. are set to 0 km.w.e. (which is outside the bounds of 0.5 km.w.e. to 14.0 km.w.e.) when the mountain file is loaded.

    _interpolator = sciint.RegularGridInterpolator(
        (
            constants._SLANT_DEPTHS,
            np.cos(np.radians(constants.ANGLES_FOR_S_FLUXES))[::-1],
        ),
        np.log(u_intensities)[:, ::-1],
        method="linear",
        bounds_error=False,
        fill_value=-np.inf,
    )

    return _interpolator


# Calculate the underground intensities


def calc_u_intensities(
    method=None,
    output=None,
    file_name="",
    force=False,
    **kwargs,
):
    """
    Calculate underground intensities in units of [(cm^2 s sr)^-1].

    Parameters
    ----------
    method : str, optional (default: "sd" if overburden is flat and "dd" if overburden is mountain)
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

    model : str in {"daemonflux", "mceq"}, optional (default: "mceq")
        The model to use to calculate surface fluxes. MCeq provides the option to specify primary, interaction, and density model keyword arguments, while daemonflux uses "gsf" as the primary model and "ddm" as the interaction model. The default model will be changed to "daemonflux" in v3.1.0. This parameter is case-insensitive.

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
    u_intensities : NumPy ndarray
        An array containing the underground intensities. For flat overburdens (if method is "sd", "eq", or "tr"), the length will be that of angles, and the intensities will be in units of [(cm^2 s sr)^-1]. For mountain overburdens (if method is "dd"), the shape will be (len(constants.mountain.zenith), len(constants.mountain.azimuthal)), and the intensities will be in units of [(cm^2 s sr km.w.e.)^-1].
    """

    # Check values

    constants._check_constants()

    if type(method) == str:
        method = method.lower()

    assert method in [
        None,
        "sd",
        "eq",
        "tr",
        "dd",
    ], '"{0}" is not a valid method. Use "sd", "tr", "eq", or "dd". See the Tutorial at {1} for an explanation.'.format(
        method,
        "https://github.com/wjwoodley/mute/blob/main/docs/Tutorial.md#calculating-underground-intensities",
    )

    if method is None:

        if constants.get_overburden() == "flat":
            method = "sd"
        else:
            method = "dd"

    if output is None:
        output = constants.get_output()

    # Get the keyword arguments

    u_fluxes = kwargs.get("u_fluxes", None)
    s_fluxes = kwargs.get("s_fluxes", None)
    survival_probability_tensor = kwargs.get("survival_probability_tensor", None)
    angles = kwargs.get("angles", constants.angles)
    depths = kwargs.get("depths", constants.slant_depths)
    E_th = kwargs.get("E_th", 0)
    model = kwargs.get("model", "mceq")

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
            s_fluxes = surface.load_s_fluxes_from_file(**kwargs)
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

    required_shape = (
        len(constants._SLANT_DEPTHS),
        len(constants.ENERGIES),
        len(constants.ANGLES_FOR_S_FLUXES),
    )

    assert (
        u_fluxes.shape == required_shape
    ), f"The underground flux tensor does not have the correct shape. The shape must be {required_shape}, not {u_fluxes.shape}. Set full_tensor to True in the underground.calc_u_fluxes() function in order to return a tensor of shape {required_shape}. Alternatively, try passing a surface flux matrix and / or a survival probability tensor in as parameters."

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
            output,
            file_name,
            force,
        )

    # Calculate the double-differential underground intensities

    elif method == "dd":

        return _calc_u_intensities_dd(
            u_fluxes,
            E_th_i,
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

    if constants.get_verbose() > 1:
        print("Calculating underground intensities.")

    u_intensities = scii.simpson(
        u_fluxes[:, E_th_i:, :], x=constants.ENERGIES[E_th_i:], axis=1
    )

    # Create an interpolator for the underground intensities
    # Inverse the angles axis because cos(theta) decreases as theta increases

    u_intensities_interpolator = sciint.RectBivariateSpline(
        constants._SLANT_DEPTHS,
        np.cos(np.radians(constants.ANGLES_FOR_S_FLUXES))[::-1],
        np.log(u_intensities)[:, ::-1],
        kx=1,
        ky=1,
    )

    # Interpolate the underground intensities to the set slant depths and zenith angles
    # Set grid = False so cos(theta) does not have to be strictly increasing
    # Take the exponent because the logarithm of the intensities was used to create the interpolator

    interp_u_intensities = np.exp(
        u_intensities_interpolator(
            constants.slant_depths, np.cos(np.radians(constants.angles)), grid=False
        )
    )

    # Interpolate intensities to user-input angles
    # If the user has not specified angles, return the whole matrix with no interpolation
    # Check to make sure both arrays have the same number of values, and that those values are equal at each index

    if not (
        len(angles) == len(constants.angles) and np.allclose(angles, constants.angles)
    ):

        interp_u_intensities = np.exp(
            np.interp(
                1 / np.cos(np.radians(angles)),
                1 / np.cos(np.radians(constants.angles)),
                np.log(interp_u_intensities),
            )
        )

    # Multiply by cos(theta) if calculating vertical-equivalent intensities

    if method == "eq":

        interp_u_intensities = np.cos(np.radians(angles)) * interp_u_intensities

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
                "{0}_underground_intensities_{1}.txt".format(
                    constants.get_lab(), method
                ),
            )

        file_out = open(file_name, "w")

        for j in range(len(angles)):

            file_out.write(
                "{0:1.5f} {1:1.14e}\n".format(angles[j], interp_u_intensities[j])
            )

        file_out.close()

        if constants.get_verbose() > 1:
            print(f"Underground intensities written to {file_name}.")

    return interp_u_intensities


# Calculate true vertical underground intensities


def _calc_u_intensities_tr(
    u_fluxes,
    depths,
    E_th_i,
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

    u_intensities = scii.simpson(
        u_fluxes[:, E_th_i:, :], x=constants.ENERGIES[E_th_i:], axis=1
    )

    u_intensities_tr = u_intensities[:, 0]

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
                "{0}_underground_intensities_tr.txt".format(constants.get_lab()),
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
                u_fluxes[x, E_th_i:, j], x=constants.ENERGIES[E_th_i:]
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
            np.exp(_interpolator(np.column_stack([depth, np.cos(np.radians(angle))])))
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
                "{0}_underground_intensities_dd.txt".format(constants.get_lab()),
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


# Calculate underground angular distributions


def calc_u_ang_dist(kind, output=None, file_name="", force=False, **kwargs):
    """
    Calculate an underground angular distribution in units of [(cm^2 s)^-1] for "zenith" or [(cm^2 s rad)^-1] for "azimuthal". This function can only be used for mountain overburdens.

    Parameters
    ----------
    kind : str in {"zenith", "azimuthal"}
        The kind of angular distribution to calculate. Options:
        zenith    = Zenith angular distribution (includes integration over the azimuthal angle)
        azimuthal = Azimuthal angular distribution (includes integration over the zenith angle)

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

    model : str in {"daemonflux", "mceq"}, optional (default: "mceq")
        The model to use to calculate surface fluxes. MCeq provides the option to specify primary, interaction, and density model keyword arguments, while daemonflux uses "gsf" as the primary model and "ddm" as the interaction model. The default model will be changed to "daemonflux" in v3.1.0. This parameter is case-insensitive.

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
    u_ang_dist : NumPy ndarray
        An array containing the underground energy spectrum. The length will be that of constants.mountain.zenith for "zenith" or constants.mountain.azimuthal for "azimuthal", and the spectrum will be in units of [(cm^2 s)^-1] for "zenith" or [(cm^2 s rad)^-1] for "azimuthal".
    """

    # Check values

    constants._check_constants()

    assert (
        constants.get_overburden() == "mountain"
    ), 'Angular distributions can only be calculated for labs under mountains. Please set the overburden to mountain with mtc.set_overburden("mountain").'

    assert kind in ["zenith", "azimuthal"], 'kind must be "zenith" or "azimuthal".'

    if output is None:
        output = constants.get_output()

    # Get the underground intensity matrix

    u_intensities = calc_u_intensities(
        method="dd", output=output, force=force, **kwargs
    )

    # Calculate the zenith angular distribution

    if kind == "zenith":

        # Integrate over phi for zenith angular distribution

        u_ang_dist = scii.simpson(
            u_intensities, x=np.radians(constants.mountain.azimuthal), axis=1
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
                    "{0}_underground_zenith_angular_distribution.txt".format(
                        constants.get_lab()
                    ),
                )

            file_out = open(file_name, "w")

            for j in range(len(constants.mountain.zenith)):

                file_out.write(
                    "{0:1.5f} {1:1.14e}\n".format(
                        constants.mountain.zenith[j], u_ang_dist[j]
                    )
                )

            file_out.close()

            if constants.get_verbose() > 1:
                print(
                    f"Underground zenith angular distribution written to {file_name}."
                )

        return u_ang_dist

    # Calculate the azimuthal angular distribution

    elif kind == "azimuthal":

        # Integrate over cos(theta) for azimuthal angular distribution

        u_ang_dist = abs(
            scii.simpson(
                u_intensities, x=np.cos(np.radians(constants.mountain.zenith)), axis=0
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
                    "{0}_underground_azimuthal_angular_distribution.txt".format(
                        constants.get_lab()
                    ),
                )

            file_out = open(file_name, "w")

            for az in range(len(constants.mountain.azimuthal)):

                file_out.write(
                    "{0:1.5f} {1:1.14e}\n".format(
                        constants.mountain.azimuthal[az], u_ang_dist[az]
                    )
                )

            file_out.close()

            if constants.get_verbose() > 1:
                print(
                    f"Underground azimuthal angular distribution written to {file_name}."
                )

        return u_ang_dist

    else:

        raise ValueError(f'Unknown kind {kind}. Please use "zenith" or "azimuthal".')


# Calculate underground energy distributions


def calc_u_e_spect(output=None, file_name="", force=False, **kwargs):
    """
    Calculate an underground energy spectrum in units of [(cm^2 s sr MeV)^-1].

    Parameters
    ----------
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

    E_th : float, optional (default: 0)
        An energy threshold in [MeV].

    model : str in {"daemonflux", "mceq"}, optional (default: "mceq")
        The model to use to calculate surface fluxes. MCeq provides the option to specify primary, interaction, and density model keyword arguments, while daemonflux uses "gsf" as the primary model and "ddm" as the interaction model. The default model will be changed to "daemonflux" in v3.1.0. This parameter is case-insensitive.

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
    u_e_spect : NumPy ndarray
        An array containing the underground energy spectrum. The length will be that of constants.ENERGIES, and the spectrum will be in units of [(cm^2 s sr MeV)^-1].
    """

    # Check values

    constants._check_constants()

    if output is None:
        output = constants.get_output()

    # Get the keyword arguments

    u_fluxes = kwargs.get("u_fluxes", None)
    s_fluxes = kwargs.get("s_fluxes", None)
    survival_probability_tensor = kwargs.get("survival_probability_tensor", None)
    angles = kwargs.get("angles", constants.angles)
    depths = kwargs.get("depths", constants.slant_depths)
    E_th = kwargs.get("E_th", 0)

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

            s_fluxes = surface.load_s_fluxes_from_file(**kwargs)

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

    required_shape = (
        len(constants._SLANT_DEPTHS),
        len(constants.ENERGIES),
        len(constants.ANGLES_FOR_S_FLUXES),
    )

    assert (
        u_fluxes.shape == required_shape
    ), f"The underground flux tensor does not have the correct shape. The shape must be {required_shape}, not {u_fluxes.shape}. Set full_tensor to True in the underground.calc_u_fluxes() function in order to return a tensor of shape {required_shape}. Alternatively, try passing a surface flux matrix and / or a survival probability tensor in as parameters."

    # Calculate the underground energy spectrum for flat overburdens

    if constants.get_overburden() == "flat":

        interp_result = np.zeros(
            (len(constants.ENERGIES[E_th_i:]), len(constants.ANGLES_FOR_S_FLUXES))
        )

        with np.errstate(all="ignore"):

            # Create an interpolator function

            interpolator = sciint.RegularGridInterpolator(
                (
                    constants._SLANT_DEPTHS,
                    constants.ENERGIES,
                    constants.ANGLES_FOR_S_FLUXES,
                ),
                np.log(u_fluxes),
                method="linear",
                bounds_error=False,
                fill_value=-np.inf,
            )

            # Loop over the zenith angles to specify the slant depth

            for j in range(len(constants.ANGLES_FOR_S_FLUXES)):

                interp_result[:, j] = np.nan_to_num(
                    np.exp(
                        interpolator(
                            (
                                constants.get_vertical_depth()
                                / np.cos(np.radians(constants.ANGLES_FOR_S_FLUXES[j])),
                                constants.ENERGIES[E_th_i:],
                                constants.ANGLES_FOR_S_FLUXES[j],
                            )
                        )
                    )
                )

        # Integrate over the zenith and azimuthal angles to calculate the energy spectrum

        u_e_spect = (
            2
            * np.pi
            * abs(
                scii.simpson(
                    interp_result,
                    x=np.cos(np.radians(constants.ANGLES_FOR_S_FLUXES)),
                    axis=1,
                )
            )
        )

    # Calculate the underground energy spectrum for mountains

    elif constants.get_overburden() == "mountain":

        interp_result = np.zeros(
            (
                len(constants.mountain.zenith),
                len(constants.mountain.azimuthal),
                len(constants.ENERGIES[E_th_i:]),
            )
        )

        with np.errstate(all="ignore"):

            # Create an interpolator function

            interpolator = sciint.RegularGridInterpolator(
                (
                    constants._SLANT_DEPTHS,
                    constants.ENERGIES,
                    constants.ANGLES_FOR_S_FLUXES,
                ),
                np.log(u_fluxes),
                method="linear",
                bounds_error=False,
                fill_value=-np.inf,
            )

            # Loop over the zenith and azimuthal angles to specify the slant depth index
            # Interpolate the underground fluxes to the slant depths and zenith angles
            # Implicitly reshape the resulting tensor so it is in terms of zenith and azimuthal angle

            for j in range(len(constants.mountain.zenith)):

                for az in range(len(constants.mountain.azimuthal)):

                    interp_result[j, az, :] = np.nan_to_num(
                        np.exp(
                            interpolator(
                                (
                                    constants.mountain.slant_depths[j, az],
                                    constants.ENERGIES[E_th_i:],
                                    constants.mountain.zenith[j],
                                )
                            )
                        )
                    )

        # Integrate over the zenith and azimuthal angles to calculate the energy spectrum

        u_e_spect = scii.simpson(
            abs(
                scii.simpson(
                    interp_result,
                    x=np.cos(np.radians(constants.mountain.zenith)),
                    axis=0,
                )
            ),
            x=np.radians(constants.mountain.azimuthal),
            axis=0,
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
                "{0}_underground_energy_spectrum.txt".format(constants.get_lab()),
            )

        file_out = open(file_name, "w")

        for i in range(len(constants.ENERGIES[E_th_i:])):

            file_out.write(
                "{0:1.14f} {1:1.14e}\n".format(
                    constants.ENERGIES[E_th_i + i], u_e_spect[i]
                )
            )

        file_out.close()

        if constants.get_verbose() > 1:

            print(f"Underground energy spectrum written to {file_name}.")

    # Return the underground energy spectrum

    return u_e_spect


# Calculate mean and median underground energies


def calc_u_mean_e(force=False, **kwargs):
    """
    Calculate mean and median underground energies and their 68% and 95% confidence intervals in units of [MeV].

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

    E_th : float, optional (default: 0)
        An energy threshold in [MeV].

    model : str in {"daemonflux", "mceq"}, optional (default: "mceq")
        The model to use to calculate surface fluxes. MCeq provides the option to specify primary, interaction, and density model keyword arguments, while daemonflux uses "gsf" as the primary model and "ddm" as the interaction model. The default model will be changed to "daemonflux" in v3.1.0. This parameter is case-insensitive.

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
    An array containing the following values:
        u_mean_e          = The mean underground energy
        u_median_e        = The median underground energy
        u_median_e_68_pos = The positive 68% confidence interval on the median underground energy
        u_median_e_68_neg = The negative 68% confidence interval on the median underground energy
        u_median_e_95_pos = The positive 95% confidence interval on the median underground energy
        u_median_e_95_neg = The negative 95% confidence interval on the median underground energy

    All values will be in units of [MeV]. Both positive and negative values are returned for the 68% and 95% confidence intervals because the intervals are asymmetric.
    """

    # Check values

    constants._check_constants()

    # Calculate the energy spectrum

    u_e_spect = calc_u_e_spect(force=force, **kwargs)
    E_th_i = len(constants.ENERGIES) - len(u_e_spect)

    # Calculate the mean as the first moment of the energy spectrum

    if constants.get_verbose() > 1:
        print("Calculating mean underground energies.")

    u_mean_e = scii.simpson(
        constants.ENERGIES[E_th_i:] * u_e_spect, x=constants.ENERGIES[E_th_i:]
    ) / scii.simpson(u_e_spect, x=constants.ENERGIES[E_th_i:])

    # Calculate the cumulative integral to find the median
    # Normalise by dividing by the last element (the maximum)

    cumulative_dist = scii.cumulative_trapezoid(u_e_spect, constants.ENERGIES[E_th_i:])
    cumulative_dist_norm = cumulative_dist / cumulative_dist[-1]

    # Find the roots of the cumulative distribution

    def cumulative_dist_roots(value_in):

        return np.exp(
            sciint.UnivariateSpline(
                np.log(constants.ENERGIES[(E_th_i + 1) :]),
                cumulative_dist_norm - value_in,
                k=3,
                s=0,
                ext="zeros",
            ).roots()
        )

    # Calculate the median and positive and negative 68% and 95% confidence intervals

    u_median_e = cumulative_dist_roots(0.5)[0]
    u_median_e_68_pos = cumulative_dist_roots(0.68)[0]
    u_median_e_68_neg = cumulative_dist_roots(1 - 0.68)[0]
    u_median_e_95_pos = cumulative_dist_roots(0.95)[0]
    u_median_e_95_neg = cumulative_dist_roots(1 - 0.95)[0]

    if constants.get_verbose() > 1:
        print("Finished calculating mean underground energies.")

    # Return the results in an array

    return np.array(
        (
            u_mean_e,
            u_median_e,
            u_median_e_68_pos,
            u_median_e_68_neg,
            u_median_e_95_pos,
            u_median_e_95_neg,
        )
    )


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

    model : str in {"daemonflux", "mceq"}, optional (default: "mceq")
        The model to use to calculate surface fluxes. MCeq provides the option to specify primary, interaction, and density model keyword arguments, while daemonflux uses "gsf" as the primary model and "ddm" as the interaction model. The default model will be changed to "daemonflux" in v3.1.0. This parameter is case-insensitive.

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
    u_tot_flux : float
        The total underground flux in units of [(cm^2 s)^-1].
    """

    # Check values

    constants._check_constants()

    # Calculate the total underground flux for flat overburdens

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
        # Therefore, take the absolute value
        # Otherwise, the answer will be negative

        u_tot_flux = (
            2 * np.pi * abs(scii.simpson(u_intensities, x=np.cos(np.radians(angles))))
        )

        return float(u_tot_flux)

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
                    u_intensities[j, :], x=np.radians(constants.mountain.azimuthal)
                )
                for j in range(len(constants.mountain.zenith))
            ][::-1],
            x=np.cos(np.radians(constants.mountain.zenith[::-1])),
        )

        return float(u_tot_flux)

    else:

        raise NotImplementedError(
            'Overburdens of type {0} are not available. The only options are "flat" and "mountain".'.format(
                constants.get_overburden()
            )
        )


def calc_depth(kind, u_tot_flux=None, force=False, **kwargs):
    """
    Calculate different kinds of depths from the loaded mountain profile in units of [km.w.e.].

    Parameters
    ----------
    kind : str in {"ev", "v", "min", "max", "avg", "mean"}
        The kind of depth to calculate. Options:
        ev   = The equivalent vertical depth obtained by fitting the input or calculated total flux to a total flux curve for standard rock
        v    = The straight vertical depth (engineering depth) at the center of the mountain profile
        min  = The minimum depth in the mountain profile file
        max  = The maximum depth in the mountain profile file
        avg  = The average depth in the mountain profile file
        mean = An alias for "avg"

    u_tot_flux : float, optional (default: None)
        A total flux value to fit to the internally calculated total flux curve. If None, this is calculated using the default parameters in calc_u_tot_flux().

    force : bool
        If True, force the calculation of new matrices and the creation of new directories if required.

    Other Parameters
    -----------------
    u_fluxes : NumPy ndarray, optional (default: calcaulated from default surface fluxes and survival probabilities)
        An underground flux tensor of shape (28, 91, 20).

    s_fluxes : NumPy ndarray, optional (default: taken from surface.load_s_fluxes_from_file())
        A surface flux matrix of shape (91, 20).

    survival_probability_tensor : NumPy ndarray, optional (default: taken from rock_2.65_1000000_survival_probabilities.npy)
        A survival probability tensor of shape (91, 28, 91). This parameter allows the total flux to be fit to a curve for a user-specified medium. If no survival probability tensor is specified, the pre-calculated tensor for standard rock will be used.

    E_th : float, optional (default: 0)
        An energy threshold in [MeV].

    model : str in {"daemonflux", "mceq"}, optional (default: "mceq")
        The model to use to calculate surface fluxes. MCeq provides the option to specify primary, interaction, and density model keyword arguments, while daemonflux uses "gsf" as the primary model and "ddm" as the interaction model. The default model will be changed to "daemonflux" in v3.1.0. This parameter is case-insensitive.

    return_error : bool, optional (default: False)
        If True, if model == "daemonflux", return the error on the surface fluxes instead of the surface fluxes. This moves the curve that u_tot_flux is fit to slightly and can be used to calculate uncertainties on the depth.

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
    depth : float
        The depth of kind kind in units of [km.w.e.].
    """

    # Check values

    constants._check_constants()

    assert (
        constants.get_overburden() == "mountain"
    ), 'Depths can only be calculated for labs under mountains. Please set the overburden to mountain with mtc.set_overburden("mountain").'

    assert kind in ["ev", "v", "min", "max", "avg", "mean"], f"Unknown kind {kind}."

    # Equivalent vertical depth (fit to standard rock total flux curve)

    if kind == "ev":

        # Calculate the total flux value to fit

        if u_tot_flux is None:
            u_tot_flux = calc_u_tot_flux(force=force, **kwargs)

        # Get the keyword arguments

        u_fluxes = kwargs.get("u_fluxes", None)
        s_fluxes = kwargs.get("s_fluxes", None)
        survival_probability_tensor = kwargs.get("survival_probability_tensor", None)

        # Calculate the underground fluxes
        # Get the surface flux matrix and survival probability tensor
        # If the calculation functions have to be run through the getter functions, the global output setting is used

        if u_fluxes is None:

            if constants.get_verbose() > 1:
                print("Calculating underground fluxes.")

            if s_fluxes is None:
                s_fluxes = surface.load_s_fluxes_from_file(force=force, **kwargs)

            if survival_probability_tensor is None:
                survival_probability_tensor = (
                    propagation.load_survival_probability_tensor_from_file(
                        os.path.join(
                            constants.get_directory(),
                            "survival_probabilities",
                            "rock_2.65_1000000_survival_probabilities.npy",
                        )
                    )
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

        required_shape = (
            len(constants._SLANT_DEPTHS),
            len(constants.ENERGIES),
            len(constants.ANGLES_FOR_S_FLUXES),
        )

        assert (
            u_fluxes.shape == required_shape
        ), f"The underground flux tensor does not have the correct shape. The shape must be {required_shape}, not {u_fluxes.shape}. Set full_tensor to True in the underground.calc_u_fluxes() function in order to return a tensor of shape {required_shape}. Alternatively, try passing a surface flux matrix and / or a survival probability tensor in as parameters."

        # Calculate a total flux curve to fit to
        # Loop over depths to construct the curve
        # The fit must be done to a flat overburden curve (this is the meaning of "ev" or "'equivalent vertical' depth")
        # Because the overburden is set to "mountain", the total flux calculation must be recreated for flat overburdens for each depth

        # Calculate underground intensities to integrate over
        # See _calc_u_intensities_sd() for details

        u_intensities = scii.simpson(u_fluxes, x=constants.ENERGIES, axis=1)
        u_intensities_interpolator = sciint.RectBivariateSpline(
            constants._SLANT_DEPTHS,
            np.cos(np.radians(constants.ANGLES_FOR_S_FLUXES))[::-1],
            np.log(u_intensities)[:, ::-1],
            kx=1,
            ky=1,
        )

        # Define and loop over depths

        depths = np.arange(0.5, 8.1, 0.5)

        u_tot_fluxes = np.zeros(len(depths))

        for xi in range(len(depths)):

            slant_depths = constants._SLANT_DEPTHS[
                constants._SLANT_DEPTHS >= depths[xi]
            ]
            angles = np.degrees(np.arccos(depths[xi] / slant_depths))

            interp_u_intensities = np.exp(
                u_intensities_interpolator(
                    slant_depths, np.cos(np.radians(angles)), grid=False
                )
            )

            # Calculate total underground flux

            u_tot_fluxes[xi] = (
                2
                * np.pi
                * abs(scii.simpson(interp_u_intensities, x=np.cos(np.radians(angles))))
            )

        # Minimise the depth

        def minimise(X, fluxes, root):

            spline = sciint.UnivariateSpline(
                depths[1:-1], np.log(fluxes[1:-1]), k=1, s=0
            )

            return (spline(X) - np.log(root)) ** 2

        # Calculate fitted depth values

        h = scio.minimize(minimise, 1, args=(u_tot_fluxes, u_tot_flux)).x[0]

    # Vertical depth (engineering depth)

    elif kind == "v":

        h = constants.mountain.slant_depths[0, 0]

    # Minimum depth

    elif kind == "min":

        h = np.min(constants.mountain.slant_depths)

    # Maximum depth

    elif kind == "max":

        h = np.max(constants.mountain.slant_depths)

    # Average depth

    elif kind == "avg" or kind == "mean":

        h = np.mean(constants.mountain.slant_depths)

    return float(h)
