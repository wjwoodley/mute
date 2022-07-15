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
from tqdm import tqdm

import mute.constants as constants

try:

    import proposal as pp

except ImportError:

    pass

# Create the propagator


def _create_propagator(force):

    """This function creates the propagator object in PROPOSAL for use in propagation._propagation_loop()"""

    # Check values

    constants._check_constants(force=force)

    # The creation of the propagator can be very slow the first time
    # The propagator is used in every iteration of the doubly-nested propagation loop
    # Make it a global variable so it only has to be created once

    global propagator

    if constants.get_verbose() > 1:
        print("Creating propagator.")

    # Propagator arguments

    mu = pp.particle.MuMinusDef()
    cuts = pp.EnergyCutSettings(500, 0.05, True)

    if constants.get_medium() == "rock":

        medium = pp.medium.StandardRock()

    elif constants.get_medium() == "water":

        medium = pp.medium.Water()

    elif constants.get_medium() == "ice":

        medium = pp.medium.Ice()

    elif constants.get_medium() == "air":

        medium = pp.medium.Air()

    else:

        raise NotImplementedError(
            "Medium type {0} not implemented. The only options are {1}, {2}, {3}, and {4}.".format(
                constants.get_medium(), "rock", "water", "ice", "air"
            )
        )

    args = {"particle_def": mu, "target": medium, "interpolate": True, "cuts": cuts}

    # Initialise standard cross-sections, then specify and set parametrisation models

    cross_sections = pp.crosssection.make_std_crosssection(**args)

    brems_param = pp.parametrization.bremsstrahlung.KelnerKokoulinPetrukhin(lpm=False)
    epair_param = pp.parametrization.pairproduction.KelnerKokoulinPetrukhin(lpm=False)
    ionis_param = pp.parametrization.ionization.BetheBlochRossi(energy_cuts=cuts)
    shado_param = pp.parametrization.photonuclear.ShadowButkevichMikheyev()
    photo_param = pp.parametrization.photonuclear.AbramowiczLevinLevyMaor97(
        shadow_effect=shado_param
    )

    cross_sections[0] = pp.crosssection.make_crosssection(brems_param, **args)
    cross_sections[1] = pp.crosssection.make_crosssection(epair_param, **args)
    cross_sections[2] = pp.crosssection.make_crosssection(ionis_param, **args)
    cross_sections[3] = pp.crosssection.make_crosssection(photo_param, **args)

    # Propagation utility

    collection = pp.PropagationUtilityCollection()

    collection.interaction = pp.make_interaction(cross_sections, True)
    collection.displacement = pp.make_displacement(cross_sections, True)
    collection.time = pp.make_time(cross_sections, mu, True)
    collection.decay = pp.make_decay(cross_sections, mu, True)

    pp.PropagationUtilityCollection.cont_rand = False

    utility = pp.PropagationUtility(collection=collection)

    # Other settings

    pp.do_exact_time = False

    # Set up geometry

    detector = pp.geometry.Sphere(
        position=pp.Cartesian3D(0, 0, 0), radius=10000000, inner_radius=0
    )
    density_distr = pp.density_distribution.density_homogeneous(
        mass_density=constants.get_density()
    )

    propagator = pp.Propagator(mu, [(detector, utility, density_distr)])

    if constants.get_verbose() > 1:
        print("Finished creating propagator.")

    return propagator


# Propagation function


def _propagation_loop(energy, slant_depth, force):

    """This function propagates n_muon muons, looping over the energies and slant depths, and returns the muons' underground energies in [MeV]."""

    # Check values

    constants._check_constants(force=force)

    n_muon = constants.get_n_muon()

    # Convert the slant depth from [km.w.e.] to [cm]

    convert_to_cm = 1e5 * 0.997 / constants.get_density()

    # Initialise a list of underground energies

    u_energies_ix = []

    # Define the initial state of the muon

    mu_initial = pp.particle.ParticleState()
    mu_initial.energy = energy + constants._MU_MASS
    mu_initial.position = pp.Cartesian3D(0, 0, 0)
    mu_initial.direction = pp.Cartesian3D(0, 0, -1)

    # Propagate n_muon muons

    for _ in range(n_muon):

        # Propagate the muons

        track = propagator.propagate(mu_initial, slant_depth * convert_to_cm)

        # Test whether or not the muon has energy left (has not lost all of its energy or has not decayed)
        # If it does, record its energy
        # If it does not, ignore this muon and proceed with the next loop iteration

        if (
            track.track_energies()[-1] != constants._MU_MASS
            and track.track_types()[-1] != pp.particle.Interaction_Type.decay
        ):

            # Store the final underground energy of the muon

            u_energies_ix.append(track.track_energies()[-1])

    # Return the underground energies for the muon

    return u_energies_ix


# Propagate the muons and return underground energies


def propagate_muons(seed=0, job_array_number=0, output=None, force=False):

    """
    Propagate muons for the default surface energy grid and slant depths.

    The default surface energy grid is given by constants.ENERGIES, and the default slant depths are given by constants._SLANT_DEPTHS.

    Parameters
    ----------
    seed : int, optional (default: 0)
        The random seed for use in the PROPOSAL propagator.

    job_array_number : int, optional (default: 0)
        The job array number from a high-statistics run on a computer cluster. This is set so the underground energy files from each job in the job array will be named differently.

    output : bool, optional (default: taken from constants.get_output())
        If True, an output file will be created to store the results.

    force : bool, optional (default: False)
        If True, force the creation of new directories if required.

    Returns
    -------
    u_energies : NumPy ndarray
        A two-dimensional array containing lists of underground energies for muons that survived the propagation. The shape of the array will be (91, 28), and the underground energies will be in units of [MeV].
    """

    # Check values

    constants._check_constants(force=force)

    if output is None:
        output = constants.get_output()

    assert type(job_array_number) == int, "job_array_number must be an integer."

    # Create the propagator once

    _create_propagator(force=force)

    # Set the random seed

    pp.RandomGenerator.get().set_seed(seed)

    # Initialise a matrix of underground energies

    u_energies = np.zeros(
        (len(constants.ENERGIES), len(constants._SLANT_DEPTHS)), dtype=np.ndarray
    )

    # Run the propagation function

    if constants.get_verbose() >= 1:

        print(
            "Propagating {0} muons.".format(
                constants.get_n_muon()
                * len(constants.ENERGIES)
                * len(constants._SLANT_DEPTHS)
            )
        )

    for i in (
        tqdm(range(len(constants.ENERGIES)))
        if constants.get_verbose() >= 1
        else range(len(constants.ENERGIES))
    ):

        for x in range(len(constants._SLANT_DEPTHS)):

            u_energies[i, x] = _propagation_loop(
                constants.ENERGIES[i], constants._SLANT_DEPTHS[x], force=force
            )

    if constants.get_verbose() >= 1:
        print("Finished propagation.")

    # Write the results to the file

    if output:

        constants._check_directory(
            os.path.join(constants.get_directory(), "underground_energies"), force=force
        )

        file_name = os.path.join(
            constants.get_directory(),
            "underground_energies",
            "{0}_{1}_{2}_Underground_Energies_{3}.npy".format(
                constants.get_medium(),
                constants.get_density(),
                constants.get_n_muon(),
                job_array_number,
            ),
        )

        np.save(file_name, u_energies)

        if constants.get_verbose() > 1:
            print(f"Underground energies written to {file_name}.")

    return u_energies


# Load underground energies


def _load_u_energies_from_files(file_name_pattern, n_job=1, force=False):

    """
    Load the underground energies resulting from the PROPOSAL Monte Carlo from a file or collection of files stored in data/underground_energies.

    Parameters
    ----------
    file_name_pattern : str
        The file name pattern for the file(s) that the underground energy data is stored in. This should end before an underscore so the function can append the job array number. For example, pass "underground_energies" for a set of files beginning with "underground_energies_0.npy".

    n_job : int, optional (default: 1)
        The number of jobs that were run on the computer cluster. Set this to the number of files the underground energies are spread across.

    force : bool, optional (default: False)
        If True, force the propagation of muons and the creation of new directories if required.

    Returns
    -------
    u_energies : NumPy ndarray
        A two-dimensional array containing lists of underground energies for muons that survived the propagation. The shape of the array will be (91, 28), and the underground energies will be in units of [MeV].
    """

    # Check values

    constants._check_constants(force=force)

    # Check that the directory exists

    if not os.path.exists(
        os.path.join(constants.get_directory(), "underground_energies")
    ):

        if constants.get_verbose() >= 1:

            print(
                f"{constants.get_directory}/underground_energies does not exist. Underground energies not loaded."
            )

        return

    # Test if the file exists

    if not os.path.isfile(f"{file_name_pattern}_0.npy"):

        if constants.get_verbose() >= 1:

            print(
                f"{file_name_pattern}_0.npy does not exist. Underground energies not loaded."
            )

        return

    # Fill a u_energies array with empty lists that will be able to be extended

    u_energies = np.empty(
        (len(constants.ENERGIES), len(constants._SLANT_DEPTHS)), dtype=object
    )

    for i in np.ndindex(u_energies.shape):
        u_energies[i] = []

    # Loop over all output files and add the contents to u_energies

    for job_array_number in (
        tqdm(range(n_job)) if constants.get_verbose() >= 1 else range(n_job)
    ):

        u_energies += np.load(
            f"{file_name_pattern}_{job_array_number}.npy", allow_pickle=True
        )

    if constants.get_verbose() > 1:
        print("Loaded underground energies.")

    return u_energies


# Calculate survival probabilities


def calc_survival_probability_tensor(
    seed=0, file_name_pattern=None, n_job=1, output=None, file_name="", force=False
):

    """
    Calculate survival probabilities in units of [(MeV^2 km.w.e.)^-1] for the default surface energy grid and slant depths.

    The default surface energy grid is given by constants.ENERGIES, and the default slant depths are given by constants._SLANT_DEPTHS. If the propagation of muons has already been done, this will load the underground energies file, unless force is set to True.

    Parameters
    ----------
    seed : int, optional (default: 0)
        The random seed for use in the PROPOSAL propagator.

    file_name_pattern : str, optional
        The file name pattern for the file(s) that the underground energy data is stored in. This should end before an underscore so the function can append the job array number. For example, pass "underground_energies" for a set of files beginning with "underground_energies_0.npy".

    n_job : int, optional (default: 1)
        The number of jobs that were run on the computer cluster. Set this to the number of files the underground energy data is spread across.

    output : bool, optional (default: taken from constants.get_output())
        If True, an output file will be created to store the results.

    file_name : str, optional (default: constructed from set global propagation constants)
        The name of the file in which to store the results. If output is False or None, this is ignored.

    force : bool, optional (default: False)
        If True, this will force the muons to be propagated in the Monte Carlo whether an underground energies file already exists or not.

    Returns
    -------
    survival_probability_tensor : NumPy ndarray
        A three-dimensional array containing the survival probabilities. The shape of the array will be (91, 28, 91), and the survival probabilities will be in units of [(MeV^2 km.w.e.)^-1].
    """

    # Check values

    if output is None:
        output = constants.get_output()

    # Construct a default file name pattern

    file_name_pattern_default = os.path.join(
        constants.get_directory(),
        "underground_energies",
        "{0}_{1}_{2}_Underground_Energies".format(
            constants.get_medium(),
            constants.get_density(),
            int(constants.get_n_muon() / n_job),
        ),
    )

    # Check if propagate_muons() should be forced or not
    # If not, check if underground energy files exist
    # If not, ask if muons should be propagated

    if force:

        u_energies = propagate_muons(seed=seed, output=output, force=force)

    else:

        # Check if the user has specified underground energies to load
        # If not, look for the default file name pattern and check if it exists
        # If so, load the underground energies
        # If not, ask if muons should be propagated

        if file_name_pattern is not None:

            u_energies = _load_u_energies_from_files(
                file_name_pattern=os.path.join(
                    constants.get_directory(), "underground_energies", file_name_pattern
                ),
                n_job=n_job,
                force=force,
            )

        elif os.path.isfile(f"{file_name_pattern_default}_0.npy"):

            u_energies = _load_u_energies_from_files(
                file_name_pattern=file_name_pattern_default, n_job=n_job, force=force
            )

        else:

            answer = input(
                "No underground energy file currently exists for the set medium, density, or number of muons. Would you like to create one (y/n)?: "
            )

            if answer.lower() == "y":

                u_energies = propagate_muons(seed=seed, output=output, force=force)

            else:

                print("Underground energies not calculated.")
                print("Survival probabilities not calculated.")

                return

    # Check that the underground energies were loaded

    if u_energies is None:

        raise Exception(
            "Survival probabilities not calculated. The underground energies were not loaded or calculated correctly."
        )

    # Calculate the survival probabilities
    # Zeroth index = Surface energy
    # First index  = Slant depth
    # Second index = Underground energy

    survival_probability_tensor = np.zeros(
        (len(constants.ENERGIES), len(constants._SLANT_DEPTHS), len(constants.ENERGIES))
    )

    if constants.get_verbose() > 1:
        print("Calculating survival probabilities.")

    # Loop over the surface energies and slant depths
    # Histogram the underground energies to count how many are in each underground energy bin
    # Divide the counts by the number of muons
    # This counts how many muons survived per bin out of the total number that was thrown

    for i in range(len(constants.ENERGIES)):

        for x in range(len(constants._SLANT_DEPTHS)):

            survival_probability_tensor[i, x, :] = np.histogram(
                np.array(u_energies[i, x]), bins=constants._E_BINS
            )[0] / float(constants.get_n_muon())

    if constants.get_verbose() > 1:
        print("Finished calculating survival probabilities.")

    # Write the results to the file

    if output:

        constants._check_directory(
            os.path.join(constants.get_directory(), "survival_probabilities"),
            force=force,
        )

        if file_name == "" or not isinstance(file_name, str):

            file_name = os.path.join(
                constants.get_directory(),
                "survival_probabilities",
                "{0}_{1}_{2}_Survival_Probabilities.txt".format(
                    constants.get_medium(),
                    constants.get_density(),
                    constants.get_n_muon(),
                ),
            )

        file_out = open(file_name, "w")

        for i in range(len(constants.ENERGIES)):

            for x in range(len(constants._SLANT_DEPTHS)):

                for u in range(len(constants.ENERGIES)):

                    file_out.write(
                        "{0:1.14f} {1:1.5f} {2:1.14f} {3:1.14e}\n".format(
                            constants.ENERGIES[i],
                            constants._SLANT_DEPTHS[x],
                            constants.ENERGIES[u],
                            survival_probability_tensor[i, x, u],
                        )
                    )

        file_out.close()

        if constants.get_verbose() > 1:
            print(f"Survival probabilities written to {file_name}.")

    return survival_probability_tensor


def load_survival_probability_tensor_from_file(file_name="", force=False):

    """
    Retrieve a survival probability tensor in units of [(MeV^2 km.w.e.)^-1] stored in data/survival_probabilities based on the set global propagation constants.

    The function searches for a file name that matches the set medium, density, and number of muons. If the file does not exist, prompt the user to run propagation.calc_survival_probability_tensor().

    Parameters
    ----------
    file_name : str, optional (default: constructed from set global propagation constants)
        The name of the file in which to store the results. If output is False or None, this is ignored.

    force : bool
        If True, force the calculation of a new survival probability tensor and the creation of new directories if required.

    Returns
    -------
    survival_probability_tensor : NumPy ndarray
        A three-dimensional array containing the survival probabilities. The shape of the array will be (91, 28, 91), and the survival probabilities will be in units of [(MeV^2 km.w.e.)^-1].
    """

    # Check values

    constants._check_constants(force=force)

    # Check if a survival probability tensor has already been loaded
    # If so, return the tensor that has already been loaded
    # If not, load a tensor based on the set global propagation constants

    if (
        constants._current_survival_probability_tensor is not None
        and constants._survival_probability_tensor_configuration
        == {
            "medium": constants.get_medium(),
            "density": constants.get_density(),
            "n_muon": constants.get_n_muon(),
        }
    ):

        if constants.get_verbose() > 1:
            print(
                "Survival probabilities already loaded for {0} with density {1} gcm^-3 and {2} muons.".format(
                    constants._survival_probability_tensor_configuration["medium"],
                    constants._survival_probability_tensor_configuration["density"],
                    constants._survival_probability_tensor_configuration["n_muon"],
                )
            )

        return constants._current_survival_probability_tensor

    # Define a function to update the survival probability tensor cache

    def update_cache(survival_probability_tensor):

        constants._survival_probability_tensor_configuration = {
            "medium": constants.get_medium(),
            "density": constants.get_density(),
            "n_muon": constants.get_n_muon(),
        }
        constants._current_survival_probability_tensor = survival_probability_tensor

    # Define a function to run if there is no survival probability file

    def no_file(force):

        # If the file does not exist, ask the user if they want to run PROPOSAL to create it

        if not force:

            answer = input(
                "No survival probability matrix currently exists for the set medium, density, or number of muons. Would you like to create one (y/n)?: "
            )

        if force or answer.lower() == "y":

            survival_probability_tensor = calc_survival_probability_tensor(force=force)

            update_cache(survival_probability_tensor)

            return survival_probability_tensor

        else:

            print("Survival probabilities not calculated.")

            return

    # Construct a file name based on the set medium, density, and number of muons if one has not been provided

    if file_name == "" or not isinstance(file_name, str):

        file_name = os.path.join(
            constants.get_directory(),
            "survival_probabilities",
            "{0}_{1}_{2}_Survival_Probabilities.txt".format(
                constants.get_medium(), constants.get_density(), constants.get_n_muon()
            ),
        )

    # Check if the file exists

    if os.path.isfile(file_name):

        if constants.get_verbose() > 1:
            print(f"Loading survival probabilities from {file_name}.")

        # Check that the file has the correct numbers of energies and slant depths

        file = open(file_name, "r")
        n_lines = len(file.read().splitlines())
        file.close()

        # If so, read in the survival probabilities in from it

        if n_lines == len(constants.ENERGIES) * len(constants._SLANT_DEPTHS) * len(
            constants.ENERGIES
        ):

            survival_probability_tensor = np.reshape(
                np.loadtxt(file_name)[:, 3],
                (
                    len(constants.ENERGIES),
                    len(constants._SLANT_DEPTHS),
                    len(constants.ENERGIES),
                ),
            )

            if constants.get_verbose() > 1:
                print("Loaded survival probabilities.")

            update_cache(survival_probability_tensor)

            return survival_probability_tensor

        # If the file does not have the correct numbers of energies and slant depths, run the no_file() function

        else:

            return no_file(force=force)

    # If the file does not exist, run the no_file() function

    else:

        return no_file(force=force)
