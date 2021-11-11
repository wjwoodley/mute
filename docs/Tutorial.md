# Tutorial

## Short Example

The nearly-minimal code below will calculate an array of true vertical underground intensities using DDM and SIBYLL-2.3d for an array of 23 default slant depths and will plot them against the slant depths.

```
# Import packages

import mute.constants as mtc
import mute.underground as mtu

import matplotlib.pyplot as plt

# Set the constants

mtc.set_overburden("flat")
mtc.set_medium("rock")
mtc.set_density(2.65)

# Calculate true vertical intensities for the default slant depths and atmosphere

intensities_DDM = mtu.calc_u_intensities_tr(interaction_model = "DDM")
intensities_SIBYLL = mtu.calc_u_intensities_tr(interaction_model = "SIBYLL-2.3d")

# Plot the results

plt.plot(mtc.slant_depths, intensities_DDM, color = "blue", label = "DDM")
plt.plot(mtc.slant_depths, intensities_SIBYLL, color = "magenta", label = "SIBYLL-2.3d")
plt.semilogy()
plt.legend()
plt.show()
```

**NOTE:** MUTE currently provides calculations only for labs under flat overburdens. Non-flat overburdens will be implemented in the future.

## Importing all MUTE Modules

The MUTE constants and functions are split between four modules. They can be imported as follows:

```
import mute.constants as mtc
import mute.surface as mts
import mute.propagation as mtp
import mute.underground as mtu
```

When MUTE is imported for the first time, the data files supplied with the GitHub release will be downloaded to the directory MUTE is installed to. These data files contain some surface flux matrices for the default models and atmosphere as well as a survival probability tensor for standard rock and a propagation of 100 000 muons.

## Changing the Constants

The globally-set constants are stored in the ``constants`` module, and are set using setter functions. The following piece of code sets all of the constants to their default values.

```
mtc.set_verbose(2)
mtc.set_output(True)
# mtc.set_directory("mute/data")
mtc.set_lab("Default")
mtc.set_overburden("flat")
mtc.set_vertical_depth(1)
mtc.set_medium("rock")
mtc.set_density(2.65)
mtc.set_n_muon(1000)
```

The docstrings for each function describe what values they can take.

It is usually best to leave the directory as the default value, which is given by ``os.path.join(os.path.dirname(__file__), "data")``. Data files downloaded automatically from GitHub will be stored where MUTE is installed (the location of this directory can be found by running ``pip show mute``), as will output data files created while using MUTE. The directory can be changed for ease of use on a computing cluster (see the example in [``examples/cluster``](../examples/cluster) or the example at the end of this tutorial).

When the vertical depth, h, is set, the values of ``mtc.slant_depths`` and ``mtc.angles`` change according to h = X/cos(θ), where X is the slant depth, and θ is the zenith angle. The default values are stored in ``mtc.SLANT_DEPTHS`` and ``mtc.ANGLES`` respectively. The default slant depths are an array of 23 depths between 1 and 12 km.w.e. going up in steps of 0.5 km.w.e.:

```
array([ 1. ,  1.5,  2. ,  2.5,  3. ,  3.5,  4. ,  4.5,  5. ,  5.5,  6. ,
        6.5,  7. ,  7.5,  8. ,  8.5,  9. ,  9.5, 10. , 10.5, 11. , 11.5,
       12. ])
```

For a lab 3.1 km.w.e. underground (or underwater), the vertical depth can be set with ``mtc.set_vertical_depth(3.1)``. Then, printing ``mtc.slant_depths`` gives:

```
array([ 3.5,  4. ,  4.5,  5. ,  5.5,  6. ,  6.5,  7. ,  7.5,  8. ,  8.5,
        9. ,  9.5, 10. , 10.5, 11. , 11.5, 12. ])
```

The angles in ``mtc.angles`` are set according to these depths as well. MUTE uses the default values to calculate the underground flux matrix, then interpolates to the angles in ``mtc.angles`` (or those specified by the ``angles`` or ``depths`` input parameters in the underground flux or intensity functions).

## Units

The following are the default units used throughout MUTE:

* **Energies:** MeV
* **Angles:** Degrees
* **Depths:** km.w.e.
* **Fluxes:** (cm^2 s sr MeV)^-1
* **Intensities:** (cm^2 s sr)^-1
* **Total Fluxes:** (cm^2 s)^-1
* **Densities:** g cm^-3

## Calculating Underground Intensities

Underground intensities can be calculated using the function ``mtu.calc_u_intensities()``. Calling it like this will use the defaults for all of the optional keyword arguments:

```
mtu.calc_u_intensities(angles = None, location = "USStd", month = None, interaction_model = "SIBYLL-2.3c", primary_model = "GSF", atmosphere = "CORSIKA", output = None, force = False)
```

The values for ``angles`` will be taken from ``mtc.angles``, and the value for ``output`` will be taken from ``mtc.get_output()``.

This will return an array of underground muon intensities for the 23 default angles (or otherwise, depending on the value set for the vertical depth). The angles can be changed by passing angles in degrees (float or array-like) into the function. To return the underground intensities for 100 angles between 0 and 85 degrees, for example, one can do:

```
angles = np.linspace(0, 85, 100)

mtu.calc_u_intensities(angles)
```

The function will interpolate over the default energy and angle grids to get the underground intensities at these angles.

If the surface flux or survival probability matrices that are required to calculate the underground fluxes are not found, MUTE will ask if they should be created. To skip this prompt and force them to be created if they are needed, set ``force`` to ``True``. No matter what ``force`` is set to, if a matrix is output to a file, it will always overwrite any existing file with the same name.

### True Vertical Underground Intensities

True vertical underground intensities are calculated with the function ``mtu.calc_u_intensities_tr()``. The input parameters are the same as those for the ``mtu.calc_u_intensities()`` function with the exception of the first parameter, which is ``depths`` instead of ``angles``. The slant depths can be changed by passing slant depths in km.w.e. (float or array-like) into the function. To return the true vertical underground intensities for 20 slant depths between 3 and 10 km.w.e., for example, one can do:

```
depths = np.linspace(3, 10, 20)

mtu.calc_u_intensities_tr(depths)
```

### Changing the Atmosphere

The default atmosphere is [US Standard Atmosphere](https://ntrs.nasa.gov/citations/19770009539). This is set in the following way:

```
mtu.calc_u_intensities(location = "USStd", month = None, atmosphere = "CORSIKA")
```

In order to calculate the underground intensities (and thus the surface fluxes) at a given location, ``month`` must be set to one of the months of the year as a string (``"January"``, ``"February"``, etc.), and ``atmosphere`` must be set to ``"MSIS00"`` (the [NRLMSISE-00 model](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2002JA009430) is used to calculate seasonal variations of the atmsophere). The following will calculate underground intensities for Tokyo in July:

```
mtu.calc_u_intensities(location = "Tokyo", month = "July", atmosphere = "MSIS00")
```

The following locations are available:

- ``"SouthPole"``
- ``"Karlsruhe"``
- ``"Geneva"``
- ``"Tokyo"``
- ``"SanGrasso"``
- ``"TelAviv"``
- ``"KSC"``
- ``"SoudanMine"``
- ``"Tsukuba"``
- ``"LynnLake"``
- ``"PeaceRiver"``
- ``"FtSumner"``

Additional locations specified by (longitude, latitude, altitude) coordinates can be added by editing the MCEq source code.

### Changing the Models

The default hadronic interaction model (which describes the hadronic interactions in the atmosphere that lead to the creation of pions and kaons, which then decay into muons) is ``"SIBYLL-2.3c"``, and is set with the ``interaction_model`` keyword argument. The following hadronic interaction models are available:

- ``"DDM"``
- ``"SIBYLL-2.3d"``
- ``"SIBYLL-2.3c"``
- ``"SIBYLL-2.3"``
- ``"SIBYLL-2.1"``
- ``"EPOS-LHC"``
- ``"QGSJet-II-04"``
- ``"QGSJet-II-03"``
- ``"QGSJet-01c"``
- ``"DPMJET-III-3.0.6"``
- ``"DPMJET-III-19.1"``
- ``"SIBYLL-2.3c_pp"``

Note that ``"DDM"`` and ``"SIBYLL-2.3d"`` are only available using the files provided on GitHub; MCEq cannot currently calculate new matrices for these models.

To calculate underground intensities for EPOS-LHC, for example, one can do:

```
mtu.calc_u_intensities(interaction_model = "EPOS-LHC")
```

For more information, see the [MCEq Documentation](https://mceq.readthedocs.io/en/latest/tutorial.html#changing-hadronic-interaction-models). While MCEq will remove the dashes and points, in order for MUTE to find the correct file in memory when calculating the underground fluxes, the dashes, dots, and capitalisation must match exactly what was used when the function was first ran.

The default primary cosmic ray flux model is ``"GSF"``, for GlobalSplineFitBeta, and is set with the ``primary_model`` keyword argument. The following primary models are available:

- ``"GSF"``: GlobalSplineFitBeta
- ``"GH"``: GaisserHonda
- ``"HG"``: HillasGaisser2012

To calculate underground intensities for GaisserHonda, for example, one can do:

```
mtu.calc_u_intensities(primary_model = "GH")
```

Other models will be added in the future. For more information, see the [MCEq Documentation](https://mceq.readthedocs.io/en/latest/tutorial.html#changing-cosmic-ray-flux-model).

### Vertical-Equivalent Underground Intensities

Vertical equivalent underground intensities can be calculated using the ``mtu.calc_u_intensities_eq()`` function, which takes the same arguments as the ``mtu.calc_u_intensities()`` function. This function multiplies the intensities calculated with ``mtu.calc_u_intensities()`` by cosine of the input angles and returns the result in an array.

## Calculating Underground Fluxes

Underground fluxes can be calculated using the ``mtu.calc_u_fluxes()`` function. This function is also ran by MUTE when calculating underground intensities. Its arguments are the same as ``mtu.calc_u_intensities()``:

```
mtu.calc_u_fluxes(angles = None, location = "USStd", month = None, interaction_model = "SIBYLL-2.3c", primary_model = "GSF", atmosphere = "CORSIKA", output = None, force = False)
```

When calculating underground fluxes and intensities, ``mtu.calc_u_fluxes()`` calls the functions ``mts.load_s_fluxes_from_file()``, which searches for a surface flux file that has a name matching the set atmosphere and models, and ``mtp.load_survival_probability_tensor_from_file()``, which searches for a survival probability file that has a name matching the set global parameters. If it finds the required files, it will load the surface fluxes and / or survival probabilities from those files. If it does not, it will try to calculate new surface flux and / or survival probability matrices.

The ``mtu.calc_u_fluxes()`` function returns a tuple with two elements. The first element is the underground flux matrix for the angles specified by the ``angles`` parameter in the function. The second element of the tuple is the underground flux matrix for 0 degrees. This is used by MUTE to calculate true vertical intensities, but is otherwise unneeded. Therefore, whenever underground fluxes are being calculated, the function should be indexed with ``[0]`` to get the matrix.

## Calculating Surface Fluxes

[MCEq](https://github.com/afedynitch/MCEq) is used to calculate surface muon fluxes. Surface flux matrices can be calculated using the ``mts.calc_s_fluxes()`` function, which takes the same arguments as ``mtu.calc_u_intensities()``, described above, with the exception of ``angles`` (the angles used are always the default 45 angles given by ``mtc.ANGLES_FOR_S_FLUXES``). The default call is:

```
mts.calc_s_fluxes(location = "USStd", month = None, interaction_model = "SIBYLL-2.3c", primary_model = "GSF", atmosphere = "CORSIKA", output = None, force = False)
```

To calculate surface fluxes for Gran Sasso in January using the HillasGaisser2012 primary flux model, for example, one can do:

```
mts.calc_s_fluxes(location = "SanGrasso", month = "January", primary_model = "HG", atmosphere = "MSIS00")
```

MCEq may return negative fluxes for the highest (EeV-scale) energies, but the incredibly small orders of magnitude of these fluxes makes a negligible difference in any results.

## Calculating Survival Probabilities

The main functions in ``mute.propagation`` used for the propagation of muons through matter and the calculation of survival probability tensors are:

1. ``propagate_muons(seed = 0, job_array_number = 0, output = None, force = False)``
2. ``load_u_energies_from_files(n_job = 1, force = False)``
3. ``calc_survival_probability_tensor(seed = 0, output = None, force = False)``

For general purposes, the third function can be used. It will check if underground (or underwater, etc.) muon energies for the set global constants already exist. If they do, it will load them and calculate survival probabilities. If they do not, it will call ``mtp.propagate_muons()``, which begins the Monte Carlo simulation using [PROPOSAL](https://github.com/tudo-astroparticlephysics/PROPOSAL).

### Using MUTE on a Computing Cluster

For high statistics simulations on a computing cluster, the first function is the most useful. The MUTE code to set up a job array of Monte Carlo simulations for 1000 muons per energy-slant depth bin (per job) can be written as follows:

```
# Import packages

import mute.constants as mtc
import mute.propagation as mtp

import argparse

# Parse the job array number from the command line

parser = argparse.ArgumentParser()
parser.add_argument("a", type = int)
args = parser.parse_args()

# Set the constants

mtc.set_verbose(0)
mtc.set_output(True)
mtc.set_directory("mute/data")
mtc.set_lab("Example")
mtc.set_n_muon(1000)

# Propagate the muons

mtp.propagate_muons(seed = args.a, job_array_number = args.a, force = True)
```

To prevent any unnecessary output, the verbosity can be set to ``0`` (though it might be a good idea to keep it at the default ``2`` to help figure out what went wrong if the job fails). To make it easier to access and ``sftp`` out the output Monte Carlo underground energy files, the directory can be changed to the working directory of the code or elsewhere. The ``force`` parameter in ``mtc.propagate_muons()`` should be set to ``True`` to force the creation of any required directories or files so the program does not hang, waiting for input until the job times out.

The job array number can be read in from the command line using ``argparse`` and the seed can be changed in the propagation function to ensure the Monte Carlo results are different for each job. If the job array consisted of 100 jobs, the survival probabilities can then be calculated as below, setting the ``n_job`` input parameter to ``100``:

```
import mute.constants as mtc
import mute.propagation as mtp

mtc.set_lab("Example")
mtc.set_n_muon(1000)

mtp.load_u_energies_from_files(n_job = 100)
mtp.calc_survival_probability_tensor()
```

Here, ``seed`` does not need to be set in ``mtp.calc_survival_probability_tensor()`` because this function will not invoke ``mtp.propagate_muons()``, since the underground energies have already been loaded.
