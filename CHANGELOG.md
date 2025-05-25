# **Changelog**

## [3.0.0](https://github.com/wjwoodley/mute/releases/tag/3.0.0) - 24 May 2025

### **New Features**

* **daemonflux:** [daemonflux](https://github.com/mceq-project/daemonflux) can now be used from within the MUTE framework to calculate suface fluxes. Because either daemonflux or MCEq can be used, a new argument ``model`` has been added to ``surface`` and ``underground`` functions to specify which is to be used (default: ``mceq``; in the next release, the default will become ``daemonflux``).
* **New Functions:** A number of new functions can now be used.
  * **``mts.calc_s_e_spect()``:** Calculates surface energy spectra.
  * **``mtu.calc_u_ang_dist()``:** Calculates underground zenith or azimuthal angular distributions (for mountain overburdens only, since the distributions are trivial for flat overburdens).
  * **``mtu.calc_u_e_spect()``:** Calculates underground energy spectra.
  * **``mtu.calc_u_mean_e()``:** Calculates underground mean and median energies, as well as upper and lower 68% and 95% confidence intervals for the median.
  * **``mtu.calc_depth()``:** Calculates various depths (equivalent vertical, straight vertical, minimum, maximum, and average) for mountain overburdens.
  * A number of new examples have been added to [``/examples``](examples) to demonstrate use of these functions.
* **Provided Mountain Maps:** Mountain maps for various laboratories are now provided (Y2L, Super-Kamiokande, KamLAND, LSM, LNGS, and CJPL). See [``Tutorial_Labs``](/docs/Tutorial_Labs.md) for information on use.
* **New Propagation Media:** A number of new propagation medium options, corresponding to different underground and underwater laboratories, are now available in addition to standard rock. ANTARES water is now available as a propagation medium. The ice propagation medium is now an alias for (fresh) water. New transfer tensors for these media are provided in the supplied data files. See [``Tutorial_Labs``](/docs/Tutorial_Labs.md) for information on use. Air is no longer available as a propagation medium.
* **Location Coordinates:** Locations for MCEq can now be set by passing in (latitude, longitude) tuples (in degrees) into the ``location`` parameter in any function.
* **Loading Mountains:** New parameters are available when loading mountains, including the ``scale`` parameter, which can be used to scale all slant depths in the loaded mountain map by some ratio (useful for calculating uncertainties, for example).
* **Loading Surface Fluxes with Files:** Surface fluxes can now be loaded from user-named files by passing the file name into the ``file_name`` parameter in the ``mts.load_s_fluxes_from_file()`` function.

### **Bug Fixes**

* **Intensity Calculation:** A bug in the calculation of intensities for mountains was corrected using a new calculation method. The new results differ from those returned by MUTE v2 by less than 7% at the deepest depths.
* **Loading Transfer Tensors:** Transfer tensors can now be loaded successfully when passing a file name to ``mtp.load_survival_probability_tensor_from_file()``.
* **NumPy and SciPy:** MUTE function calls to NumPy and SciPy have been updated to correspond with the latest versions of these libraries.

### **Other Changes**

* **Number of Muons:** The default value for ``n_muon`` is now ``1e6`` instead of ``1e5``. New transfer tensors for 1e6 muons are provided in the supplied data files.
* **Lowercase Letters:** Most elements of all file names, as well as some arguments, are now lowercase. Interaction models also now have their dashes and periods removed if passed into functions with dashes and periods. This is backwards-compatible. For example, user-input ``"SIBYLL-2.3d"`` is now turned into ``"sibyll23d"`` in all cases. Note that the location arguments, including ``"USStd"``, are case-sensitive, as are month arguments, including ``None``, as these arguments are passed directly into MCEq. New file names in lowercase are provided in the supplied data files, and old files have been removed (though will persist in the data directory if MUTE is updated over a previous installation).
* **Vertical Depths for Mountains:** An exception is now thrown when attempting to set a vertical depth while the overburden type is set to mountain. This alerts the user that they are possibly in a different overburden mode from intended.
* **Changing the Overburden Type:** When the overburden type is changed, the set and loaded constants specific to flat or mountain overburdens are now reset (meaning, if, for example, the overburden type is ``"flat"`` and the vertical depth is set to 3.0 km.w.e., then the overburden type is changed to ``"mountain"``, then changed back to ``"flat"``, the vertical depth will be reset to 0.5 km.w.e.).
* **Transfer Tensor Files:** Transfer tensors are now stored in ``.npy`` files instead of plain text files in order to save disk space.
* **Documentation:** The tutorial has now been split up into different files for better organisation and easier navigation.
* **ASCII Art:** The ``__init__.py`` file now prints ASCII art when a MUTE module is imported for the first time.
* **Default Interaction Model:** In the [previous release of MUTE](https://github.com/wjwoodley/mute/releases/tag/2.0.1), it was stated that the default interaction model would be changed from SIBYLL-2.3c to DDM. Instead, the default interaction model remains SIBYLL-2.3c. Although SIBYLL-2.3d and DDM are now available to be used in MCEq, SIBYLL-2.3c is still the most recent model for which all features in MUTE are available, particularly calculating surface fluxes with the NRLMSISE-00 atmosphere in order to specify location and month.

### **Deprecations and Upcoming Changes**

* **Energy Units:** In the next release (v3.1.0), the energy units in all cases will change from [MeV] to [GeV].
* **Setting Densities:** The ``mtc.set_density()`` and ``mtc.get_density()`` functions have been renamed ``mtc.set_reference_density()`` and ``mtc.get_reference_density()`` respectively. The former will be removed in the next release (v3.1.0).
* **Default Surface Flux Model:** In the next release (v3.1.0), the default surface flux model set by the ``model`` parameter will be changed to ``daemonflux`` (it is currently ``"mceq"``).
* **Hillas-Gaisser Model Parameters:** The ``"hg"`` primary model option has been split into ``"h3a"`` and ``"h4a"`` to refer to the three- and four-component Hillas-Gaisser primary flux models. The former is now deprecated and will be removed in the next release (v3.1.0). This is to conform to naming conventions within the cosmic ray community.
* **File Name Parameters:** In the next release (v3.1.0), file name parameters previously called ``file_name`` will be called ``input_file`` and ``output_file`` to remove ambiguity.

## [2.0.1](https://github.com/wjwoodley/mute/releases/tag/2.0.1) - 17 November 2022

### **Bug Fixes**

* **PROPOSAL Version:** The suggested PROPOSAL version has been updated to v7.4.2. The previously-suggested v7.1.1 sometimes returns an error when calculating new transfer tensors.

### **Other Changes**

* **Intensities Method:** The ``method`` argument is now an optional argument in the ``mtu.calc_u_intensities()`` function. The default for flat overburdens is ``"sd"``, and the default for mountains is ``"dd"``. Other methods can still be specified as normal, and this should not change any existing scripts. See the [Tutorial](docs/Tutorial.md#calculating-underground-intensities) or the function docstrings for more information.
* **Default Interaction Model:** In the next release of MUTE, the default interaction model will be changed from SIBYLL-2.3c to DDM.

## [2.0.0](https://github.com/wjwoodley/mute/releases/tag/2.0.0) - 15 July 2022

### **New Features**

* **Mountains:** Calculations can now be done for non-flat overburdens. The ``calc_u_fluxes()`` function returns a three-dimensional tensor when the overburden type is set to ``"mountain"`` (or when the new ``full_tensor`` argument is set to ``True``). [An example](examples/example_mountain_calculations.ipynb) has been added to show this.
* **Intensity Functions:** Intensities are now calculated with the ``mtu.calc_u_intensities()`` function by passing the ``method`` argument. The previous functions, ``mtu.calc_u_intensities_tr()`` and ``mtu.calc_u_intensities_eq()`` have been removed.
* **Surface Calculations:** The functions ``mts.calc_s_intensities()`` and ``mts.calc_s_tot_fluxes()`` have been added for calculations of surface intensities and total fluxes.
* **Calculations from Pre-Calculated Matrices:** Underground fluxes, intensities, and total fluxes can now be calculated from underground flux matrices, surface flux matrices, and survival probability tensors that are already defined in the code by passing them directly into the function with the ``u_fluxes``, ``s_fluxes``, and ``survival_probability_tensor`` arguments. The old interface of specifying the models and atmospheres is also still available.
* **Energy Threshold:** An energy threshold can now be specified by setting the ``E_th`` argument in ``mtu.calc_u_intensities()`` and ``mtu.calc_u_tot_fluxes`` to a value in MeV.
* **File Names:** Custom input and output file names can now be used in all function by specifying the optional ``file_name`` argument. If the argument is not given, the default file name will be used. If ``output`` is set to ``False`` (either globally with ``mtc.set_output(False)`` or in the function call), no output will be written to the file, and the file name will be ignored. The previously-used ``file_name`` argument in ``mtp.calc_survival_probability_tensor()`` has been renamed ``file_name_pattern``.
* **Primary Models:** The SimplePowerlaw27 model can now be set with ``"PL27"``.
* **Months Variables:** Constants called ``MONTHS`` and ``MONTHS_SNAMES`` have been added for ease of looping through the months of the year. [An example](examples/example_seasonal_variations.ipynb) has been added to show this.

### **Bug Fixes**

* **Downloading Data Files:** The zip folder containing the data files that is downloaded from GitHub is now properly deleted after unzipping.
* **Joining Paths:** File paths are now joined using ``os.path.join()``, rather than using string concatenation.
* **Closing Files:** Some underground calculation functions previously did not properly close the output files after writing the results.
* **Survival Probability Cache:** The survival probability tensor is now only reloaded or re-calculated if the global propagation constants are changed, rather than every time a function is called.
* **Default Tracking:** Default tracking is now turned off in MCEq, improving computation time of surface fluxes.

### **Other Changes**

* **Argument Order:** The order of arguments in the surface and underground functions has been changed to make them more sensical and consistent. See the [Tutorial](docs/Tutorial.md) or the docstrings of individual functions for more information.
* **Slant Depths:** The default slant depths now go from 0.5 to 14 km.w.e. instead of 1 to 12 km.w.e. Values for slant depths outside this range are possible to calculate with the default transfer tensors by setting ``mtc.shallow_extrapolation`` to ``True``, but the results are not guaranteed to be good.
* **Underground Fluxes Return:** The return of ``mtu.calc_u_fluxes()`` is now one single tensor, not a tuple of two matrices.
* **Angles in Underground Fluxes:** ``mtu.calc_u_fluxes()`` no longer takes an ``angles`` argument. The tensor is always returned for the default angles given by ``mtc._ANGLES``.
* **Energy Grid:** The energy grid can now only be accessed with ``mtc.ENERGIES``, not ``mtc.energies``.
* **Surface Output Files:** The ``surface_fluxes`` directory in ``mute/data/`` has been renamed ``surface`` for more generality. The surface flux file names now contain the month as well (or ``None``) if no month is set.
* **Print Grid Functions:** The functions that display the grids used in the output files have been removed.

## [1.0.1](https://github.com/wjwoodley/mute/releases/tag/1.0.1) - 24 December 2021

### Bug Fixes

* **Total Flux Function:** Unnecessary lines in the ``mtu.calc_u_tot_flux()`` function which were causing an error have been removed.
* **Vertical Depth Restrictions:** The restrictions in the ``mtc.set_vertical_depth()`` function for the depth to be between 1 km.w.e. and 12 km.w.e. were removed, due to lack of ensured consistency with a loaded survival probability tensor.
* **Surface Flux Test File:** Fixed a bug where the directory was not being set properly to find the test data file.

### Other Changes

* **Total Flux Test File:** A test file [``test_total_flux.py``](mute/tests/test_total_flux.py) has been added to test the calculation of the total fluxes.
* **Zenodo DOI:** The Zenodo DOI has been included in the [README.md](README.md) for citation of different versions of the code.

## [1.0.0](https://github.com/wjwoodley/mute/releases/tag/1.0.0) - 19 December 2021

### New Features

* **Energy Cuts:** Energies above ~100 PeV are now excluded to make the Monte Carlo and calculations more efficient.
* **Primary Models:** Zatsepin-Sokolskaya models can now be set with ``"ZS"`` and ``"ZSP"``. Additional primary models can now be specified by passing a tuple for the ``primary_model`` argument in the ``mts.calc_s_fluxes()`` function, as is done in MCEq.
* **Air:** Air can now be used as a medium of propagation (with ``mtc.set_medium("air")``).
* **File Names:** An argument has been added to specify alternative file name patterns when loading underground energies from a file with ``mtp.calc_survival_probability_tensor()``.
* **Spline:** A spline function has been added to make calculations of surface fluxes more efficient, and to improve calculations of underground fluxes. The surface fluxes are now calculated for a smaller number of zenith angles, making the MCEq calculations quicker as well.

### Bug Fixes

* **Number of Muons:** Fixed consistency issue with definition of ``n_muon`` between ``mtp._load_u_energies_from_files()`` and ``mtp.calc_survival_probability_tensor()``.
* **Data Directory:** The data files from GitHub are now downloaded to the directory set with ``mtc.set_directory()``.
* **Test Files:** Fixed surface fluxes and underground intensity test files.

### Other Changes

* **Integrals:** Integration is now done using the more accurate ``scipy.integrate.simpson()`` function, instead of ``np.trapz()``.
* **Surface Fluxes Test:** A ``test`` argument has been added to ``mts.calc_s_fluxes()`` which runs only three zenith angles when testing the surface fluxes.
* **Disk Storage:** Underground energies are now stored as binary files instead of literal ASCII files, saving disk space and increasing efficiency of both outputting and loading in underground energies (30 MB to 18 MB, and 20 minutes to 4 minutes).
* **Memory:** The underground energies are now loaded in a more memory-efficient way. The intermediate ``mtu.load_u_energies_from_files()`` function no longer needs to be run.