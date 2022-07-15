# Changelog

## [2.0.0](https://github.com/wjwoodley/mute/releases/tag/2.0.0) - 14 July 2022

### **New Features**

* **Mountains:** Calculations can now be done for non-flat overburdens. The ``calc_u_fluxes()`` function returns a three-dimensional tensor when the overburden type is set to ``"mountain"`` (or when the new ``full_tensor`` argument is set to ``True``). [An example](examples/example_mountain_calculations.ipynb) has been added to show this.
* **Intensity Functions:** Intensities are now calculated with the ``mtu.calc_u_intensities()`` function by passing the ``method`` argument. The previous functions, ``mtu.calc_u_intensities_tr()`` and ``mtu.calc_u_intensities_eq()`` have been removed.
* **Surface Calculations:** The functions ``mts.calc_s_intensities()`` and ``mts.calc_s_tot_fluxes()`` have been added for calculations of surface intensities and total fluxes.
* **Calculations from Pre-Calculated Matrices:** Underground fluxes, intensities, and total fluxes can now be calculated from surface flux matrices and survival probability tensors that are already defined in the code by passing them directly into the function with the ``s_fluxes`` and ``survival_probability_tensor`` arguments. The old interface of specifying the models and atmospheres is also still available.
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
* **Energy Grid:** The energy grid can now only be accessed with ``mtc._ENERGIES``, not ``mtc.energies``.
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