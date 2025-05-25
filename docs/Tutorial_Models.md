# **Tutorial - Changing the Models**

## Table of Contents

1. [Primary Model](#primary-model)
2. [Interaction Model](#interaction-model)
3. [Atmospheric Model](#atmospheric-model)

## Primary Model

The default primary cosmic ray flux model is GlobalSplineFitBeta, set with the ``primary_model`` keyword argument as ``"gsf"``. The following primary models are available to be set with a string:

| Model                            | How to Use in MUTE    | Alternatively                          | Reference File |
|:---------------------------------|:----------------------|:---------------------------------------|:---------------|
| **GlobalSplineFitBeta**          | ``"gsf"``             | ``(pm.GlobalSplineFitBeta, None)``     | ``"surface_fluxes_USStd_None_sibyll23c_gsf.txt"``
| **HillasGaisser2012 (H3a)**      | ``"hg"`` or ``"h3a"`` | ``(pm.HillasGaisser2012, "H3a")``      | ``"surface_fluxes_USStd_None_sibyll23c_h3a.txt"``
| **HillasGaisser2012 (H4a)**      | ``"h4a"``             | ``(pm.HillasGaisser2012, "H4a")``      | ``"surface_fluxes_USStd_None_sibyll23c_h4a.txt"``
| **GaisserHonda**                 | ``"gh"``              | ``(pm.GaisserHonda, None)``            | ``"surface_fluxes_USStd_None_sibyll23c_gh.txt"``
| **GaisserStanevTilav (3-gen)**   | ``"gst3"``            | ``(pm.GaisserStanevTilav, "3-gen")``   | ``"surface_fluxes_USStd_None_sibyll23c_gst3.txt"``
| **GaisserStanevTilav (4-gen)**   | ``"gst4"``            | ``(pm.GaisserStanevTilav, "4-gen")``   | ``"surface_fluxes_USStd_None_sibyll23c_gst4.txt"``
| **ZatsepinSokolskaya (Default)** | ``"zs"``              | ``(pm.ZatsepinSokolskaya, "default")`` | ``"surface_fluxes_USStd_None_sibyll23c_zs.txt"``
| **ZatsepinSokolskaya (PAMELA)**  | ``"zsp"``             | ``(pm.ZatsepinSokolskaya, "pamela")``  | ``"surface_fluxes_USStd_None_sibyll23c_zsp.txt"``
| **SimplePowerlaw27**             | ``"pl27"``            | ``(pm.SimplePowerlaw27, None)``        | ``"surface_fluxes_USStd_None_sibyll23c_pl27.txt"``

Note that the default primary model is ``"sibyll23c"``, though this can be changed (see [below](#interaction-model)).

To calculate underground intensities for GaisserHonda, for example, one can do:

```python
mtu.calc_u_intensities(method = "sd", primary_model = "gh")
```

Alternatively, in the ``mts.calc_s_fluxes()`` function, the primary model may be set using a tuple. This gives access to the rest of the models available in MCEq. For example:

```python
import crflux.models as pm

s_fluxes = mts.calc_s_fluxes(primary_model = (pm.GaisserStanevTilav, "3-gen"))

mtu.calc_u_intensities(method = "sd", s_fluxes = s_fluxes)
```

This option is only available in the ``mts.calc_s_fluxes()`` function. The other loading and calculation functions require the primary model to be specified with one of the strings in the list above, as they will search for files with names that contain the strings.

For more information, see the <a href="https://mceq.readthedocs.io/en/latest/tutorial.html#changing-cosmic-ray-flux-model" target="_blank">MCEq Documentation</a> and the <a href="https://crfluxmodels.readthedocs.io/en/latest/index.html">crflux Documentation</a>. For an example, see [``/examples/example_primary_flux_models.ipynb``](../examples/example_primary_flux_models.ipynb).

## Interaction Model

The default hadronic interaction model is SIBYLL-2.3c, set with the ``interaction_model`` keyword argument as ``"sibyll23c"``. The following hadronic interaction models are available:

| Model                          | How to Use in MUTE      | Reference File |
|:-------------------------------|:------------------------|:---------------|
| **DDM**                        | ``"ddm"``               | ``"surface_fluxes_USStd_None_ddm_gsf.txt"``
| **DDM Positive Error**         | ``"ddm_err_pos"``       | ``"surface_fluxes_USStd_None_ddm_err_pos_gsf.txt"``
| **DDM Negative Error**         | ``"ddm_err_neg"``       | ``"surface_fluxes_USStd_None_ddm_err_neg_gsf.txt"``
| **SIBYLL-2.3d**                | ``"sibyll23d"``         | ``"surface_fluxes_USStd_None_sibyll23d_gsf.txt"``
| **SIBYLL-2.3d Positive Error** | ``"sibyll23d_err_pos"`` | ``"surface_fluxes_USStd_None_sibyll23d_err_pos_gsf.txt"``
| **SIBYLL-2.3d Negative Error** | ``"sibyll23d_err_neg"`` | ``"surface_fluxes_USStd_None_sibyll23d_err_neg_gsf.txt"``
| **SIBYLL-2.3c**                | ``"sibyll23c"``         | ``"surface_fluxes_USStd_None_sibyll23c_gsf.txt"``
| **SIBYLL-2.3c03**              | ``"sibyll23c03"``       | ``"surface_fluxes_USStd_None_sibyll23c03_gsf.txt"``
| **SIBYLL-2.3pp**               | ``"sibyll23pp"``        | ``"surface_fluxes_USStd_None_sibyll23pp_gsf.txt"``
| **SIBYLL-2.3**                 | ``"sibyll23"``          | ``"surface_fluxes_USStd_None_sibyll23_gsf.txt"``
| **SIBYLL-2.1**                 | ``"sibyll21"``          | ``"surface_fluxes_USStd_None_sibyll21_gsf.txt"``
| **EPOS-LHC**                   | ``"eposlhc"``           | ``"surface_fluxes_USStd_None_eposlhc_gsf.txt"``
| **QGSJet-II-04**               | ``"qgsjetii04"``        | ``"surface_fluxes_USStd_None_qgsjetii04_gsf.txt"``
| **QGSJet-II-03**               | ``"qgsjetii03"``        | ``"surface_fluxes_USStd_None_qgsjetii03_gsf.txt"``
| **QGSJet-01c**                 | ``"qgsjet01c"``         | ``"surface_fluxes_USStd_None_qgsjet01c_gsf.txt"``
| **DPMJET-III-3.0.6**           | ``"dpmjetiii306"``      | ``"surface_fluxes_USStd_None_dpmjetiii306_gsf.txt"``
| **DPMJET-III-19.1**            | ``"dpmjetiii191"``      | ``"surface_fluxes_USStd_None_dpmjetiii191_gsf.txt"``

Note that the uncertainties for ``"ddm"`` and ``"sibyll23d"`` are only available using the files provided by MUTE from GitHub for the default primary model ``"gsf"``; MCEq cannot calculate new matrices for these uncertainties. Additionally, new ``"ddm"`` matrices cannot be calculated through MUTE, but can be generated using MCEq (see [``DDM_example.ipynb``](https://github.com/mceq-project/mceq-examples/blob/main/DDM_example.ipynb)) and passed into any ``mute.underground`` function using the ``s_fluxes`` parameter. Note as well that the default primary model is ``"gsf"``, though this can be changed (see [above](#primary-model))

To calculate underground intensities for EPOS-LHC, for example, one can do:

```python
mtu.calc_u_intensities(method = "sd", interaction_model = "eposlhc")
```

For more information, see the <a href="https://mceq.readthedocs.io/en/latest/tutorial.html#changing-hadronic-interaction-models" target="_blank">MCEq Documentation</a>.

## Atmospheric Model

The default atmosphere is <a href="https://ntrs.nasa.gov/citations/19770009539" target="_blank">US Standard Atmosphere</a>. Using the default is equivalent to running:

```python
mtu.calc_u_intensities(method = "sd", model = "mceq", atmosphere = "corsika", location = "USStd", month = None)
```

In order to calculate underground fluxes and intensities (and thus the surface fluxes) at a given location, ``atmosphere`` must be set to ``"msis00"`` (the <a href="https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2002JA009430" target="_blank">NRLMSISE-00 Model</a> is used to calculate seasonal variations of the atmsophere). ``month`` must also be set to one of the months of the year as a string (``"January"``, ``"February"``, etc.).

Note that the MCEq models ``"ddm"`` and ``"sibyll23d"`` are not available for different atmospheric models; the most recent working interaction model is ``"sibyll23c"``.

The following locations are available:

| Location                   | How to Use in MUTE | Coordinates |
|:---------------------------|:-------------------|:------------|
| **US Standard Atmosphere** | ``"USStd"``        | N/A
| **South Pole**             | ``"SouthPole"``    | ``(-90.00, 0.0)``
| **Karlsruhe**              | ``"Karlsruhe"``    | ``(49.00, 8.4)``
| **Geneva**                 | ``"Geneva"``       | ``(46.20, 6.1)``
| **Tokyo**                  | ``"Tokyo"``        | ``(35.00, 139.0)``
| **Gran Sasso**             | ``"SanGrasso"``    | ``(42.40, 13.5)``
| **Tel Aviv**               | ``"TelAviv"``      | ``(32.10, 34.8)``
| **Kennedy Space Centre**   | ``"KSC"``          | ``(32.10, -80.7)``
| **Soudan Mine**            | ``"SoudanMine"``   | ``(47.80, -92.2)``
| **Tsukuba**                | ``"Tsukuba"``      | ``(36.20, 140.1)``
| **Lynn Lake**              | ``"LynnLake"``     | ``(56.90, -101.1)``
| **Peace River**            | ``"PeaceRiver"``   | ``(56.15, -117.2)``
| **Ft Sumner**              | ``"FtSumner"``     | ``(34.50, -104.2)``
| **Lake Baikal**            | ``"LakeBaikal"``   | ``(51.60, 103.91)``
| **P-ONE**                  | ``"P-ONE"``        | ``(47.90, -127.7)``
| **KM3NeT-ARCA**            | ``"KM3NeT-ARCA"``  | ``(36.6827, 15.1322)``
| **KM3NeT-ORCA**            | ``"KM3NeT-ORCA"``  | ``(42.80, 6.0)``
| **TRIDENT**                | ``"TRIDENT"``      | ``(17.30, 114.0)``

Note that because these location strings are passed into MCEq directly, they are case-sensitive, unlike other model parameters.

As of [v3.0.0](https://github.com/wjwoodley/mute/releases/tag/3.0.0), arbitrary locations can also be passed into calculation functions using a tuple for the ``location`` argument. The following are two equivalent ways of calculating underground intensities for Tokyo in July:

```python
mtu.calc_u_intensities(method = "sd", model = "mceq", atmosphere = "msis00", location = "Tokyo", month = "July")
mtu.calc_u_intensities(method = "sd", model = "mceq", atmosphere = "msis00", location = (35.00, 139.0), month = "July")
```

Additional location strings specified by ``(longitude, latitude, altitude)`` coordinates can also be added by editing the MCEq source code at line <a href="https://github.com/mceq-project/MCEq/blob/ff1c56ab02f7eafaee2dcb5eb6263e3d1978c6bf/MCEq/geometry/nrlmsise00_mceq.py#L29" target="_blank">29 of `nrlmsise00_mceq.py`</a> and lines <a href="https://github.com/mceq-project/MCEq/blob/ff1c56ab02f7eafaee2dcb5eb6263e3d1978c6bf/MCEq/geometry/density_profiles.py#L1438" target="_blank">1438</a> and <a href="https://github.com/mceq-project/MCEq/blob/ff1c56ab02f7eafaee2dcb5eb6263e3d1978c6bf/MCEq/geometry/density_profiles.py#L642" target="_blank">642 of `density_profiles.py`</a>.

Month names are given in the variables ``mtc.MONTHS`` and ``mtc.MONTHS_SNAMES``, which are, respectively:

```python
["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]
```
```python
["Jan.", "Feb.", "Mar.", "Apr.", "May.", "Jun.", "Jul.", "Aug.", "Sep.", "Oct.", "Nov.", "Dec."]
```

For an example of calculations of seasonal variations at the surface and underground, see [``/examples/example_seasonal_variations.ipynb``](../examples/example_seasonal_variations.ipynb).