# **Tutorial - Modelling Labs**

## Table of Contents

1. [Laboratory Locations](#laboratory-locations)
1. [Rock Types](#rock-types)
2. [Mountain Maps](#mountain-maps)

## Laboratory Locations

The following table provides the locations of the labs used for the total flux plot in [``/examples/example_total_flux_plot.ipynb``](../examples/example_total_flux_plot.ipynb) for reference. These coordinates can be passed to the ``location`` parameter when the ``atmosphere`` parameter is set to ``msis00`` (see [``Tutorial_Models``](Tutorial_Models.md#atmospheric-model)) to specify the atmosphere for the specific location of the lab.

| Laboratory  | Coordinates        | How to Use in MUTE |
|:------------|:-------------------|:-------------------|
| **WIPP**    | (32.372, -103.794) | ``location = (32.372, -103.794)``
| **Y2L**     | (38.010, 128.543)  | ``location = (38.010, 128.543)``
| **Soudan**  | (47.823, -92.237)  | ``location = (47.823, -92.237)``
| **Kamioka** | (36.423, 137.315)  | ``location = (36.423, 137.315)``
| **Boulby**  | (54.553, -0.825)   | ``location = (54.553, -0.825)``
| **SUPL**    | (-37.070, 142.810) | ``location = (-37.070, 142.810)``
| **LNGS**    | (42.400, 13.500)   | ``location = (42.400, 13.500)``
| **LSM**     | (45.179, 6.689)    | ``location = (45.179, 6.689)``
| **SURF**    | (44.353, -103.744) | ``location = (44.353, -103.744)``
| **SNOLAB**  | (46.472, -81.187)  | ``location = (46.472, -81.187)``
| **CJPL-I**  | (28.153, 101.711)  | ``location = (28.153, 101.711)``

## Rock Types

New in MUTE [v3.0.0](https://github.com/wjwoodley/mute/releases/tag/3.0.0) is the ability to model various types of rock in addition to standard rock. Rock types can be set by setting the reference density with the ``mtc.set_density()`` function to a specific value. The density is associated with $\langle Z\rangle$ and $\langle A\rangle$ values. The available rock types for given laboratories are listed in the table below.

| Lab or Medium     | Density | $\langle Z\rangle$ | $\langle A\rangle$ | Medium        | How to Use in MUTE | Reference File |
|:------------------|:--------|:-------------------|:-------------------|:--------------|:-------------------|:---------------|
| **Standard Rock** | 2.65    | 11.00              | 22.00              | Rock [<a href="https://github.com/tudo-astroparticlephysics/PROPOSAL/blob/c3b01b7a7006674c9f83db8ea20f882f7169029b/src/PROPOSAL/detail/PROPOSAL/medium/Medium.cxx#L309" target="_blank">Link</a>] | ``mtc.set_medium("rock")``<br>``mtc.set_reference_density(2.65)`` | ``"rock_2.65_1000000_survival_probabilities.npy"`` |
| **WIPP**          | 2.3     | 14.00              | 29.25              | Salt [<a href="https://github.com/tudo-astroparticlephysics/PROPOSAL/blob/c3b01b7a7006674c9f83db8ea20f882f7169029b/src/PROPOSAL/detail/PROPOSAL/medium/Medium.cxx#L279" target="_blank">Link</a>] | ``mtc.set_medium("salt")``<br>``mtc.set_reference_density(2.3)`` | ``"salt_2.3_1000000_survival_probabilities.npy"`` |
| **Y2L** [1]       | 2.7     | 11.79              | 23.79              | Y2L Rock      | ``mtc.set_medium("y2l_rock")``<br>``mtc.set_reference_density(2.7)`` | ``"y2l_rock_2.7_1000000_survival_probabilities.npy"`` |
| **Soudan**        | 2.85    | 12.32              | 24.90              | Rock          | `mtc.set_medium("rock")`<br>`mtc.set_reference_density(2.85)` | ``"rock_2.85_1000000_survival_probabilities.npy"`` |
| **Kamioka**       | 2.7     | 11.31              | 22.76              | Rock          | ``mtc.set_medium("rock")``<br>``mtc.set_reference_density(2.7)`` | ``"rock_2.7_1000000_survival_probabilities.npy"`` |
| **Boulby**        | 2.62    | 11.70              | 23.60              | Rock          | ``mtc.set_medium("rock")``<br>``mtc.set_reference_density(2.62)`` | ``"rock_2.62_1000000_survival_probabilities.npy"`` |
| **SUPL** [2]      | -       | -                  | -                  | Rock          | ``mtc.set_medium("rock")``<br>``mtc.set_reference_density(2.65)`` | ``"rock_2.65_1000000_survival_probabilities.npy"`` |
| **LNGS**          | 2.72    | 11.42              | 22.83              | Rock          | ``mtc.set_medium("rock")``<br>``mtc.set_reference_density(2.72)`` | `"rock_2.72_1000000_survival_probabilities.npy"` |
| **LSM**           | 2.73    | 11.74              | 23.48              | Fréjus Rock [<a href="https://github.com/tudo-astroparticlephysics/PROPOSAL/blob/c3b01b7a7006674c9f83db8ea20f882f7169029b/src/PROPOSAL/detail/PROPOSAL/medium/Medium.cxx#L321" target="_blank">Link</a>]| ``mtc.set_medium("frejus_rock")``<br>``mtc.set_reference_density(2.73)`` | ``"frejus_rock_2.73_1000000_survival_probabilities.npy"`` |
| **SURF**          | 2.86    | 12.01              | 23.98              | Rock          | ``mtc.set_medium("rock")``<br>``mtc.set_reference_density(2.86)`` | ``"rock_2.86_1000000_survival_probabilities.npy"`` |
| **SNOLAB**        | 2.83    | 12.02              | 24.22              | Rock          | ``mtc.set_medium("rock")``<br>``mtc.set_reference_density(2.83)`` | ``"rock_2.83_1000000_survival_probabilities.npy"`` |
| **CJPL-I**        | 2.8     | 12.15              | 24.30              | Rock          | ``mtc.set_medium("rock")``<br>``mtc.set_reference_density(2.8)`` | ``"rock_2.8_1000000_survival_probabilities.npy"`` |
| **Water or Ice**  | 0.997   | 7.33               | 14.43              | Water [<a href="https://github.com/tudo-astroparticlephysics/PROPOSAL/blob/c3b01b7a7006674c9f83db8ea20f882f7169029b/src/PROPOSAL/detail/PROPOSAL/medium/Medium.cxx#L237" target="_blank">Link</a>] | ``mtc.set_medium("water")``<br>``mtc.set_reference_density(0.997)`` | ``"water_0.997_1000000_survival_probabilities.npy"`` |
| **ANTARES Water** | 1.03975 | 7.42               | 14.76              | ANTARES Water [<a href="https://github.com/tudo-astroparticlephysics/PROPOSAL/blob/c3b01b7a7006674c9f83db8ea20f882f7169029b/src/PROPOSAL/detail/PROPOSAL/medium/Medium.cxx#L426" target="_blank">Link</a>] | ``mtc.set_medium("antares_water")``<br>``mtc.set_reference_density(1.03975)`` | ``"antares_water_1.03975_1000000_survival_probabilities.npy"`` |

```
[1]: Because the rocks above Y2L and Kamioka both have an average density of 2.7 gcm^-3, Y2L rock must be
     specified by "y2l_rock" rather than "rock" like all other rock types. This does not imply a different
     set of Sternheimer parameters; "y2l_rock" uses the Sternheimer parameters of standard rock.

[2]: Because <Z> and <A> values are not available for SUPL, the transfer tensor for standard rock is used.
```

Setting the medium as described above tells MUTE to reference the listed reference file, which contains a pre-calculated and supplied survival probability tensor that was calculated for the specified $\langle Z\rangle$ and $\langle A\rangle$ values. These values for each medium as listed in Table II of <https://inspirehep.net/literature/2799258> and in the table above were set in the PROPOSAL source code for component defintions in line <a href="https://github.com/tudo-astroparticlephysics/PROPOSAL/blob/c3b01b7a7006674c9f83db8ea20f882f7169029b/src/PROPOSAL/detail/PROPOSAL/medium/Components.cxx#L285" target="_blank">285 of ``Components.cxx``</a> in order to generate these transfer tensors.

Details on chemical composition and Sternheimer parameters for each medium as listed in Table VI of <https://inspirehep.net/literature/2799258> are given in the PROPOSAL source code for medium definitions <a href="https://github.com/tudo-astroparticlephysics/PROPOSAL/blob/c3b01b7a7006674c9f83db8ea20f882f7169029b/src/PROPOSAL/detail/PROPOSAL/medium/Medium.cxx" target="_blank">here</a>.

## Mountain Maps

Also new in MUTE [v3.0.0](https://github.com/wjwoodley/mute/releases/tag/3.0.0) is the supply of mountain map files for multiple laboratories under mountains. Each mountain map is centered around a specific detector in the given lab in order to provide an origin to the coordinate system used in the map. Maps are typically produced from satallite data centered on the (latitude, longitude) coordinates of the detector, or from detector data. As MUTE requires the mountain map files to be in a specific format in order to be compatible with the ``mtc.load_mountain()`` function, and as some maps are difficult to find or not available in the literature, this is a convenient way of quickly loading the maps for calculations for one or multiple laboratories. The available maps and how to load them are listed in the table below.

| Lab | Detector | How to Use in MUTE | Reference File | Citations |
|:----|:---------|:-------------------|:---------------|:----------|
| **Y2L** [1] | COSINE-100       |``mtc.set_medium("y2l_rock")``<br>``mtc.set_reference_density(2.7)``<br>``mtc.load_mountain("Y2L")``|``y2l_mountain.txt``| |
| **Kamioka** | Super-Kamiokande |``mtc.set_medium("rock")``<br>``mtc.set_reference_density(2.7)``<br>``mtc.load_mountain("SuperK")``|``superk_mountain.txt``|<a href="https://inspirehep.net/literature/824640" target="_blank">S. Abe et al. (KamLAND), Phys. Rev. C <b>81</b>, 025807 (2010).</a><br><br><a href="https://inspirehep.net/literature/607144" target="_blank">Y. Fukuda et al. (Super-Kamiokande), Nucl. Instrum. Meth. A <b>501</b>, 418 (2003).</a><br><br>Digital Map 50 m Grid (Elevation), Geographical Survey Institute of Japan (1997). |
| **Kamioka** | KamLAND          |``mtc.set_medium("rock")``<br>``mtc.set_reference_density(2.7)``<br>``mtc.load_mountain("KamLAND")``|``kamland_mountain.txt``|<a href="https://inspirehep.net/literature/824640" target="_blank">S. Abe et al. (KamLAND), Phys. Rev. C <b>81</b>, 025807 (2010).</a><br><br>Digital Map 50 m Grid (Elevation), Geographical Survey Institute of Japan (1997). |
| **LNGS**    | LVD              |``mtc.set_medium("rock")``<br>``mtc.set_reference_density(2.72)``<br>``mtc.load_mountain("LNGS")``|``lngs_mountain.txt``|<a href="https://inspirehep.net/literature/471316" target="_blank">M. Aglietta et al. (LVD), Phys. Rev. D <b>58</b>, 092005 (1998).</a>|
| **LSM**     | Fréjus           |``mtc.set_medium("frejus_rock")``<br>``mtc.set_reference_density(2.73)``<br>``mtc.load_mountain("LSM")``|``lsm_mountain.txt``|<a href="https://inspirehep.net/literature/26080" target="_blank">C. Berger et al. (Fréjus), Phys. Rev. D <b>40</b>, 2163 (1989).</a> |
| **CJPL-I**  | JNE              |``mtc.set_medium("rock")``<br>``mtc.set_reference_density(2.8)``<br>``mtc.load_mountain("CJPL")``|``cjpl_mountain.txt``|<a href="https://inspirehep.net/literature/1809695" target="_blank">Z. Guo et al. (JNE), Chin. Phys. C <b>45</b>, 025001 (2021).</a> |

```
[1]: The slant depth values for the first 27 azimuthal bins of theta = 0 in this map have been averaged to
     a value of 1.67094153 km.w.e. (previously a minimum of 0.66349473 km.w.e. and a maximum of
     0.68703341 km.w.e.) due to issues with these bins. This change lowers the total underground flux by
     0.02%, which is now (4.72379 ± 0.10664)e-7 cm^-2 s^-1, while it was previously
     (4.72495 ± 0.10665)e-7 cm^-2 s^-1.
```

We request on behalf of the collaborations that have kindly shared their maps that you please take care to cite the citations provided in the table if you make use of any of these maps in your work.

Note that, for these given mountain maps, the string is case-insensitive, meaning ``mtc.load_mountain("Y2L")`` and ``mtc.load_mountain("y2l")`` will load the same map. For user-provided paths, the string is not case-insensitive.