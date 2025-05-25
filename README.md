# MUTE

[![DOI](https://zenodo.org/badge/DOI/10.5282/zenodo.5791812.svg)](https://doi.org/10.5281/zenodo.5791812)

MUTE (**MU**on in**T**ensity cod**E**) is a computational tool for calculating atmospheric muon fluxes and intensities underground. It makes use of the state-of-the-art codes [daemonflux](https://github.com/mceq-project/daemonflux) and [MCEq](https://github.com/afedynitch/MCEq), to calculate surface fluxes, and [PROPOSAL](https://github.com/tudo-astroparticlephysics/PROPOSAL), to simulate the propagation of muons through rock and water.

## Installation

MUTE can be installed via pip:

```
pip install mute
```

This will install most of the requirements, including MCEq.

### Additional Requirements

In order to generate custom survival probability tensors, PROPOSAL should be installed. Because it requires compilation, it needs to be installed separately (see detailed installation instructions [here](https://github.com/tudo-astroparticlephysics/PROPOSAL/blob/master/INSTALL.md)):

```
pip install proposal==7.6.2
```

MUTE has been tested with PROPOSAL v7.6.2; earlier or later versions are not guaranteed to work. Environment set-up help for PROPOSAL installation on a Mac can be found [here](docs/Installing_PROPOSAL_on_a_Mac.md).

## Getting Started

![Total Underground Flux vs Vertical Depth](docs/total_underground_flux.png)

The code for the above plot is found in [``/examples/example_total_flux_plot.ipynb``](examples/example_total_flux_plot.ipynb).

For a basic example and a detailed description of how to use MUTE, see the [Tutorial](docs/Tutorial.md). For further examples, see the [examples](examples).

## Citation

Please cite https://inspirehep.net/literature/1927720 for MUTE v1-2 and https://inspirehep.net/literature/2799258 for MUTE v3.

The current version of MUTE can be cited with the [Zenodo DOI](https://zenodo.org/record/5791812).

The citations for the models and propagation tools used by MUTE can be found in the [MCEq](https://github.com/afedynitch/MCEq#please-cite-our-work) and [PROPOSAL](https://github.com/tudo-astroparticlephysics/PROPOSAL#how-to-cite-proposal) documentation.

### Authors

[William Woodley](mailto:wwoodley@ualberta.ca)

### Contributors

Anatoli Fedynitch ([@afedynitch](https://github.com/afedynitch)) and Marie-Cécile Piro ([@MarieCecile](https://github.com/MarieCecile))

## License

MUTE is licensed under the BSD 3-Clause License (see [LICENSE](LICENSE)).