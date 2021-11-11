# MUTE

MUTE (**MU**on in**T**ensity cod**E**) is a computational tool for calculating atmospheric muon fluxes and intensities underground. It makes use of the state-of-the-art codes [MCEq](https://github.com/afedynitch/MCEq), to calculate surface fluxes, and [PROPOSAL](https://github.com/tudo-astroparticlephysics/PROPOSAL), to simulate the propagation of muons through rock and water.

## Installation

MUTE can be installed via pip:

```
pip install mute
```

This will install most of the requirements, including MCEq.

### Additional Requirements

In order to generate custom survival probability tensors, PROPOSAL should be installed. Because it requires compilation, it needs to be installed separately (see detailed installation instructions [here](https://github.com/tudo-astroparticlephysics/PROPOSAL/blob/master/INSTALL.md)):

```
pip install proposal==7.1.1
```

MUTE has been tested with PROPOSAL v7.1.1; earlier or later versions are not guaranteed to work. Environment set-up help for PROPOSAL installation on a Mac can be found [here](docs/Installing_PROPOSAL_on_a_Mac.md).

## Getting Started

For a basic example and a detailed description of how to use MUTE, see the [Tutorial](docs/Tutorial.md). For further examples, see the [examples](examples).

## Citation

Please cite https://inspirehep.net/literature/1927720.

The citations for the models and propagation tools used by MUTE can be found in the [MCEq](https://github.com/afedynitch/MCEq#please-cite-our-work) and [PROPOSAL](https://github.com/tudo-astroparticlephysics/PROPOSAL#how-to-cite-proposal) documentation.

### Authors

[William Woodley](mailto:wwoodley@ualberta.ca)

### Contributors

Anatoli Fedynitch and Marie-Cécile Piro

## License

MUTE is licensed under the BSD 3-Clause License (see [LICENSE](LICENSE)).
