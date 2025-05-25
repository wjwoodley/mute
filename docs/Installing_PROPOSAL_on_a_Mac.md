# Installing PROPOSAL on a Mac

Tested with Mac OS Version 15.4.1.

## Requirements

The following are required to set up a proper environment to install PROPOSAL:

* ``python3``
* ``brew``
* ``cmake``
* ``g++``
* ``xcode``

Update Homebrew and install ``gcc`` (to compile PROPOSAL) and ``xcode`` (for command line tools):

```
brew update
brew upgrade
brew install gcc
xcode-select --install
```

Unlock Conan:

```
conan remove --locks
```

## Set the Python Environment up

Use Python version 3.6 or higher. For example:

```
brew install pyenv
pyenv install 3.8.10
pyenv global 3.8.10
eval "$(pyenv init -)"
python -V
``` 

## Install PROPOSAL

Continue with the installation of PROPOSAL with ``pip``. MUTE has been tested with PROPOSAL v7.6.2; earlier or later versions are not guaranteed to work.

```
python3 -m pip install proposal==7.6.2
```

Verify that PROPOSAL imports properly in Python:

```
python3

>>> import proposal
>>>
```

If this works, PROPOSAL has been installed.