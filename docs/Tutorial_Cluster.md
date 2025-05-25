# **Tutorial - Using MUTE on a Computer Cluster**

For high-statistics simulations on a computer cluster, the ``mtp.propagation_muons()`` function is the most useful. The MUTE code to set up a job array of Monte Carlo simulations for 1000 muons per energy-slant depth bin (per job) can be written as follows:

```python
# Import packages

import argparse

import mute.constants as mtc
import mute.propagation as mtp

# Parse the job array number from the command line

parser = argparse.ArgumentParser()
parser.add_argument("job_array_number", type = int)
args = parser.parse_args()

# Set the constants

mtc.set_verbose(0)
mtc.set_output(True)
mtc.set_directory("mute/data")
mtc.set_lab("Example")
mtc.set_n_muon(1000)

# Propagate the muons

mtp.propagate_muons(seed = args.job_array_number, job_array_number = args.job_array_number, force = True)
```

To prevent any unnecessary output, the verbosity can be set to ``0`` (though it might be a good idea to keep it at the default ``2`` to help figure out what went wrong if the job fails). To make it easier to access and ``sftp`` out the output Monte Carlo underground energy files, the directory can be changed to the working directory of the code or elsewhere. The ``force`` parameter in ``mtc.propagate_muons()`` should be set to ``True`` to force the creation of any required directories or files so the program does not hang, waiting for user input until the job times out.

The seed, which can be changed in the propagation function to ensure the Monte Carlo results are different for each job, and the job array number can be both read in from the command line using ``argparse``. If the job array consisted of 100 jobs, the survival probabilities can then be calculated as below, setting the ``n_job`` argument in the ``mtp.calc_survival_probability_tensor()`` function to ``100``:

```python
import mute.constants as mtc
import mute.propagation as mtp

mtc.set_lab("Example")
mtc.set_n_muon(100000)

mtp.calc_survival_probability_tensor(n_job = 100)
```

Here, ``seed`` does not need to be set in ``mtp.calc_survival_probability_tensor()`` because this function will not invoke ``mtp.propagate_muons()``, since the underground energies have already been loaded. Note also that the number of muons was set to 1000 when running the propagation, but is set to 100000 in the code just above. By setting ``n_job`` to 100, MUTE will recognise that the 100000 muons were split evenly between 100 jobs of 1000 muons each, and will search for underground energy files corresponding to 1000 muons.