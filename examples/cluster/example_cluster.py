# This is an example of a simple script that can be sent to a computer cluster, such as Compute Canada, to run a high number of muons. It can be submitted to slurm with sbatch example_submission.sh.
# With the number of job arrays set to 100 in example_submission.sh and the number of muons set to 1000 in this file, this will give a total of 1e5 muons.
# This will output 100 files of underground energies which can be read into MUTE with the function mtp.load_u_energies_from_files() with the parameter n_job set to 100.
# Set the force parameter in mtp.propagate_muons() to True to force the creation of any required directories or files (so the program does not wait for your input until the job times out).
# Make sure proposal is installed on the cluster.

# Import packages

import mute.constants as mtc
import mute.propagation as mtp

import argparse

# Parse the job array number from the command line

parser = argparse.ArgumentParser()
parser.add_argument("a", type=int)
args = parser.parse_args()

# Set the constants

mtc.set_verbose(0)
mtc.set_output(True)
mtc.set_directory("mute/data")
mtc.set_lab("Example")
mtc.set_n_muon(1000)

# Propagate the muons

mtp.propagate_muons(seed=args.a, job_array_number=args.a, force=True)