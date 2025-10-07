### IMPORTS AND STATIC STUFF ###
from bepipocket.normal_run import normal_run
from bepipocket.bepipocket_run import bepipocket_run
#from bepipocket.discopocket_run import discopocket_run
from pathlib import Path
import argparse 

### COMMAND LINE ARGUMENTS ###
parser = argparse.ArgumentParser("Create CHAI-1 models with different modelling mode")
parser.add_argument("-i", required=True, action="store", dest="fasta_file", type=Path, help="Input directory of fasta input files for Chai-1 to run.")
parser.add_argument("-ir", action="store", dest="restraint_file", default=None, type=Path, help="Input directory for restraint files for Chai-1 to run. Only required if running a restraint mode.")
parser.add_argument("-o", required=True, action="store", dest="out_dir", type=Path, help="Structure output directory")
#parser.add_argument("-pred", action="store", choices=["normal", "restraint", "bepipocket", "discopocket"], required=True, dest="pred", help="Chai-1 structure modelling mode to run.") #TODO
parser.add_argument("-pred", action="store", choices=["normal", "restraint", "bepipocket"], required=True, dest="pred", help="Chai-1 structure modelling mode to run.")
parser.add_argument("-nr_runs", action="store", dest="nr_runs", type=int, default=6, help="Number of runs for Chai-1 structure modelling.")
parser.add_argument("-bp3scores", action="store", dest="bp3scores", default=None, type=Path, help="Path to dict containing precomputed BepiPred-3.0 scores. {FGKAJ...:array([0.4,,0.3,0.5,0.6,0.8...])..}.")
parser.add_argument("-patch_mode", action="store_true", dest="patch_mode", help="Use surface patch with high epitope propensity for restraint runs for BepiPocket or DiscoPocket. Note: This approach has not been benchmarked yet).")
#parser.add_argument("-disco3scores", action="store", dest="disco3scores", default=None, type=Path, help="Path to dict containing precomputed discotope-3 scores")
parser.add_argument("-msa_directory", action="store", dest="msa_directory", default=None, type=Path, help="Look for MSA .pqt files with sequence hash filenames mathcing query sequences in this directory.")

# set variables
args = parser.parse_args()
fasta_file = args.fasta_file
out_dir = args.out_dir
pred = args.pred
restraint_file = args.restraint_file
nr_runs= args.nr_runs

bp3scores = args.bp3scores
#disco3scores = args.disco3scores #TODO
patch_mode = args.patch_mode
msa_directory = args.msa_directory

# chai-1 normal prediction mode 
if pred == "normal":
    normal_run(fasta_file, out_dir, overwrite_earlier_jobcontent=False, seeds=nr_runs, msa_directory=msa_directory)

# chai-1 restraint prediction mode (User defined restraints, as described in Chai-1 documentation)
elif pred == "restraint":
    normal_run(fasta_file, out_dir, restraint_file=restraint_file, overwrite_earlier_jobcontent=False, nr_runs=nr_runs, msa_directory=msa_directory)

# chai-1 bepipocket (use BepiPred-3.0 to guide antibody-epitope restraints)
elif pred == "bepipocket":
    bepipocket_run(fasta_file, out_dir, patch_mode=patch_mode, nr_runs=nr_runs, bp3_score_lookup=bp3scores, msa_directory=msa_directory)

# chai-1 discopocket (use DiscoTope-3.0 to guide antibody-epitope restraints)
#TODO 
# elif pred == "discopocket":
#     abag_disco_score_key = out_dir.name.split("_discotopemap")[0]
#     discopocket_chairun(fasta_file, out_dir, patch_mode=patch_mode, nr_runs=seeds,
#                          discotope3_score_lookup = disco3scores, msa_directory=msa_directory, abag_disco_score_key=abag_disco_score_key)
