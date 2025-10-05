### IMPORTS AND STATIC STUFF ###
from chaiabag.normal_chairun import normal_chairun
from chaiabag.surfmap_chairun import surfmap_chairun
from chaiabag.bepipredmap_chairun import bepipredmap_chairun
from chaiabag.discotopemap_chairun import discotopemap_chairun

from pathlib import Path
import argparse 

### COMMAND LINE ARGUMENTS ###
parser = argparse.ArgumentParser("Create CHAI-1 models with different modelling mode")
parser.add_argument("-i", required=True, action="store", dest="fasta_file", type=Path, help="Input directory of fasta input files for Chai-1 to run.")
parser.add_argument("-ir", action="store", dest="restraint_file", default=None, type=Path, help="Input directory for restraint files for Chai-1 to run. Only required if running a restraint mode.")
parser.add_argument("-o", required=True, action="store", dest="out_dir", type=Path, help="Structure output directory")
parser.add_argument("-pred", action="store", choices=["normal", "restraint", "surfmap", "bepipredmap", "discotopemap"], required=True, dest="pred", help="Chai-1 structure modelling mode to run.")
parser.add_argument("-seeds", action="store", dest="seeds", type=int, default=6, help="Number of seeds to run Chai-1 structure modelling mode to run.")
parser.add_argument("-surfmap_maxruns", action="store", dest="surfmap_maxruns", type=int, default=50, help="Maximum number of Chai-1 surfmap runs to run.")
parser.add_argument("-bp3scores", action="store", dest="bp3scores", default=None, type=Path, help="Path to dict containing precomputed BepiPred-3.0 scores. {FGKAJ...:array([0.4,,0.3,0.5,0.6,0.8...])..}")
parser.add_argument("-disco3scores", action="store", dest="disco3scores", default=None, type=Path, help="Path to dict containing discotope-3 scores")
parser.add_argument("-patch_mode", action="store_true", dest="patch_mode", help="Use surface patches for when running surfmap or bepipredmap mode")
parser.add_argument("-maxcover_mode", action="store_true", dest="maxcover_mode", help="Skip BepiPred-3.0 predicted epitope residues that already been covered.")
parser.add_argument("-msa_directory", action="store", dest="msa_directory", default=None, type=Path, help="Look for MSA .pqt files with sequence hash filenames mathcing query sequences in this directory.")

# set variables
args = parser.parse_args()
fasta_file = args.fasta_file
out_dir = args.out_dir
pred = args.pred
restraint_file = args.restraint_file
seeds= args.seeds
surfmap_maxruns = args.surfmap_maxruns
bp3scores = args.bp3scores
disco3scores = args.disco3scores
patch_mode = args.patch_mode
maxcover_mode = args.maxcover_mode
msa_directory = args.msa_directory

# chai-1 normal prediction mode 
if pred == "normal":
    normal_chairun(fasta_file, out_dir, overwrite_earlier_jobcontent=False, seeds=seeds, msa_directory=msa_directory)
# chai-1 restraint prediction mode
elif pred == "restraint":
    normal_chairun(fasta_file, out_dir, restraint_file=restraint_file, overwrite_earlier_jobcontent=False, seeds=seeds, msa_directory=msa_directory)
#chai-1 surfmap prediciton mode (using DSSP comp. surf. acces. to cover antigen surface)
elif pred == "surfmap":
    surfmap_chairun(fasta_file, out_dir, overwrite_earlier_jobcontent=True, surfmap_maxruns=surfmap_maxruns, patch_mode=patch_mode)
# chai1 bepipred3map mode (using BepiPred-3.0 to cover antigen surface)
elif pred == "bepipredmap":
    bepipredmap_chairun(fasta_file, out_dir, patch_mode=patch_mode, bepipredmap_runs=seeds, bp3_score_lookup=bp3scores, maxcover_mode=maxcover_mode, msa_directory=msa_directory)
elif pred == "discotopemap":
    abag_disco_score_key = out_dir.name.split("_discotopemap")[0]
    discotopemap_chairun(fasta_file, out_dir, patch_mode=patch_mode, discotopemap_runs=seeds,
                         discotope3_score_lookup = disco3scores, maxcover_mode=maxcover_mode, msa_directory=msa_directory, abag_disco_score_key=abag_disco_score_key)