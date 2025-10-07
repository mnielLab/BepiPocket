### IMPORTS ###

from pathlib import Path
import sys
import pdb
MODULE_DIR = str( Path( Path(__file__).parent.resolve() ) )
sys.path.append(MODULE_DIR)
from chai_lab.chai1 import run_inference

### FUNCTIONS ###

def normal_run(fasta_path, outdir, restraint_file=None, seeds=6, num_trunk_recycles=4, num_diffn_timesteps=200, overwrite_earlier_jobcontent=False, msa_directory=None):

    """
    num_trunk_recycles=4, num_diffn_timesteps=200 are the default settings from Chai-1
    """

    for seed in range(seeds):
        out_path = outdir / f"seed{seed}"

        # delete results from earlier run (if wanting to re-run evaluation)
        if out_path.is_dir() and overwrite_earlier_jobcontent:
            for f in out_path.glob("*"): f.unlink()

        # check how many files are in seed (should be 10: 5 .cif (structure) + 5 .npz (confidence))
        if out_path.is_dir():
            nr_score_files = len( list(out_path.glob("*.npz")) )
            nr_structure_files = len( list(out_path.glob("*.cif")) )
            if nr_structure_files == 5 and nr_score_files == 5: 
                seed_file_check = True

        else: seed_file_check = False

        if seed_file_check:
            print(f"Skipping. Found all structure and conf. files for {seed} at {str(out_path)}.")
            continue

        # create restraint file
        if restraint_file is None:
            if msa_directory is not None: msa_directory = Path(msa_directory)
            run_inference(fasta_file=fasta_path, output_dir=out_path,
                          num_trunk_recycles=num_trunk_recycles,
                          num_diffn_timesteps=num_diffn_timesteps,
                          seed=seed, use_esm_embeddings=True, msa_directory=msa_directory)
            
        else:
            print(f"Doing CHAI job with restraints: {restraint_file}")
            if msa_directory is not None: msa_directory = Path(msa_directory)
            run_inference(fasta_file=fasta_path, output_dir=out_path,
                          constraint_path=restraint_file,
                          num_trunk_recycles=num_trunk_recycles,
                          num_diffn_timesteps=num_diffn_timesteps,
                          seed=seed, use_esm_embeddings=True, msa_directory=msa_directory)
        
    # write done file 
    outfile = open(outdir / "done.txt", "w")
    outfile.close()