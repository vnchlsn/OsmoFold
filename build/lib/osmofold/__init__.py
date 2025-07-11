# osmofold/__init__.py
from .osmofold_local import get_tfe, sasa_to_rasa, protein_unfolded_dG, protein_folded_dG, protein_ddG_folding
from .parallel import batch_process_pdbs
from .osmofold_lite import read_fasta, protein_ddG_folding_lite, protein_folded_dG_lite, protein_unfolded_dG_lite
