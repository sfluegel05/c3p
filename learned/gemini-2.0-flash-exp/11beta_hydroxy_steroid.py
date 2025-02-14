"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 11beta-hydroxy steroid based on its SMILES string.
    An 11beta-hydroxy steroid has a steroid core with a beta-configured hydroxyl group at position 11.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an 11beta-hydroxy steroid, False otherwise
        str: Reason for the classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the 11-beta hydroxyl group SMARTS pattern, using the specific stereo bond.
    # This pattern focuses on the carbon at the 11-position and its neighbours,
    # and not the whole steroid ring system, which is often hard to generalize.
    # Note that [C@H] indicates beta (up) config
    hydroxy_11beta_pattern = Chem.MolFromSmarts("[C]([C])[C@H](O)[C]")
    
    if not mol.HasSubstructMatch(hydroxy_11beta_pattern):
        return False, "No 11-beta-hydroxy group found"

    # Define a relaxed steroid core with just 4 fused rings with specific number of C atoms in each
    steroid_core_pattern = Chem.MolFromSmarts("[C]1[C]2[C]3[C]4[C]1[C]5[C]2[C]([C]3)[C]45")

    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "No steroid core structure detected"


    return True, "Molecule has a steroid core with an 11-beta-hydroxyl group"