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
    # Note that [C@H] indicates beta (up) config
    hydroxy_11beta_pattern = Chem.MolFromSmarts("[C]([C])[C@H](O)[C]")
    
    if not mol.HasSubstructMatch(hydroxy_11beta_pattern):
        return False, "No 11-beta-hydroxy group found"
    
    # Count the number of carbon atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 17:
         return False, f"Too few carbons, expected 17 in core, got {carbon_count}"
    
    # Check for four fused rings
    four_ring_pattern = Chem.MolFromSmarts("[C]1~[C]~[C]~[C]~1~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]")
    if not mol.HasSubstructMatch(four_ring_pattern):
        return False, "No steroid core structure detected (no 4 fused rings)"

    return True, "Molecule has a steroid core with an 11-beta-hydroxyl group"