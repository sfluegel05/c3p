"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5beta-steroid based on its SMILES string.
    
    A 3-oxo-5beta-steroid has:
        - a steroid core
        - a ketone at the 3 position
        - a beta configuration at position 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-5beta-steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for the steroid core with 3-oxo and 5beta
    # This pattern identifies the four fused rings, the 3-oxo group, and the beta hydrogen at position 5.
    # The pattern has the 5 ring with atom numbers 1, 2, 3, 4, 10
    # The 6 ring directly attached to ring 5 is identified with atoms numbered 5, 6, 7, 8, 9 and 10
    # Atoms 1 and 5 have a stereochemistry associated with them, with the pattern specifying that the hydrogen at 5 is above the plane of the ring, and thus a beta configuration.
    combined_pattern = Chem.MolFromSmarts("[C]1[C](=O)[C][C@]2([H])[C]1[C]3[C]([H])([H])[C]([H])([H])[C]4([H])[C@]([H])([C]([H])([H])[C]3[H])[C]2[C]4[H]")

    # Check for the combined 3-oxo and 5beta group
    if not mol.HasSubstructMatch(combined_pattern):
        return False, "Molecule does not contain a 3-oxo group with a beta configuration at position 5"
    

    return True, "Molecule is a 3-oxo-5beta-steroid"