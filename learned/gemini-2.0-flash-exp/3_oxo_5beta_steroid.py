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

    # Define SMARTS patterns
    # 3-oxo group (ketone at position 3) and 5beta configuration. Explicit hydrogens are important
    combined_pattern = Chem.MolFromSmarts("[C]1[C](=O)[C][C@]2([H])[C]1[C]")

    # Check for the combined 3-oxo and 5beta group
    if not mol.HasSubstructMatch(combined_pattern):
        return False, "Molecule does not contain a 3-oxo group with a beta configuration at position 5"
    

    return True, "Molecule is a 3-oxo-5beta-steroid"