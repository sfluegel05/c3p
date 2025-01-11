"""
Classifies: CHEBI:38131 lactol
"""
from rdkit import Chem

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol based on its SMILES string.
    A lactol is characterized as a cyclic hemiacetal, where an -OH group
    and an ether oxygen are bonded to the same carbon within a ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lactol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Define the lactol pattern - cyclic hemiacetal [C;R]([OH])[O;R]
    lactol_pattern = Chem.MolFromSmarts('[C;R]([OH])[O;R]')
    if mol.HasSubstructMatch(lactol_pattern):
        return True, "Contains cyclic hemiacetal structure"
    else:
        return False, "Does not contain cyclic hemiacetal structure"

# Example usage:
# is_lactol("O1C2=C(C(C[C@@]1...=C2)O)O")