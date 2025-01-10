"""
Classifies: CHEBI:38131 lactol
"""
from rdkit import Chem

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol based on its SMILES string.
    A lactol is characterized as a cyclic hemiacetal, featuring a carbon atom
    bonded to both a hydroxyl group and an ether oxygen within a ring.

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

    # New lactol pattern
    # This pattern attempts to detect a carbon in a ring system bonded to an -OH group and having adjacent ether like cyclic linkage.
    lactol_smart_pattern = '[C;R]([OH])[O;R]'

    lactol_pattern = Chem.MolFromSmarts(lactol_smart_pattern)
    
    if lactol_pattern and mol.HasSubstructMatch(lactol_pattern):
        return True, "Contains cyclic hemiacetal (lactol) structure"
    else:
        return False, "Does not contain cyclic hemiacetal (lactol) structure"

# Example usage:
# result, reason = is_lactol("O1C2=C(C(C[C@@]1(C3=CC=CC=C3)O)=O)C(=CC(=C2)O)O")
# print(result, reason)