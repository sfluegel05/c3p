"""
Classifies: CHEBI:38131 lactol
"""
from rdkit import Chem

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol based on its SMILES string.
    A lactol is characterized as a cyclic hemiacetal, which involves a -OH group
    and an ether oxygen bonded to the same carbon within a ring.

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

    # Define the lactol pattern - cyclic hemiacetal pattern with better specificity
    # Pattern: Carbon in a ring bonded to an hydroxyl group and an ether oxygen
    lactol_pattern = Chem.MolFromSmarts('[C&R]([O&R][C&R])([OH&R])') 
    
    if mol.HasSubstructMatch(lactol_pattern):
        return True, "Contains cyclic hemiacetal (lactol) structure"
    else:
        return False, "Does not contain cyclic hemiacetal (lactol) structure"

# Example usage:
# result, reason = is_lactol("O1C2=C(C(C[C@@]1(C3=CC=...=C2)O)O")
# print(result, reason)