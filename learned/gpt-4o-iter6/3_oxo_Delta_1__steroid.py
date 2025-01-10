"""
Classifies: CHEBI:20156 3-oxo-Delta(1) steroid
"""
from rdkit import Chem

def is_3_oxo_Delta_1_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(1) steroid based on its SMILES string.
    A 3-oxo-Delta(1) steroid is characterized by a steroid structure with
    a double bond between positions 1 and 2 and a ketone group at the 3-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-Delta(1) steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Simplified SMARTS for 3-oxo steroid with Delta(1) bond
    steroid_core_pattern = Chem.MolFromSmarts("C1=CC[C@@H]2[C@H](C=O)CC[C@]2(C1)")
    
    # Check if the molecule contains the substructure
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Does not match the 3-oxo-Delta(1) steroid pattern"

    return True, "Matches the 3-oxo-Delta(1) steroid pattern"