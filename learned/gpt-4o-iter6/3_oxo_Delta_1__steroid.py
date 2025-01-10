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

    # SMARTS pattern for steroid structure with Delta(1) double bond 
    steroid_pattern = Chem.MolFromSmarts("C1(=C)[C@@H]2CC[C@]3([C@H]4CCC(=O)C=C[C@]4(C)[C@@H]3CC[C@@H]12)") 
    
    # Check if the steroid structure with Delta(1) double bond and 3-oxo group is present
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Structure does not match a 3-oxo-Delta(1) steroid pattern"
    
    return True, "Matches the 3-oxo-Delta(1) steroid pattern"