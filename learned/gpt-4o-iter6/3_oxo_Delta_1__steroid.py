"""
Classifies: CHEBI:20156 3-oxo-Delta(1) steroid
"""
from rdkit import Chem

def is_3_oxo_Delta_1_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(1) steroid based on its SMILES string.
    A 3-oxo-Delta(1) steroid is characterized by a steroid nucleus with
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
    
    # SMARTS for 3-oxo steroid with Delta(1) bond
    # Tetracyclic system with a ketone group at position 3 and a double bond between positions 1 and 2.
    stereochemistry_core = Chem.MolFromSmarts(
        "C1=CCC2C(C1)=CC=C3C(=O)CC4C3(CCC4)"
    )
    
    # Validate presence of the steroid backbone, double bond between 1-2, and 3-ketone
    if not mol.HasSubstructMatch(stereochemistry_core):
        return False, "Does not have required tetracyclic nucleus, Delta(1) bond, or 3-ketone"

    return True, "Matches the 3-oxo-Delta(1) steroid pattern"