"""
Classifies: CHEBI:20156 3-oxo-Delta(1) steroid
"""
from rdkit import Chem

def is_3_oxo_Delta_1_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(1) steroid based on its SMILES string.
    A 3-oxo-Delta(1) steroid is characterized by a double bond between positions 1 and 2
    and a ketone group at the 3-position on the steroid backbone.

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

    # SMARTS pattern for double bond between positions 1 and 2 in a steroid backbone
    delta_1_pattern = Chem.MolFromSmarts("C1=CCCCC1")

    # SMARTS pattern for ketone group at the 3-position
    ketone_3_pattern = Chem.MolFromSmarts("C(=O)[C;R]") 

    if not mol.HasSubstructMatch(delta_1_pattern):
        return False, "No double bond between positions 1 and 2"
        
    if not mol.HasSubstructMatch(ketone_3_pattern):
        return False, "No ketone group at the 3-position"

    return True, "Matches the 3-oxo-Delta(1) steroid pattern"