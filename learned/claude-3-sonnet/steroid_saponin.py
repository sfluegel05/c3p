"""
Classifies: CHEBI:61655 steroid saponin
"""
"""
Classifies: CHEBI:50047 steroid saponin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdmolfiles

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.
    A steroid saponin is any saponin derived from a hydroxysteroid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid saponin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for cyclopentanoperhydrophenanthrene scaffold (core steroid structure)
    steroid_core_pattern = rdmolfiles.MolFromSmarts("[C&r5,r6]1[C&r5,r6][C&r5,r6][C&r5,r6][C&r5,r6][C&r5,r6]1[C&r6]2[C&r5,r6][C&r5,r6][C&r5,r6][C&r5,r6][C&r6]2")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "No steroid core found"
    
    # ... (rest of the code remains the same) ...