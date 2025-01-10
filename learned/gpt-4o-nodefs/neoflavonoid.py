"""
Classifies: CHEBI:71971 neoflavonoid
"""
from rdkit import Chem

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a neoflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a phenylcoumarin pattern (backbone for many neoflavonoids)
    phenylcoumarin_pattern = Chem.MolFromSmarts('O=C1OC2=C([C]=C(C=2)C3=CC=CC=C3)C1')
    if not mol.HasSubstructMatch(phenylcoumarin_pattern):
        return False, "No phenylcoumarin backbone found"

    # Check for the variety of functional groups prevalent in neoflavonoids
    if not (mol.HasSubstructMatch(Chem.MolFromSmarts('O')) or  # Presence of hydroxyl groups
            mol.HasSubstructMatch(Chem.MolFromSmarts('OC')) or # Methoxy groups
            mol.HasSubstructMatch(Chem.MolFromSmarts('COC'))): # More complex ethers
        return False, "Lacks characteristic functional groups of neoflavonoids"
    
    # Further checks might include confirmations for additional specific substitutions
    # like additional oxygen groups, branching alkenyl chains etc.

    return True, "Contains phenylcoumarin backbone with characteristic functional groups of neoflavonoids"