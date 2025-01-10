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

    # Various patterns associated with neoflavonoids (broader than previously)
    # Chroman and furan backbone-like structures could be seen in neoflavonoids
    chroman_like_patterns = [
        # Pattern for central chroman structure or similar
        Chem.MolFromSmarts('O1C=2C(C=CC(O)=C2C3=CC=CC=C3)=CO1'),
        # Furan-like structure in some neoflavonoids
        Chem.MolFromSmarts('O1C=2C(C=CC(O)=C2)CO1')
    ]

    # Checking for any of these potential backbones
    if not any(mol.HasSubstructMatch(pattern) for pattern in chroman_like_patterns):
        return False, "No chroman-like or furan-like pattern found typical of neoflavonoids"
    
    # Check for functional groups commonly seen in neoflavonoids
    # Hydroxyl, methoxy, ester, ketone, etc.
    # More diversified search for functional groups than previous
    functional_group_patterns = [
        # Hydroxyl groups, methoxy, ethers
        Chem.MolFromSmarts('[OX2H]'),      # Hydroxyl
        Chem.MolFromSmarts('CO'),          # Methoxy/Ether
        Chem.MolFromSmarts('[CX3](=O)O'),  # Ester
        Chem.MolFromSmarts('C=O')          # Carbonyl
    ]

    if not any(mol.HasSubstructMatch(fg_pattern) for fg_pattern in functional_group_patterns):
        return False, "Lacks characteristic oxygenated functional groups of neoflavonoids"

    return True, "Contains chroman-like or furan-like backbone with characteristic functional groups of neoflavonoids"