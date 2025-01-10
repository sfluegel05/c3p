"""
Classifies: CHEBI:71971 neoflavonoid
"""
from rdkit import Chem

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid based on its SMILES string.
    A neoflavonoid is typically characterized by a core benzopyran with a phenyl group at the 4-position and various functional group substitutions.

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

    # Define potential neoflavonoid core patterns
    benzopyran_coumarin_patterns = [
        # Neoflavonoid pattern: benzopyran core with phenyl at 4-position
        Chem.MolFromSmarts('O1C=CC2=CC=CC=C2C=1C3=CC=CC=C3'),  # Basic benzopyran with phenyl group
        Chem.MolFromSmarts('O=C1OC=CC2=CC=CC=C12'),  # Coumarin core
    ]

    # Checking for neoflavonoid core structures
    if not any(mol.HasSubstructMatch(pattern) for pattern in benzopyran_coumarin_patterns):
        return False, "No core benzopyran or coumarin backbone found typical of neoflavonoids"
    
    # Check for functional groups
    functional_group_patterns = [
        Chem.MolFromSmarts('[OX2H]'),      # Hydroxyl
        Chem.MolFromSmarts('CO'),          # Methoxy/Ether
        Chem.MolFromSmarts('[CX3](=O)O'),  # Ester
        Chem.MolFromSmarts('C=O')          # Carbonyl
    ]

    # Ensure molecule has several of the functional groups for neoflavonoids
    if not any(mol.HasSubstructMatch(fg_pattern) for fg_pattern in functional_group_patterns):
        return False, "Lacks characteristic oxygenated functional groups of neoflavonoids"

    return True, "Contains core benzopyran or coumarin backbone with characteristic functional groups of neoflavonoids"