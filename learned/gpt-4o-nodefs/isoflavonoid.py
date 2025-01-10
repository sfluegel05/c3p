"""
Classifies: CHEBI:50753 isoflavonoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    Isoflavonoids are characterized by a 3-phenylchromen-4-one structure.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Look for 1-benzopyran-4-one pattern
    # The SMARTS pattern considers variability in isoflavonoid derivatives.
    isoflavonoid_patterns = [
        Chem.MolFromSmarts("O=C1C=CC2=CC=CC=C2OC1"),  # Basic 1-benzopyran-4-one
        Chem.MolFromSmarts("C1=CC=C2C=C(C=O)C=CO2C1"),  # 3-phenylchromen-4-one
    ]

    for pattern in isoflavonoid_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a core 1-benzopyran-4-one or 3-phenylchromen-4-one structure characteristic of isoflavonoids"

    # Note for potential misidentification due to missing patterns
    return False, "No core 1-benzopyran-4-one or 3-phenylchromen-4-one structure found"