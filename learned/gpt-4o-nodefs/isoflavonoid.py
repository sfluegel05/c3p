"""
Classifies: CHEBI:50753 isoflavonoid
"""
from rdkit import Chem

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    Isoflavonoids contain a core flavonoid structure with various substitutions.
    
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

    # Isoflavonoids are characterized by a chromenone derivative structure
    # Isoflavonoid patterns based on literature and known structures
    isoflavonoid_patterns = [
        Chem.MolFromSmarts("O=C1C=CC2=CC=C(O)C=C2OC1"),  # Basic isoflavonoid chromenone structure
        Chem.MolFromSmarts("O=C1C=CC2=CC=CC=C2O1"),      # Allow more diversity in core structure
        Chem.MolFromSmarts("C1=CC=C(C2=O)C=CC2=CO1"),      # 3-phenylchromen-4-one basic form
    ]

    for pattern in isoflavonoid_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a core isoflavonoid structure characteristic of isoflavonoids"

    # Inform about the potential misclassification
    return False, "No core isoflavonoid structure found"