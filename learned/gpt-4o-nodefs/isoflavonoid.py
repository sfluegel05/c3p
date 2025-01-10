"""
Classifies: CHEBI:50753 isoflavonoid
"""
from rdkit import Chem

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    Isoflavonoids are characterized by a chromen-4-one structure with a phenyl group at the C2 position.
    
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

    # Updated pattern for isoflavonoid core structure with variability
    isoflavonoid_patterns = [
        Chem.MolFromSmarts("O=C1C=CC2=C(O)C=CC=C2O1"),  # More generalized chromenone structure
        Chem.MolFromSmarts("C1=CC(=O)C2=C(C1)C=CC(O)=C2"), # Basic 3-phenylchromen-4-one form
        Chem.MolFromSmarts("c1coc2c1c(=O)cc(O)c2"), # Common isoflavonoid features
        Chem.MolFromSmarts("Oc1ccc2C=CC(=O)C=C2c1"),   # Variations with phenolic hydroxyl groups
    ]

    # Allow for variations by checking several patterns
    for pattern in isoflavonoid_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a core isoflavonoid structure characteristic of isoflavonoids"

    # Inform about the potential misclassification
    return False, "No core isoflavonoid structure found"