"""
Classifies: CHEBI:50753 isoflavonoid
"""
from rdkit import Chem

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    An isoflavonoid is defined as any 1-benzopyran with an aryl substituent at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a simplified core pattern for isoflavonoids
    # 1-Benzopyran with an aryl (phenyl or similar aromatic ring) at position 3
    isoflavonoid_patterns = [
        Chem.MolFromSmarts("c1cc2occc2c3ccccc3c1"),  # Basic 1-benzopyran system with aryl
        Chem.MolFromSmarts("c1cc2occc2c(c1)-c3ccccc3"),  # Basic form with alternative aryl linker
        Chem.MolFromSmarts("c1cc2occc2c(c1C3)C3"), # Open form with different tautomers
        # Consider more specific patterns if needed
    ]

    # Check if the molecule matches any of the isoflavonoid patterns
    for pattern in isoflavonoid_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains 1-benzopyran with aryl substituent at position 3, matching isoflavonoid structure"
    
    return False, "Does not match isoflavonoid structure: 1-benzopyran with aryl at position 3"