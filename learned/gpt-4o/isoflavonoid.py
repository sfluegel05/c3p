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

    # Define isoflavonoid core pattern:
    # Flexible core benzopyran with a phenyl group attached at position 3
    isoflavonoid_patterns = [
        Chem.MolFromSmarts("c1cc(-c2coc3c(c2)c(cc(c3)O)=O)ccc1"),  # The initial pattern
        Chem.MolFromSmarts("c1coc2c(c1)c(ccc2-c3ccccc3)=O"),      # Alternate tautomeric structure
        Chem.MolFromSmarts("Oc1ccc2nc(=O)c(cc2c1)-c3ccccc3")       # Open form possibility
        # Additional patterns can be added to match different stereochemistry or tautomeric forms
    ]

    # Check if it matches any of the isoflavonoid patterns
    for pattern in isoflavonoid_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains 1-benzopyran with aryl substituent at position 3, matching isoflavonoid structure"
    
    return False, "Does not match isoflavonoid structure: 1-benzopyran with aryl at position 3"