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

    # Define core patterns for isoflavonoids, accommodating both possible orientations and substitutions
    isoflavonoid_patterns = [
        Chem.MolFromSmarts("c1cc-2coc3c(c2c(c1)O)cc(=O)c4c(cccc4)-c3"),  # 1-benzopyran with aryl
        Chem.MolFromSmarts("c1cc-2c(c3c(cccc3)cc(=O)c4c2oc(c1C4)*)Oc"),  # Expanded for variety in aromatics
        Chem.MolFromSmarts("c1c-2c(c3c(coc3)cc(=O)c4c2oc(c1)C4)-c5ccccc5") # Open form and tautomers
        # Additional robust patterns to cover stereochemistry
    ]

    # Check if the molecule matches any of the isoflavonoid patterns
    for pattern in isoflavonoid_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains 1-benzopyran with aryl substituent at position 3, matching isoflavonoid structure"
    
    return False, "Does not match isoflavonoid structure: 1-benzopyran with aryl at position 3"