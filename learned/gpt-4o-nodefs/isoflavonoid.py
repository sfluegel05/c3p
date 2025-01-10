"""
Classifies: CHEBI:50753 isoflavonoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    Isoflavonoids are characterized by a 1-benzopyran-4-one structure.
    
    NOTE: This is a heuristic approach and may not capture all isoflavonoids accurately.
    
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

    # Look for 1-benzopyran-4-one pattern
    isoflavonoid_pattern = Chem.MolFromSmarts("O=C1C=CC2=CC=CC=C2OC1")
    if not mol.HasSubstructMatch(isoflavonoid_pattern):
        return False, "No 1-benzopyran-4-one core structure found"

    # Additional checks could be performed based on other isoflavonoid features, such as substituents
    # These can include methoxy/hydroxy groups or various phenolic substructures

    return True, "Contains the core 1-benzopyran-4-one structure characteristic of isoflavonoids"