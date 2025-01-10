"""
Classifies: CHEBI:47916 flavonoid
"""
from rdkit import Chem

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Recognize several common flavonoid patterns
    # Basic flavonoid-like backbones and common substructures

    patterns = [
        # Flavone and flavonol basic pattern
        Chem.MolFromSmarts('c1cc(O)ccc1-c2c(=O)cc(-c3cc(O)ccc3)oc2'),  
    
        # Isoflavone pattern
        Chem.MolFromSmarts('c1cc(O)ccc1-c2cc(-c3cc(O)ccc3)oc(=O)c2'),
        
        # Flavanol pattern
        Chem.MolFromSmarts('c1cc(O)ccc1-c2c(O)cc(-c3cc(O)ccc3)oc2'),

        # Anthocyanidin pattern
        Chem.MolFromSmarts('[o+]=c1c(cc(O)cc1-c2cc(O)ccc2)cc(O)c3cc(O)ccc3'),
        
        # Generalized tricyclic system with hydroxy groups
        Chem.MolFromSmarts('c1cc(O)c(O)cc1-c2cc(O)cc(-c3cc(O)c(O)cc3)o2'),
        
        # Catechin pattern
        Chem.MolFromSmarts('c1cc(O)ccc1-c2cc(O)cc(-c3cc(O)ccc3)o2')
    ]

    for pattern in patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains flavonoid-like backbone"

    return False, "No flavonoid-like backbone found"