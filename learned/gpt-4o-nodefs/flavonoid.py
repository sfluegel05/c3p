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

    # Recognize broader flavonoid motifs
    patterns = [
        # General flavone pattern with an open backbone
        Chem.MolFromSmarts('Oc1ccccc1-c2ccc(O)c(O)c2'), 
        
        # Isoflavone pattern, A-C ring structure
        Chem.MolFromSmarts('Oc1ccc2c(c1)ccc(=O)o2'), 
        
        # Flavanone patterns (saturated C ring)
        Chem.MolFromSmarts('C1CCC(=O)c2ccccc2O1'), 
        
        # Catechin/Anthocyanidin pattern, involving oxidanium
        Chem.MolFromSmarts('O[c+]1ccc(cc1)c2ccc(O)cc2'),

        # Broader 2-phenylchromen-4-one pattern
        Chem.MolFromSmarts('O=C1c2ccccc2oc3ccccc13'),

        # Glycosylated form, allow a generic attachment point
        Chem.MolFromSmarts('O1C2=CC(=CC=C2C(=O)O1)-C=3C=CC(O)=C(O)C=3'),
    ]

    for pattern in patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains flavonoid-like backbone or functional derivative"

    return False, "No flavonoid-like backbone or pattern found"