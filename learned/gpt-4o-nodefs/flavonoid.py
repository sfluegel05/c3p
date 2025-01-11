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

    # Look for structure characteristic of flavonoids: C6-C3-C6 system
    # Typically, 2 phenolic rings connected by a 3-carbon bridge 
    # forming a 6-membered heterocyclic ring flavonoid backbone.
    flavonoid_pattern = Chem.MolFromSmarts('c1cc(O)ccc1-c2ccccc2') # Basic flavonoid-like pattern
    if not mol.HasSubstructMatch(flavonoid_pattern):
        return False, "No flavonoid-like backbone found"
    
    # Check for additional hydroxy or methoxy groups which are common
    hydroxy_pattern = Chem.MolFromSmarts('[OX2H]')  # Hydroxy group
    methoxy_pattern = Chem.MolFromSmarts('OC')  # Methoxy group
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "Missing hydroxy groups found in flavonoids"
    if not mol.HasSubstructMatch(methoxy_pattern):
        return False, "Missing optional methoxy groups found in flavonoids"

    return True, "Flavonoid-like structure with phenolic rings and 3-carbon bridge"