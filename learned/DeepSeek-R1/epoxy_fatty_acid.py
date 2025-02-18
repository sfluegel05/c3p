"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
"""
Classifies: CHEBI_90255 epoxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.
    An epoxy fatty acid contains a carboxylic acid group and an epoxide ring.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an epoxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for carboxylic acid group (COOH)
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group"
    
    # Check for epoxide ring (O-C-O in a 3-membered ring)
    epoxide_pattern = Chem.MolFromSmarts("[O]1CCO1")
    if not mol.HasSubstructMatch(epoxide_pattern):
        return False, "No epoxide ring found"
    
    # Basic check for chain length: at least 12 carbons (typical for fatty acids)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 12:
        return False, "Carbon chain too short for fatty acid"
    
    return True, "Contains carboxylic acid and epoxide ring"