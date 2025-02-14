"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
"""
Classifies: CHEBI:38165 3-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.
    A 3-hydroxy fatty acid is any fatty acid with a hydroxy functional group in the beta- or 3-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_pattern)
    if len(carboxylic_matches) < 1:
        return False, "No carboxylic acid group found"
    
    # Check for long carbon chain (fatty acid) with possible cyclopropane rings
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~*~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 1:
        return False, "No long carbon chain found"
    
    # Look for hydroxyl group(s) along the carbon chain
    hydroxy_pattern = Chem.MolFromSmarts("[CX4]([CX3]([OX2H1])[CX4,CX3])")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if len(hydroxy_matches) < 1:
        return False, "No hydroxyl group found along the carbon chain"
    
    # Allow for methoxy groups
    methoxy_pattern = Chem.MolFromSmarts("[CX4]([CX3]([OX2C]))")
    methoxy_matches = mol.GetSubstructMatches(methoxy_pattern)
    
    # Check for branched chains
    branched_pattern = Chem.MolFromSmarts("[CX4]([CX3])([CX3])([CX3,CX4])")
    branched_matches = mol.GetSubstructMatches(branched_pattern)
    
    return True, "Contains a carboxylic acid group, a long carbon chain, and one or more hydroxyl groups along the chain"