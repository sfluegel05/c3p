"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
"""
Classifies: CHEBI:XXXXX short-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid based on its SMILES string.
    A short-chain fatty acid is an aliphatic monocarboxylic acid with a chain length less than C6.
    Oxygen-containing substituents (e.g., hydroxyl and keto groups) are allowed, but other heteroatoms are not permitted.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a short-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group (monocarboxylic acid)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) != 1:
        return False, f"Found {len(carboxylic_acid_matches)} carboxylic acid groups, need exactly 1"
    
    # Check for aromaticity (should be aliphatic)
    if mol.GetRingInfo().NumAromaticRings() > 0:
        return False, "Contains aromatic rings, not aliphatic"
    
    # Check for rings (should be acyclic)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains rings, should be acyclic"
    
    # Check that molecule contains only C, H, and O atoms
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in ('C', 'H', 'O'):
            return False, f"Contains heteroatom {atom.GetSymbol()}, not permitted"
    
    # Count total number of carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons >= 6:
        return False, f"Contains {num_carbons} carbons, chain length must be less than C6"
    
    # Allowed functional groups (oxygen-containing groups like hydroxyl and keto groups)
    # Since we already ensured only C, H, and O atoms are present, and there is only one carboxylic acid group
    # Additional oxygen-containing groups are acceptable

    return True, "Is an aliphatic monocarboxylic acid with fewer than 6 carbons and only C, H, O atoms"