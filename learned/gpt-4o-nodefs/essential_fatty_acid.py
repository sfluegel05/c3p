"""
Classifies: CHEBI:59549 essential fatty acid
"""
from rdkit import Chem

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is an essential fatty acid based on its SMILES string.
    An essential fatty acid typically has a long carbon chain with multiple cis double bonds,
    a terminal carboxylic acid group, and potentially functional groups that do not detract
    from the fatty acid definition.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as an essential fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for terminal carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Missing terminal carboxylic acid group (C(=O)O)"
    
    # Check for multiple cis double bonds (excluding potentially aromatic)
    cis_double_bond_pattern = Chem.MolFromSmarts("C/C=C/C")
    cis_double_bond_matches = mol.GetSubstructMatches(cis_double_bond_pattern)
    if len(cis_double_bond_matches) < 3:  # Adjusted minimum number of cis double bonds
        return False, f"Insufficient number of cis double bonds, found {len(cis_double_bond_matches)}"
    
    # Check the overall length of the carbon chain
    carbon_chain_pattern = Chem.MolFromSmarts("C" * 12)  # Assume minimum 12 continuous carbons
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "Insufficient carbon chain length consistent with essential fatty acids"

    # Include check for phosphocholine group
    phosphocholine_pattern = Chem.MolFromSmarts("[N+](C)(C)[C]")
    if mol.HasSubstructMatch(phosphocholine_pattern):
        return True, "Matches essential fatty acid criteria including phosphocholine functionality"

    return True, "Matches essential fatty acid criteria with multiple cis double bonds and a carboxylic acid group"