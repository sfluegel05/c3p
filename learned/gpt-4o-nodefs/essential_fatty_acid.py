"""
Classifies: CHEBI:59549 essential fatty acid
"""
from rdkit import Chem

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is an essential fatty acid based on its SMILES string.
    An essential fatty acid typically has a long carbon chain with multiple cis double bonds
    and a terminal carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as an essential fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylic acid group at the end of the chain
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Missing terminal carboxylic acid group (C(=O)O)"
    
    # Check for multiple cis double bonds
    cis_double_bond_pattern = Chem.MolFromSmarts("C/C=C/C")
    cis_double_bond_matches = mol.GetSubstructMatches(cis_double_bond_pattern)
    if len(cis_double_bond_matches) < 2:
        return False, f"Insufficient number of cis double bonds, found {len(cis_double_bond_matches)}"
    
    # Count carbon atoms to ensure long chain length
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18 or c_count > 36:
        return False, f"Carbon chain length {c_count} is outside typical range for essential fatty acids (18 to 36 carbons)"
    
    return True, "Matches essential fatty acid criteria with multiple cis double bonds and a carboxylic acid group"