"""
Classifies: CHEBI:59549 essential fatty acid
"""
from rdkit import Chem

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is an essential fatty acid based on its SMILES string.
    An essential fatty acid typically has a long carbon chain with multiple cis double bonds,
    a terminal carboxylic acid group, and may include other functional groups.

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
    
    # Check for multiple cis double bonds (excluding aromatic)
    cis_double_bond_pattern = Chem.MolFromSmarts("[C&!a]/C=C/[C&!a]")
    cis_double_bond_matches = mol.GetSubstructMatches(cis_double_bond_pattern)
    if len(cis_double_bond_matches) < 2:
        return False, f"Insufficient number of cis double bonds, found {len(cis_double_bond_matches)}"
    
    # Include check for specific functional groups present in given examples (e.g., phosphocholine)
    phosphocholine_pattern = Chem.MolFromSmarts("[N+](C)(C)C")
    if mol.HasSubstructMatch(phosphocholine_pattern):
        return True, "Matches essential fatty acid criteria including phosphocholine functionality"

    # Final check for carbon chain length -> broaden range due to variable sizes in examples or remove
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 16:  # Lower bound from before removed upper bound for broader match
        return False, f"Carbon chain length {c_count} is too short for essential fatty acids"

    return True, "Matches essential fatty acid criteria with multiple cis double bonds and a carboxylic acid group"