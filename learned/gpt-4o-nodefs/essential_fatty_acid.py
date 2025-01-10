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
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1,-]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Missing terminal carboxylic acid group (C(=O)OH)"
    
    # Check for multiple cis configuration double bonds
    cis_double_bond_pattern = Chem.MolFromSmarts("C/C=C/C")
    cis_double_bond_matches = mol.GetSubstructMatches(cis_double_bond_pattern)
    if len(cis_double_bond_matches) < 2:  # Adjusted minimum number of cis double bonds
        return False, f"Insufficient number of cis double bonds, found {len(cis_double_bond_matches)}"
    
    # Check the carbon chain length
    carbon_atom_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_atom_count < 16:  # Essential fatty acids have at least 16 carbons
        return False, "Insufficient carbon chain length for essential fatty acids"

    # Include recognition of phosphocholine group in more complex structures
    phosphocholine_pattern = Chem.MolFromSmarts("P(=O)([O-])[O-]C[N+](C)(C)C")
    if mol.HasSubstructMatch(phosphocholine_pattern):
        return True, "Matches essential fatty acid criteria including phosphocholine functionality"

    return True, "Matches essential fatty acid criteria with multiple cis double bonds and a carboxylic acid group"