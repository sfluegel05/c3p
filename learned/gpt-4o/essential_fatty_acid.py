"""
Classifies: CHEBI:59549 essential fatty acid
"""
from rdkit import Chem

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is an essential fatty acid based on its SMILES string.
    Essential fatty acids are long-chain polyunsaturated fatty acids with multiple
    cis double bonds and a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an essential fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Look for carboxylic acid group (COOH) at one end
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"

    # Check for multiple double bonds (C=C) along the chain
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bonds = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bonds) < 2:
        return False, "Too few double bonds; polyunsaturation required"

    # Verify the presence of cis configuration double bonds
    # Look for patterns of the form C/C=C/C or C\C=C\C, indicating cis double bonds
    cis_double_bond_pattern = Chem.MolFromSmarts("C/C=C/C")
    cis_double_bonds = mol.GetSubstructMatches(cis_double_bond_pattern)
    if len(cis_double_bonds) < len(double_bonds):
        # If insufficient cis double bonds detected, verify using alternative SMARTS
        alt_cis_double_bond_pattern = Chem.MolFromSmarts(r"C\C=C\C") 
        alt_cis_double_bonds = mol.GetSubstructMatches(alt_cis_double_bond_pattern)
        if len(cis_double_bonds) + len(alt_cis_double_bonds) < len(double_bonds):
            return False, "Not enough cis double bonds for essential fatty acid"

    # Verify the total number of carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 18:
        return False, "Too few carbons for essential fatty acid"

    return True, "Contains carboxylic acid and sufficient cis double bonds for essential fatty acid"