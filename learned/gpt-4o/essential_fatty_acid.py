"""
Classifies: CHEBI:59549 essential fatty acid
"""
from rdkit import Chem

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is an essential fatty acid based on its SMILES string.
    Essential fatty acids are a sub-set of polyunsaturated fatty acids required in the diet.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an essential fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for terminal carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No terminal carboxylic acid found"
    
    # Check for multiple cis double bonds pattern
    # A simple pattern looking for cis config can be represented as: C\C=C\C 
    cis_double_bond_pattern = Chem.MolFromSmarts("C/C=C/C")
    cis_double_bond_matches = mol.GetSubstructMatches(cis_double_bond_pattern)
    
    # Ensure there are multiple such double bonds
    if len(cis_double_bond_matches) < 2:
        return False, f"Only {len(cis_double_bond_matches)} cis double bonds, need 2 or more for PUFA"
    
    # Count total number of double bonds and ensure it's greater than or equal to the number found in common EFAs
    total_double_bond_pattern = Chem.MolFromSmarts("C=C")
    total_double_bond_matches = mol.GetSubstructMatches(total_double_bond_pattern)
    if len(total_double_bond_matches) < 2:
        return False, f"Only {len(total_double_bond_matches)} total double bonds, not sufficient for polyunsaturation"
    
    return True, "Matches essential fatty acid pattern with polyunsaturation and terminal carboxylic acid"