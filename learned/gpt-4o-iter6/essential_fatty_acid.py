"""
Classifies: CHEBI:59549 essential fatty acid
"""
from rdkit import Chem

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is an essential fatty acid based on its SMILES string.
    An essential fatty acid is a polyunsaturated fatty acid required in the diet.

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

    # Check for terminal carboxylic acid group (-C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No terminal carboxylic acid group found"
    
    # Identify cis double bonds (Z or \), which are critical for polyunsaturation
    cis_double_bond_pattern = Chem.MolFromSmarts("C/C=C/C")
    cis_bond_count = len(mol.GetSubstructMatches(cis_double_bond_pattern))
    if cis_bond_count < 2:
        return False, f"Insufficient cis double bonds (found {cis_bond_count}), required polyunsaturation"

    # Check for sufficiently long aliphatic carbon chain
    aliphatic_chain_pattern = Chem.MolFromSmarts("CCCCCCCC")
    if not mol.HasSubstructMatch(aliphatic_chain_pattern):
        return False, "Aliphatic chain length insufficient for essential fatty acid"
    
    return True, "Contains key characteristics of essential fatty acid: carboxylic acid group, multiple cis double bonds, and long aliphatic chain"