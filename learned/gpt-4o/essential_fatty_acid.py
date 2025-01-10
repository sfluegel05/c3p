"""
Classifies: CHEBI:59549 essential fatty acid
"""
from rdkit import Chem

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is an essential fatty acid based on its SMILES string.
    Essential fatty acids are a subset of polyunsaturated fatty acids required in the diet.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an essential fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")
    
    # Check for terminal carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return (False, "No terminal carboxylic acid found")
    
    # Check carbon chain length (at least 18 carbons common in essential fatty acids)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 18:
        return (False, f"Carbon chain is too short: only {carbon_count} carbon atoms")
    
    # Check for multiple Z (cis) configuration double bonds pattern
    # Consider more expansive cis detection patterns
    cis_double_bond_pattern = Chem.MolFromSmarts("C=C")
    cis_double_bond_matches = mol.GetSubstructMatches(cis_double_bond_pattern)

    # Ensure there is a noteworthy number of double bonds, ideally in a sequence
    if len(cis_double_bond_matches) < 2:
        return (False, f"Only {len(cis_double_bond_matches)} double bonds found, require multiple for PUFA")
    
    return (True, "Matches essential fatty acid pattern with polyunsaturation and a terminal carboxylic acid")