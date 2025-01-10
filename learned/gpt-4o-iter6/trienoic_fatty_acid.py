"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
from rdkit import Chem

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    A trienoic fatty acid is a polyunsaturated fatty acid with three double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trienoic fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Count double bonds (C=C)
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) != 3:
        return False, f"Found {len(double_bond_matches)} carbon-carbon double bonds, need exactly 3"

    # Verify that the double bonds are part of a long carbon chain
    carbon_chain_pattern = Chem.MolFromSmarts("[C](=C)~[C](=C)~[C](=C)")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "Double bonds are not part of a long carbon chain; may not be a typical trienoic acid"

    return True, "Contains carboxylic acid group and three double bonds in a hydrocarbon chain"