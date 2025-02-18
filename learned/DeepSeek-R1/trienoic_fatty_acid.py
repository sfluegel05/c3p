"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
from rdkit import Chem

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    A trienoic fatty acid is a polyunsaturated fatty acid with exactly three double bonds in the carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trienoic fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for exactly one carboxylic acid group
    carboxylic_acid = Chem.MolFromSmarts("[CX3](=O)[OH]")
    matches = mol.GetSubstructMatches(carboxylic_acid)
    if len(matches) != 1:
        return False, f"Found {len(matches)} carboxylic acid groups, need exactly 1"

    # Count total double bonds and subtract the COOH group's bond
    total_double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    chain_double_bonds = total_double_bonds - len(matches)
    if chain_double_bonds != 3:
        return False, f"Found {chain_double_bonds} double bonds in carbon chain, need exactly 3"

    # Check for minimum carbon count (at least 12)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 12:
        return False, f"Only {carbon_count} carbons, insufficient for a fatty acid"

    return True, "Carboxylic acid with three double bonds in carbon chain"