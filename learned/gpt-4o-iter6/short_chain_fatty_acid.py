"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
from rdkit import Chem

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid based on its SMILES string.
    A short-chain fatty acid is defined as an aliphatic monocarboxylic acid
    with a chain length of less than C6 and no non-hydrocarbon substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a short-chain fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES, return False if invalid
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Count the number of carbon atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Check if carbon chain length is less than 6
    if carbon_count >= 6:
        return False, f"Carbon chain too long, found {carbon_count} carbons"

    # Check for non-hydrocarbon substituents (atoms other than C, H, O)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [1, 6, 8]:  # H, C, O are allowed
            return False, "Contains non-hydrocarbon substituents"

    return True, "Valid short-chain fatty acid with appropriate carbon chain length"

__metadata__ = {
    'chemical_class': {
        'name': 'short-chain fatty acid',
        'definition': 'An aliphatic monocarboxylic acid with a chain length of less than C6. If any non-hydrocarbon substituent is present, the compound is not normally regarded as a short-chain fatty acid.'
    }
}