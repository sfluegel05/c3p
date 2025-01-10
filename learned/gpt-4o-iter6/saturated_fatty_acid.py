"""
Classifies: CHEBI:26607 saturated fatty acid
"""
from rdkit import Chem

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.
    A saturated fatty acid should have: 
    - A terminal carboxyl group
    - No carbon-carbon double or triple bonds

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saturated fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group (COOH)
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"

    # Check for C=C and C#C bonds (unsaturation)
    unsaturation_pattern = Chem.MolFromSmarts("C=C | C#C")
    if mol.HasSubstructMatch(unsaturation_pattern):
        return False, "Contains unsaturation (carbon-carbon double or triple bond)"

    # Check for long carbon chain, minimum 4 carbons with connected pattern
    carbon_chain_pattern = Chem.MolFromSmarts("[CH2]-[CH2]-[CH2]-[CH2]")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "Does not contain sufficiently long carbon chain"

    return True, "Contains a long saturated carbon chain with a terminal carboxylic acid group"