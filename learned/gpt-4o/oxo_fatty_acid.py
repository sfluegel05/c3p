"""
Classifies: CHEBI:59644 oxo fatty acid
"""
from rdkit import Chem

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.
    An oxo fatty acid contains at least one aldehydic or ketonic group in addition to 
    a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oxo fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group (-C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Look for aldehyde group (O=C-H)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    has_aldehyde = mol.HasSubstructMatch(aldehyde_pattern)

    # Look for ketone group (C(=O)C)
    ketone_pattern = Chem.MolFromSmarts("[#6][CX3](=O)[#6]")
    has_ketone = mol.HasSubstructMatch(ketone_pattern)

    # A fatty acid should have a long chain of carbon atoms
    chain_length_threshold = 6  # example threshold for chain length
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Determine if at least one of aldehyde or ketone group is present with the carboxylic group
    if (has_aldehyde or has_ketone) and carbon_count > chain_length_threshold:
        return True, "Contains carboxylic acid group, long carbon chain, and at least one oxo group (aldehydic or ketonic)"

    return False, "Missing appropriate aldehydic or ketonic group"