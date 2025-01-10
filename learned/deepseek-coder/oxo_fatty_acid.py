"""
Classifies: CHEBI:59644 oxo fatty acid
"""
"""
Classifies: CHEBI:77667 oxo fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.
    An oxo fatty acid is a fatty acid containing at least one aldehydic or ketonic group in addition to the carboxylic acid group.

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

    # Check for carboxylic acid group (-C(=O)OH or -C(=O)[O-])
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for at least one additional oxo group (aldehyde or ketone)
    oxo_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    if len(oxo_matches) < 2:  # At least one oxo group besides the carboxylic acid
        return False, f"Found {len(oxo_matches) - 1} oxo groups, need at least one"

    # Check for a long carbon chain (typical of fatty acids)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 8:  # Minimum number of carbons for a fatty acid
        return False, f"Only {c_count} carbons, need at least 8 for a fatty acid"

    # Check for a hydrophobic tail (long carbon chain without polar groups)
    hydrophobic_pattern = Chem.MolFromSmarts("[CX4][CX4][CX4][CX4]")
    hydrophobic_matches = mol.GetSubstructMatches(hydrophobic_pattern)
    if len(hydrophobic_matches) < 2:  # At least two consecutive hydrophobic segments
        return False, "Insufficient hydrophobic tail for a fatty acid"

    return True, "Contains a carboxylic acid group and at least one additional oxo group with a long carbon chain"