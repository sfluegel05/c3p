"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
"""
Classifies: cyclic fatty acid (CHEBI) - Any fatty acid containing a ring structure
"""
from rdkit import Chem

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid based on its SMILES string.
    A cyclic fatty acid must contain:
    1. At least one carboxylic acid group (-COOH)
    2. At least one ring structure (any size)
    3. The carboxylic acid must be part of a carbon chain structure

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclic fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for any rings in the structure
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 0:
        return False, "No rings found in structure"

    # Basic check for carbon chain structure (at least 4 carbons in a row)
    carbon_chain_pattern = Chem.MolFromSmarts("C-C-C-C")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "Insufficient carbon chain structure"

    return True, "Contains carboxylic acid group and cyclic structure"