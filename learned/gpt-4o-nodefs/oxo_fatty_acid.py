"""
Classifies: CHEBI:59644 oxo fatty acid
"""
from rdkit import Chem

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oxo fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group (-C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Look for additional ketone group (not part of carboxylic acid)
    ketone_pattern = Chem.MolFromSmarts("C(=O)[C,c]")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "No additional ketone group found"

    # Check for sufficiently long carbon chain (optional - can adjust length requirement)
    carbon_chain_pattern = Chem.MolFromSmarts("[CX4]-[CX4]-[CX4]-[CX4]-[CX4]-[CX4]")  # 6 linear carbons as an example
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "No sufficient carbon chain detected"
    
    return True, "Contains carboxylic acid and additional ketone group typical of oxo fatty acids"