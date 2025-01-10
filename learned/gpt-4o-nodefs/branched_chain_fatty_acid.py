"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
from rdkit import Chem

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid based on its SMILES string.
    This includes having one or more branching points in a hydrocarbon chain with a terminal carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a branched-chain fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None

    # Define SMARTS patterns
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")  # Terminal carboxylic acid
    branch_pattern = Chem.MolFromSmarts("[CH3,CH2,CH](C)C")  # Simple branched chain pattern for non-ring carbon

    # Check for terminal carboxylic acid group (COOH)
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No terminal carboxylic acid group found"
    
    # Check for branching patterns in hydrocarbon chains
    if mol.HasSubstructMatch(branch_pattern):
        return True, "Contains a branched carbon chain with a terminal carboxylic acid group"

    # Consider more complex branching
    complex_branch_pattern = Chem.MolFromSmarts("[C;!R;!D1;!D2]([C;!R;D2,3])([C;!R;D2,3])")  # More complex branching pattern

    if mol.HasSubstructMatch(complex_branch_pattern):
        return True, "Contains a complex branched carbon chain with a terminal carboxylic acid group"

    # If no branching found, double check if could be due to known complex structures
    return False, "No suitable branching pattern found in the structure leading to a carboxylic group"