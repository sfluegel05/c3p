"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
from rdkit import Chem

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid based on its SMILES string.
    Branched-chain fatty acids have one or more carbon branching points in a hydrocarbon backbone with a terminal carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a branched-chain fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for a terminal carboxylic acid group (as many false positives may have complex rings or be central)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for branching based strictly on non-ring carbon atoms
    is_branched = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and not atom.IsInRing():  # Non-ring carbon
            carbon_neighbors = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 6)
            if carbon_neighbors > 2:  # More than 2 carbon neighbors indicate branching
                is_branched = True
                break
    
    if not is_branched:
        return False, "No branching found in a non-ring carbon chain"

    return True, "Contains terminal carboxylic acid group and branch points in non-ring carbon chain"