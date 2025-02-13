"""
Classifies: CHEBI:32952 amine
"""
"""
Classifies: CHEBI:32952 amine
"""
from rdkit import Chem

def is_amine(smiles: str):
    """
    Determines if a molecule is an amine based on its SMILES string.
    An amine is a compound formally derived from ammonia by replacing one, two, or three hydrogen atoms by hydrocarbyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the amine pattern: nitrogen bonded to 1, 2, or 3 carbon atoms
    amine_pattern = Chem.MolFromSmarts("[NX3;H0,H1,H2][CX4]")
    amine_matches = mol.GetSubstructMatches(amine_pattern)

    if len(amine_matches) == 0:
        return False, "No amine group found (nitrogen bonded to 1, 2, or 3 carbon atoms)"

    # Check if the nitrogen is bonded to at least one carbon (hydrocarbyl group)
    for match in amine_matches:
        nitrogen_idx = match[0]
        nitrogen_atom = mol.GetAtomWithIdx(nitrogen_idx)
        # Count the number of carbon neighbors
        carbon_neighbors = sum(1 for neighbor in nitrogen_atom.GetNeighbors() if neighbor.GetAtomicNum() == 6)
        if carbon_neighbors >= 1:
            return True, "Contains at least one nitrogen atom bonded to a hydrocarbyl group"

    return False, "No nitrogen atom bonded to a hydrocarbyl group found"