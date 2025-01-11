"""
Classifies: CHEBI:23824 diol
"""
from rdkit import Chem

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol based on its SMILES string.
    A diol is characterized by having two hydroxyl groups (-OH) attached to different carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define hydroxyl group SMARTS pattern
    oh_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
    
    # Find matches for hydroxyl groups
    hydroxyl_matches = mol.GetSubstructMatches(oh_pattern)
    
    # Check if there are exactly two hydroxyl groups
    if len(hydroxyl_matches) != 2:
        return False, f"Found {len(hydroxyl_matches)} hydroxyl groups, need exactly 2"

    # Check that each OH is on a different carbon
    carbon_indices = [match[0] for match in hydroxyl_matches]
    if len(set(carbon_indices)) != 2:
        return False, "Both hydroxyl groups attached to the same carbon atom"

    return True, "Contains two hydroxyl groups on different carbon atoms"

# Example usage with SMILES strings
examples = [
    "[C@@H](CC)(C(/C=C/C=C/C=C\\C/C=C\\C/C=C\\CCCC(=O)O)O)O",
    "C[C@H](CCC=C(C)C)[C@H]1CC[C@H]2C3=C(CC[C@]12C)[C@@]1(C)CC[C@H](O)[C@@](C)(CO)[C@@H]1CC3"
]

for smiles in examples:
    is_diol_result, reason = is_diol(smiles)
    print(f"SMILES: {smiles}\nIs Diol: {is_diol_result}\nReason: {reason}\n")