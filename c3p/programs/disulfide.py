"""
Classifies: CHEBI:48343 disulfide
"""
from rdkit import Chem

def is_disulfide(smiles: str):
    """
    Determines if a molecule is a disulfide (contains RSSR structure).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a disulfide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the disulfide bond pattern
    disulfide_pattern = Chem.MolFromSmarts("SS")

    # Find the disulfide bonds in the molecule
    matches = mol.GetSubstructMatches(disulfide_pattern)

    if not matches:
        return False, "No disulfide bonds found"

    # Check if each sulfur in the disulfide bond is connected to carbon (R group)
    for match in matches:
        s1, s2 = match
        s1_neighbors = mol.GetAtomWithIdx(s1).GetNeighbors()
        s2_neighbors = mol.GetAtomWithIdx(s2).GetNeighbors()

        if not any(neighbor.GetSymbol() == 'C' for neighbor in s1_neighbors if neighbor.GetIdx() != s2):
            return False, "Sulfur in disulfide bond is not connected to a carbon (R group)"
        if not any(neighbor.GetSymbol() == 'C' for neighbor in s2_neighbors if neighbor.GetIdx() != s1):
            return False, "Sulfur in disulfide bond is not connected to a carbon (R group)"

    return True, "Disulfide bond(s) found with appropriate R groups"

# Example usage
smiles = "NC(C(=O)NCC(=O)O)CSSCC(C(=O)NCC(=O)O)N"
result, reason = is_disulfide(smiles)
print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:48343',
                          'name': 'disulfide',
                          'definition': "Compounds of structure RSSR'.",
                          'parents': ['CHEBI:26835']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'True Disulfide bond(s) found with appropriate R groups\n',
    'num_true_positives': 15,
    'num_false_positives': 0,
    'num_true_negatives': 15,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}