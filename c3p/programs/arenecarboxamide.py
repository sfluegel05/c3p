"""
Classifies: CHEBI:22645 arenecarboxamide
"""
from rdkit import Chem

def is_arenecarboxamide(smiles: str):
    """
    Determines if a molecule is an arenecarboxamide (a monocarboxylic acid amide in which the amide linkage
    is bonded directly to an arene ring system).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an arenecarboxamide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of an amide group
    amide_pattern = Chem.MolFromSmarts('C(=O)N')
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide group found"

    # Check if the amide group is bonded directly to an arene ring
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    for amide_match in amide_matches:
        amide_carbon, amide_nitrogen = amide_match[0], amide_match[1]
        amide_carbon_atom = mol.GetAtomWithIdx(amide_carbon)
        amide_nitrogen_atom = mol.GetAtomWithIdx(amide_nitrogen)

        # Check if the carbon or nitrogen of the amide is bonded to an aromatic ring
        for neighbor in amide_carbon_atom.GetNeighbors():
            if neighbor.GetIsAromatic():
                return True, "Amide group is bonded directly to an arene ring"
        for neighbor in amide_nitrogen_atom.GetNeighbors():
            if neighbor.GetIsAromatic():
                return True, "Amide group is bonded directly to an arene ring"

    return False, "Amide group is not bonded directly to an arene ring"

# Example usage
smiles = 'CCOC1=CC=C(C=C1)NC(=O)C2=CC=C(C=C2)CNC3=C(C(=O)C3=O)NC4CCCCC4'
print(is_arenecarboxamide(smiles))  # Should return (True, "Amide group is bonded directly to an arene ring")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:22645',
                          'name': 'arenecarboxamide',
                          'definition': 'A monocarboxylic acid amide in which '
                                        'the amide linkage is bonded directly '
                                        'to an arene ring system.',
                          'parents': ['CHEBI:29347', 'CHEBI:62733']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(True, 'Amide group is bonded directly to an arene ring')\n",
    'num_true_positives': 102,
    'num_false_positives': 5,
    'num_true_negatives': 15,
    'num_false_negatives': 0,
    'precision': 0.9532710280373832,
    'recall': 1.0,
    'f1': 0.9760765550239235,
    'accuracy': None}