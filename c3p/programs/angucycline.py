"""
Classifies: CHEBI:48130 angucycline
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_angucycline(smiles: str):
    """
    Determines if a molecule is an angucycline.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an angucycline, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check for presence of benz[a]anthracene ring system (4 fused rings with specific pattern)
    benz_a_anthracene_pattern = Chem.MolFromSmarts('c1ccc2c(c1)ccc3c2ccc4c3cccc4')
    if not mol.HasSubstructMatch(benz_a_anthracene_pattern):
        return False, "No benz[a]anthracene ring system found"

    # Check for hydrolysable sugars (presence of glycosidic bonds)
    glycosidic_pattern = Chem.MolFromSmarts('O[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H]1O')
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No hydrolysable sugars found"

    return True, "Molecule is an angucycline"

# Example usage:
# result, reason = is_angucycline('O=C1C2=C3C(=CC=C2C(=O)C4=C1C=CC=C4OC)C(O)C(C)CC3=O')
# print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:48130',
                          'name': 'angucycline',
                          'definition': 'Polyketides produced by Actinomycetes '
                                        'which have structures based on the '
                                        'benz[a]anthracene ring system, '
                                        'several of which bear hydrolysable '
                                        'sugars.',
                          'parents': ['CHEBI:26188', 'CHEBI:51067']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 14,
    'num_false_negatives': 14,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}