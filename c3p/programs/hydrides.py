"""
Classifies: CHEBI:33692 hydrides
"""
from rdkit import Chem

def is_hydrides(smiles: str):
    """
    Determines if a molecule is a hydride (a chemical compound of hydrogen with other chemical elements).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydride, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule contains hydrogen
    contains_hydrogen = any(atom.GetSymbol() == 'H' for atom in mol.GetAtoms())
    if not contains_hydrogen:
        return False, "No hydrogen atoms found"

    # Check if the molecule contains other elements besides hydrogen
    contains_other_elements = any(atom.GetSymbol() != 'H' for atom in mol.GetAtoms())
    if not contains_other_elements:
        return False, "Only hydrogen atoms found"

    return True, "Contains hydrogen and other chemical elements"

# Example usage:
# print(is_hydrides("C(CCCCCCCCCCC)(CCCCC)C"))  # 6-Methylheptadecane
# print(is_hydrides("[N-][H]"))  # hydridonitrate(.1-)
# print(is_hydrides("[H][B+]([H])([H])[H]"))  # boronium


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33692',
                          'name': 'hydrides',
                          'definition': 'Hydrides are chemical compounds of '
                                        'hydrogen with other chemical '
                                        'elements.',
                          'parents': ['CHEBI:33608', 'CHEBI:37577']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[15:26:13] SMILES Parse Error: unclosed ring for input: '
             "'CC(C)C1CC\\C(C)=C\\CCC(=C)\\C=C\x01'\n"
             '[15:26:13] SMILES Parse Error: unclosed ring for input: '
             "'CC(C)C1=C/C=C(C)/CC\\C=C(C)\\CC\x01'\n"
             '[15:26:13] Explicit valence for atom # 0 B, 5, is greater than '
             'permitted\n'
             '[15:26:13] Explicit valence for atom # 0 B, 6, is greater than '
             'permitted\n'
             '[15:26:13] WARNING: not removing hydrogen atom without '
             'neighbors\n',
    'stdout': '',
    'num_true_positives': 2,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 207,
    'precision': 0.6666666666666666,
    'recall': 0.009569377990430622,
    'f1': 0.018867924528301883,
    'accuracy': None}