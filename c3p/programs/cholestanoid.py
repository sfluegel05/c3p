"""
Classifies: CHEBI:50401 cholestanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_cholestanoid(smiles: str):
    """
    Determines if a molecule is a cholestanoid (Any steroid based on a cholestane skeleton and its derivatives).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cholestanoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check for the presence of at least four rings
    if len(rings.AtomRings()) < 4:
        return False, "Less than four rings found"

    # Check for the presence of the steroid backbone (cyclopenta[a]phenanthrene structure)
    steroid_rings = [ring for ring in rings.AtomRings() if len(ring) == 6 or len(ring) == 5]
    if len(steroid_rings) < 4:
        return False, "Steroid backbone not found"

    # Check for the presence of substituents characteristic of cholestanoids
    cholestane_skeleton = Chem.MolFromSmiles("C[C@H](CCCC(C)C)C1CCC2C3CCC4CC[C@]4(C)C3CCC12C")
    if not mol.HasSubstructMatch(cholestane_skeleton):
        return False, "Cholestane skeleton not found"

    return True, "Cholestanoid structure identified"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:50401',
                          'name': 'cholestanoid',
                          'definition': 'Any steroid based on a cholestane '
                                        'skeleton and its derivatives.',
                          'parents': ['CHEBI:35341']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[00:55:49] SMILES Parse Error: syntax error while parsing: '
             'C/1(\\[C@]2([C@](CCC1)([C@](CC2)([C@@](CCCC(C)C)(C)[H])[H])C)[H])=C/C=C\x03/C[C@H](CCC3=C)F\n'
             '[00:55:49] SMILES Parse Error: Failed parsing SMILES '
             "'C/1(\\[C@]2([C@](CCC1)([C@](CC2)([C@@](CCCC(C)C)(C)[H])[H])C)[H])=C/C=C\x03/C[C@H](CCC3=C)F' "
             'for input: '
             "'C/1(\\[C@]2([C@](CCC1)([C@](CC2)([C@@](CCCC(C)C)(C)[H])[H])C)[H])=C/C=C\x03/C[C@H](CCC3=C)F'\n",
    'stdout': '',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 29,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}