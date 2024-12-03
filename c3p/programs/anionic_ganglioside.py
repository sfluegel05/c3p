"""
Classifies: CHEBI:79346 anionic ganglioside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_anionic_ganglioside(smiles: str):
    """
    Determines if a molecule is an anionic ganglioside.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anionic ganglioside, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a neuraminic acid derivative (e.g., Neu5Ac, Neu5Gc)
    neuraminic_acid_substructures = [
        "C[C@H](O)[C@@H](NC=O)C(=O)[O-]",  # Simplified Neu5Ac pattern
        "C[C@H](O)[C@@H](NC=O)C(=O)[O-]"   # Simplified Neu5Gc pattern
    ]
    
    neuraminic_acid_found = False
    for substructure in neuraminic_acid_substructures:
        sub_mol = Chem.MolFromSmarts(substructure)
        if mol.HasSubstructMatch(sub_mol):
            neuraminic_acid_found = True
            break

    if not neuraminic_acid_found:
        return False, "No neuraminic acid derivative found"

    # Check for the presence of a carbohydrate chain
    carbohydrate_substructure = Chem.MolFromSmarts("C(CO)O")
    if not mol.HasSubstructMatch(carbohydrate_substructure):
        return False, "No carbohydrate chain found"

    # Check for deprotonated carboxy groups
    carboxy_group = Chem.MolFromSmarts("C(=O)[O-]")
    carboxy_matches = mol.GetSubstructMatches(carboxy_group)
    if len(carboxy_matches) < 1:
        return False, "No deprotonated carboxy groups found"

    return True, "Anionic ganglioside detected"

# Test cases
smiles_examples = [
    "O([C@H]1[C@H]([C@H](O[C@H]([C@@H]1O)O[C@@H]2[C@H](O[C@@H](OC[C@@H]([C@@H](*)O)NC(=O)*)[C@@H]([C@H]2O)O)CO)CO)O[C@H]3[C@@H]([C@H]([C@@H](O)[C@H](O3)CO)O[C@H]4[C@@H]([C@@H](O[C@]5(O[C@]([C@@H]([C@H](C5)O)NC(C)=O)([C@@H]([C@H](O[C@]6(O[C@]([C@@H]([C@H](C6)O)NC(C)=O)([C@@H]([C@@H](CO)O)O)[H])C(=O)[O-])CO)O)[H])C([O-])=O)[C@H]([C@H](O4)CO)O)O)NC(C)=O)[C@]7(O[C@H]([C@H](NC(=O)C)[C@H](C7)O)[C@@H]([C@H](O)CO)O)C(=O)[O-]"
]

for smiles in smiles_examples:
    result, reason = is_anionic_ganglioside(smiles)
    print(f"SMILES: {smiles}\nResult: {result}\nReason: {reason}\n")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:79346',
                          'name': 'anionic ganglioside',
                          'definition': 'Any carbohydrate acid derivative '
                                        'anion obtained by deprotonation of at '
                                        'least one of the neuraminosyl carboxy '
                                        'groups of a ganglioside.',
                          'parents': ['CHEBI:63551']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'SMILES: '
              'O([C@H]1[C@H]([C@H](O[C@H]([C@@H]1O)O[C@@H]2[C@H](O[C@@H](OC[C@@H]([C@@H](*)O)NC(=O)*)[C@@H]([C@H]2O)O)CO)CO)O[C@H]3[C@@H]([C@H]([C@@H](O)[C@H](O3)CO)O[C@H]4[C@@H]([C@@H](O[C@]5(O[C@]([C@@H]([C@H](C5)O)NC(C)=O)([C@@H]([C@H](O[C@]6(O[C@]([C@@H]([C@H](C6)O)NC(C)=O)([C@@H]([C@@H](CO)O)O)[H])C(=O)[O-])CO)O)[H])C([O-])=O)[C@H]([C@H](O4)CO)O)O)NC(C)=O)[C@]7(O[C@H]([C@H](NC(=O)C)[C@H](C7)O)[C@@H]([C@H](O)CO)O)C(=O)[O-]\n'
              'Result: False\n'
              'Reason: No neuraminic acid derivative found\n'
              '\n',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 34,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}