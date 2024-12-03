"""
Classifies: CHEBI:84055 pentose phosphate
"""
from rdkit import Chem

def is_pentose_phosphate(smiles: str):
    """
    Determines if a molecule is a pentose phosphate (phosphate derivative of pentose).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pentose phosphate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphate groups (P(=O)(O)O)
    phosphate_groups = mol.GetSubstructMatches(Chem.MolFromSmarts('P(=O)(O)O'))
    if not phosphate_groups:
        return False, "No phosphate groups found"

    # Check for pentose structure (5 carbon sugar)
    pentose_found = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and len(atom.GetNeighbors()) == 4:
            carbon_neighbors = [neighbor.GetSymbol() for neighbor in atom.GetNeighbors()]
            if carbon_neighbors.count('C') == 3 and carbon_neighbors.count('O') == 1:
                pentose_found = True
                break

    if not pentose_found:
        return False, "No pentose structure found"

    return True, "Molecule is a pentose phosphate"

# Examples
print(is_pentose_phosphate("NC(=O)c1ncn([C@@H]2O[C@H](COP(O)(O)=O)[C@@H](O)[C@H]2O)c1N"))  # AICA ribonucleotide
print(is_pentose_phosphate("Nc1ncnc2n(cnc12)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@@H]3OP(O)(=O)O[C@@H]3[C@@H]2O)[C@@H](O)[C@H]1O"))  # ADP-alpha-D-ribose 1'',2''-cyclic phosphate


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:84055',
                          'name': 'pentose phosphate',
                          'definition': 'Any phospho sugar that is the '
                                        'phosphate derivative of pentose.',
                          'parents': ['CHEBI:33447', 'CHEBI:63409']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(False, 'No pentose structure found')\n"
              "(False, 'No pentose structure found')\n",
    'num_true_positives': 0,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 38,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}