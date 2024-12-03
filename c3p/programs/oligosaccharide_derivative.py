"""
Classifies: CHEBI:63563 oligosaccharide derivative
"""
from rdkit import Chem


def is_oligosaccharide_derivative(smiles: str):
    """
    Determines if a molecule is an oligosaccharide derivative.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oligosaccharide derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carbohydrate-like structure
    if not any(atom.GetSymbol() == 'O' for atom in mol.GetAtoms()):
        return False, "No oxygen atoms found"

    # Check for glycosidic bonds (O-glycosidic linkages)
    glycosidic_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if (begin_atom.GetSymbol() == 'O' and end_atom.GetSymbol() == 'C') or (begin_atom.GetSymbol() == 'C' and end_atom.GetSymbol() == 'O'):
                glycosidic_bonds += 1

    if glycosidic_bonds == 0:
        return False, "No glycosidic bonds found"

    # Check for acetamido groups (common in oligosaccharide derivatives)
    acetamido_groups = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N':
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C':
                    for sub_neighbor in neighbor.GetNeighbors():
                        if sub_neighbor.GetSymbol() == 'O':
                            acetamido_groups += 1

    if acetamido_groups == 0:
        return False, "No acetamido groups found"

    return True, "Molecule is an oligosaccharide derivative"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:63563',
                          'name': 'oligosaccharide derivative',
                          'definition': 'A carbohydrate derivative that is '
                                        'formally obtained from an '
                                        'oligosaccharide.',
                          'parents': ['CHEBI:167559', 'CHEBI:63299']},
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
    'num_true_positives': 185,
    'num_false_positives': 13,
    'num_true_negatives': 7,
    'num_false_negatives': 36,
    'precision': 0.9343434343434344,
    'recall': 0.8371040723981901,
    'f1': 0.8830548926014321,
    'accuracy': None}