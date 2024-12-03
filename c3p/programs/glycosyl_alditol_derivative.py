"""
Classifies: CHEBI:63424 glycosyl alditol derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_glycosyl_alditol_derivative(smiles: str):
    """
    Determines if a molecule is a glycosyl alditol derivative.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosyl alditol derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of glycosidic bond (C-O-C linkage)
    glycosidic_bond = False
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if (begin_atom.GetSymbol() == 'C' and end_atom.GetSymbol() == 'O') or (begin_atom.GetSymbol() == 'O' and end_atom.GetSymbol() == 'C'):
                glycosidic_bond = True
                break

    if not glycosidic_bond:
        return False, "No glycosidic bond found"

    # Check for the presence of alditol (polyhydroxylated structure)
    alditol = False
    hydroxyl_count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and atom.GetDegree() == 1:
            neighbor = atom.GetNeighbors()[0]
            if neighbor.GetSymbol() == 'C' and neighbor.GetDegree() == 2:
                hydroxyl_count += 1
                if hydroxyl_count >= 3:
                    alditol = True
                    break

    if not alditol:
        return False, "No alditol structure found"

    return True, "Glycosyl alditol derivative found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:63424',
                          'name': 'glycosyl alditol derivative',
                          'definition': 'An alditol derivative that is '
                                        'formally obtained from a glycosyl '
                                        'alditol.',
                          'parents': ['CHEBI:63423']},
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
    'num_true_positives': 2,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 21,
    'precision': 1.0,
    'recall': 0.08695652173913043,
    'f1': 0.16,
    'accuracy': None}