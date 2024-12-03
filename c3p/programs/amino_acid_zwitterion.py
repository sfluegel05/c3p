"""
Classifies: CHEBI:35238 amino acid zwitterion
"""
from rdkit import Chem

def is_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an amino acid zwitterion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino acid zwitterion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    found_pos_charge = False
    found_neg_charge = False
    found_carboxyl_group = False
    found_amino_group = False

    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() == 1 and atom.GetSymbol() == 'N':
            found_pos_charge = True
        elif atom.GetFormalCharge() == -1 and atom.GetSymbol() == 'O':
            found_neg_charge = True

    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if (atom1.GetSymbol() == 'C' and atom2.GetSymbol() == 'O' and atom2.GetFormalCharge() == -1) or \
           (atom2.GetSymbol() == 'C' and atom1.GetSymbol() == 'O' and atom1.GetFormalCharge() == -1):
            found_carboxyl_group = True
        if (atom1.GetSymbol() == 'C' and atom2.GetSymbol() == 'N' and atom2.GetFormalCharge() == 1) or \
           (atom2.GetSymbol() == 'C' and atom1.GetSymbol() == 'N' and atom1.GetFormalCharge() == 1):
            found_amino_group = True

    if found_pos_charge and found_neg_charge and found_carboxyl_group and found_amino_group:
        return True, "Molecule has both positively charged amino group and negatively charged carboxyl group"
    else:
        return False, "Molecule does not have both positively charged amino group and negatively charged carboxyl group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35238',
                          'name': 'amino acid zwitterion',
                          'definition': 'The zwitterionic form of an amino '
                                        'acid having a negatively charged '
                                        'carboxyl group and a positively '
                                        'charged amino group.',
                          'parents': ['CHEBI:27369']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 42,
    'num_false_positives': 6,
    'num_true_negatives': 14,
    'num_false_negatives': 0,
    'precision': 0.875,
    'recall': 1.0,
    'f1': 0.9333333333333333,
    'accuracy': None}