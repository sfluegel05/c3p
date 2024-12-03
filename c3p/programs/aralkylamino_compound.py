"""
Classifies: CHEBI:64365 aralkylamino compound
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_aralkylamino_compound(smiles: str):
    """
    Determines if a molecule is an aralkylamino compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aralkylamino compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if there is an aromatic ring present
    if not any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "No aromatic ring found"

    # Check for the presence of an amino group (NH2, NHR, or NR2)
    amino_group = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N':
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'H' or neighbor.GetSymbol() in ['C', 'N']:
                    amino_group = True
                    break
            if amino_group:
                break

    if not amino_group:
        return False, "No amino group found"

    # Check if the amino group is linked to an alkyl group
    amino_linked_to_alkyl = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N':
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C' and not neighbor.GetIsAromatic():
                    amino_linked_to_alkyl = True
                    break
            if amino_linked_to_alkyl:
                break

    if not amino_linked_to_alkyl:
        return False, "Amino group is not linked to an alkyl group"

    # Check if the alkyl group is linked to an aromatic ring
    alkyl_linked_to_aromatic = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and not atom.GetIsAromatic():
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIsAromatic():
                    alkyl_linked_to_aromatic = True
                    break
            if alkyl_linked_to_aromatic:
                break

    if not alkyl_linked_to_aromatic:
        return False, "Alkyl group is not linked to an aromatic ring"

    return True, "Molecule is an aralkylamino compound"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:64365',
                          'name': 'aralkylamino compound',
                          'definition': 'An organic amino compound in which an '
                                        'aminoalkyl group is linked to an '
                                        'arene.',
                          'parents': ['CHEBI:50047']},
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
    'num_true_positives': 18,
    'num_false_positives': 3,
    'num_true_negatives': 16,
    'num_false_negatives': 1,
    'precision': 0.8571428571428571,
    'recall': 0.9473684210526315,
    'f1': 0.9,
    'accuracy': None}