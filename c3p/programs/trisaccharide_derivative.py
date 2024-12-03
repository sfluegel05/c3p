"""
Classifies: CHEBI:63571 trisaccharide derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_trisaccharide_derivative(smiles: str):
    """
    Determines if a molecule is a trisaccharide derivative.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trisaccharide derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least three sugar rings (5 or 6-membered rings with oxygen)
    rings = mol.GetRingInfo()
    sugar_rings = 0
    for ring in rings.AtomRings():
        if 5 <= len(ring) <= 6:
            if any(mol.GetAtomWithIdx(i).GetSymbol() == 'O' for i in ring):
                sugar_rings += 1

    if sugar_rings < 3:
        return False, "Not enough sugar rings found"

    # Check for glycosidic linkages (O-glycosidic bonds)
    glycosidic_linkages = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if (begin_atom.GetSymbol() == 'O' and end_atom.GetSymbol() == 'C') or (begin_atom.GetSymbol() == 'C' and end_atom.GetSymbol() == 'O'):
                glycosidic_linkages += 1

    if glycosidic_linkages < 2:
        return False, "Not enough glycosidic linkages found"

    # Check for additional functional groups indicating derivatization
    derivatized = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in ['C', 'H', 'O']:
            derivatized = True
            break

    if derivatized:
        return True, "Trisaccharide derivative with additional functional groups"
    else:
        return True, "Trisaccharide without additional functional groups"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:63571',
                          'name': 'trisaccharide derivative',
                          'definition': 'An oligosaccharide derivative that is '
                                        'formally obtained from a '
                                        'trisaccharide.',
                          'parents': ['CHEBI:63563']},
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
    'num_true_positives': 41,
    'num_false_positives': 17,
    'num_true_negatives': 3,
    'num_false_negatives': 1,
    'precision': 0.7068965517241379,
    'recall': 0.9761904761904762,
    'f1': 0.8199999999999998,
    'accuracy': None}