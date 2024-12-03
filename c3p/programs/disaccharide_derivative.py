"""
Classifies: CHEBI:63353 disaccharide derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_disaccharide_derivative(smiles: str):
    """
    Determines if a molecule is a disaccharide derivative.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a disaccharide derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule contains at least two sugar rings
    ring_info = mol.GetRingInfo()
    sugar_rings = 0

    for ring in ring_info.AtomRings():
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if all(atom.GetSymbol() in ['C', 'O'] for atom in atoms) and len(ring) in [5, 6]:
            if any(atom.GetSymbol() == 'O' for atom in atoms):
                sugar_rings += 1

    if sugar_rings < 2:
        return False, "Less than two sugar rings found"

    # Check for glycosidic bond
    glycosidic_bond = False
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if (begin_atom.GetSymbol() == 'O' and end_atom.GetSymbol() == 'C') or (begin_atom.GetSymbol() == 'C' and end_atom.GetSymbol() == 'O'):
                if begin_atom.IsInRing() and end_atom.IsInRing():
                    glycosidic_bond = True
                    break

    if not glycosidic_bond:
        return False, "No glycosidic bond found"

    return True, "Disaccharide derivative"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:63353',
                          'name': 'disaccharide derivative',
                          'definition': 'A carbohydrate derivative that is '
                                        'formally obtained from a '
                                        'disaccharide.',
                          'parents': ['CHEBI:63299']},
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
    'num_true_positives': 69,
    'num_false_positives': 5,
    'num_true_negatives': 15,
    'num_false_negatives': 0,
    'precision': 0.9324324324324325,
    'recall': 1.0,
    'f1': 0.965034965034965,
    'accuracy': None}