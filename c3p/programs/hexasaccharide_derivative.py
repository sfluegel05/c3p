"""
Classifies: CHEBI:63565 hexasaccharide derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_hexasaccharide_derivative(smiles: str):
    """
    Determines if a molecule is a hexasaccharide derivative (an oligosaccharide derivative that is formally obtained from a hexasaccharide).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexasaccharide derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check for at least six rings (hexasaccharide)
    if len(rings.AtomRings()) < 6:
        return False, "Less than six rings found"

    # Check if all rings are 5 or 6-membered rings (common for sugars)
    for ring in rings.AtomRings():
        if len(ring) != 5 and len(ring) != 6:
            return False, "Ring with size other than 5 or 6 found"

    # Check if each ring contains oxygen atoms (common for sugars)
    for ring in rings.AtomRings():
        ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if not any(atom.GetSymbol() == 'O' for atom in ring_atoms):
            return False, "Ring without oxygen atom found"

    # Check for glycosidic linkages (C-O-C bonds between rings)
    glycosidic_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if (begin_atom.GetSymbol() == 'O' and end_atom.GetSymbol() == 'C') or (begin_atom.GetSymbol() == 'C' and end_atom.GetSymbol() == 'O'):
                glycosidic_bonds += 1

    if glycosidic_bonds < 5:
        return False, "Less than five glycosidic linkages found"

    return True, "Hexasaccharide derivative found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:63565',
                          'name': 'hexasaccharide derivative',
                          'definition': 'An oligosaccharide derivative that is '
                                        'formally obtained from a '
                                        'hexasaccharide.',
                          'parents': ['CHEBI:63563']},
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
    'num_true_positives': 12,
    'num_false_positives': 5,
    'num_true_negatives': 7,
    'num_false_negatives': 0,
    'precision': 0.7058823529411765,
    'recall': 1.0,
    'f1': 0.8275862068965517,
    'accuracy': None}