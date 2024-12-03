"""
Classifies: CHEBI:36233 disaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide (a compound in which two monosaccharides are joined by a glycosidic bond).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a disaccharide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycosidic bond
    glycosidic_bond_found = False
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetSymbol() == 'O' and end_atom.GetSymbol() == 'C':
                # Check if both atoms are part of a ring (indicating a glycosidic bond between sugar rings)
                if begin_atom.IsInRing() and end_atom.IsInRing():
                    glycosidic_bond_found = True
                    break

    if not glycosidic_bond_found:
        return False, "No glycosidic bond found"

    # Check for two monosaccharides
    ring_count = sum(1 for ring in mol.GetRingInfo().AtomRings() if len(ring) == 5 or len(ring) == 6)
    if ring_count != 2:
        return False, f"Expected 2 monosaccharide rings, found {ring_count}"

    return True, "Valid disaccharide with glycosidic bond"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36233',
                          'name': 'disaccharide',
                          'definition': 'A compound in which two '
                                        'monosaccharides are joined by a '
                                        'glycosidic bond.',
                          'parents': ['CHEBI:50699']},
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
    'num_true_positives': 32,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 4,
    'precision': 1.0,
    'recall': 0.8888888888888888,
    'f1': 0.9411764705882353,
    'accuracy': None}