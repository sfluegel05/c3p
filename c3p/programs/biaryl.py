"""
Classifies: CHEBI:64459 biaryl
"""
from rdkit import Chem

def is_biaryl(smiles: str):
    """
    Determines if a molecule is a biaryl (contains two aromatic rings or ring systems, joined by a single bond).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a biaryl, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all aromatic rings
    aromatic_rings = [ring for ring in mol.GetRingInfo().AtomRings() if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)]

    if len(aromatic_rings) < 2:
        return False, "Less than two aromatic rings found"

    # Check for a single bond connecting two aromatic rings
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetIsAromatic() and end_atom.GetIsAromatic():
                return True, "Biaryl structure found: two aromatic rings connected by a single bond"

    return False, "No single bond connecting two aromatic rings found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:64459',
                          'name': 'biaryl',
                          'definition': 'An organic aromatic compound whose '
                                        'structure contains two aromatic rings '
                                        'or ring systems, joined to each other '
                                        'by a single bond.',
                          'parents': ['CHEBI:33659']},
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
    'num_true_positives': 148,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}