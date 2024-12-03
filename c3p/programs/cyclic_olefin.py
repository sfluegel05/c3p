"""
Classifies: CHEBI:33642 cyclic olefin
"""
from rdkit import Chem

def is_cyclic_olefin(smiles: str):
    """
    Determines if a molecule is a cyclic olefin (cyclic hydrocarbon with any number of double bonds).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclic olefin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    ring_info = mol.GetRingInfo()
    
    # Check if the molecule is cyclic
    if not ring_info.NumRings():
        return False, "Molecule is not cyclic"

    # Check for the presence of double bonds in the rings
    for ring in ring_info.BondRings():
        if any(mol.GetBondWithIdx(bond_idx).GetBondType() == Chem.rdchem.BondType.DOUBLE for bond_idx in ring):
            return True, "Molecule is a cyclic olefin"

    return False, "No double bonds found in the cyclic structure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33642',
                          'name': 'cyclic olefin',
                          'definition': 'The inclusive term for any cyclic '
                                        'hydrocarbon having any number of '
                                        'double bonds.',
                          'parents': ['CHEBI:33641', 'CHEBI:33654']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 21,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 0,
    'precision': 0.9545454545454546,
    'recall': 1.0,
    'f1': 0.9767441860465117,
    'accuracy': None}