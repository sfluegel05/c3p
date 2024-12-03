"""
Classifies: CHEBI:18379 nitrile
"""
from rdkit import Chem

def is_nitrile(smiles: str):
    """
    Determines if a molecule is a nitrile (RC#N).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrile, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Search for nitrile group (C#N)
    nitrile_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 7:  # Nitrogen atom
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                    if bond.GetBondType() == Chem.rdchem.BondType.TRIPLE:
                        nitrile_found = True
                        break
            if nitrile_found:
                break

    if nitrile_found:
        return True, "Nitrile group (C#N) found"
    else:
        return False, "No nitrile group (C#N) found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:18379',
                          'name': 'nitrile',
                          'definition': 'A compound having the structure RC#N; '
                                        'thus a C-substituted derivative of '
                                        'hydrocyanic acid, HC#N. In systematic '
                                        'nomenclature, the suffix nitrile '
                                        'denotes the triply bound #N atom, not '
                                        'the carbon atom attached to it.',
                          'parents': ['CHEBI:23424']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[19:55:51] SMILES Parse Error: syntax error while parsing: '
             'C/C=C\x01/C2CC[N@+]3(CC[C@]4(C5=CC=CC=C5NC4=C2C(=O)OC)[C@@]13[H])[O-]\n'
             '[19:55:51] SMILES Parse Error: Failed parsing SMILES '
             "'C/C=C\x01/C2CC[N@+]3(CC[C@]4(C5=CC=CC=C5NC4=C2C(=O)OC)[C@@]13[H])[O-]' "
             'for input: '
             "'C/C=C\x01/C2CC[N@+]3(CC[C@]4(C5=CC=CC=C5NC4=C2C(=O)OC)[C@@]13[H])[O-]'\n",
    'stdout': '',
    'num_true_positives': 65,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}