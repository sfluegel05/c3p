"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid (a polyunsaturated fatty acid that contains three double bonds).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trienoic fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxyl_group = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetTotalNumHs() == 0:
            neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
            if 'O' in neighbors and 'O' in neighbors:
                carboxyl_group = True
                break

    if not carboxyl_group:
        return False, "No carboxylic acid group found"

    # Check for three double bonds
    double_bond_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            double_bond_count += 1

    if double_bond_count != 3:
        return False, f"Incorrect number of double bonds: {double_bond_count}"

    return True, "Molecule is a trienoic fatty acid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:73155',
                          'name': 'trienoic fatty acid',
                          'definition': 'Any polyunsaturated fatty acid that '
                                        'contains three double bonds.',
                          'parents': ['CHEBI:26208', 'CHEBI:53339']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[02:34:56] SMILES Parse Error: syntax error while parsing: '
             'C(CC/C=C\\C=C\x01/O[C@H]1C/C=C\\CC)CCCCC(=O)O\n'
             '[02:34:56] SMILES Parse Error: Failed parsing SMILES '
             "'C(CC/C=C\\C=C\x01/O[C@H]1C/C=C\\CC)CCCCC(=O)O' for input: "
             "'C(CC/C=C\\C=C\x01/O[C@H]1C/C=C\\CC)CCCCC(=O)O'\n",
    'stdout': '',
    'num_true_positives': 1,
    'num_false_positives': 2,
    'num_true_negatives': 15,
    'num_false_negatives': 16,
    'precision': 0.3333333333333333,
    'recall': 0.058823529411764705,
    'f1': 0.1,
    'accuracy': None}