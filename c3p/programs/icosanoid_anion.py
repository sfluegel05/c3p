"""
Classifies: CHEBI:62937 icosanoid anion
"""
from rdkit import Chem

def is_icosanoid_anion(smiles: str):
    """
    Determines if a molecule is an icosanoid anion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an icosanoid anion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylate group (COO-)
    carboxylate_group = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 3:
                oxygen_count = 0
                negative_charge_count = 0
                for neighbor in neighbors:
                    if neighbor.GetSymbol() == 'O':
                        oxygen_count += 1
                        if neighbor.GetFormalCharge() == -1:
                            negative_charge_count += 1
                if oxygen_count == 2 and negative_charge_count == 1:
                    carboxylate_group = True
                    break

    if not carboxylate_group:
        return False, "No carboxylate group found"

    # Check for chain length of at least 20 carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    if carbon_count < 20:
        return False, "Carbon chain length is less than 20"

    # Check for presence of at least three double bonds (C=C)
    double_bond_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            double_bond_count += 1

    if double_bond_count < 3:
        return False, "Less than three double bonds found"

    return True, "Molecule is an icosanoid anion"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:62937',
                          'name': 'icosanoid anion',
                          'definition': 'The carboxylic acid anion that is the '
                                        'conjugate base of an icosanoid, '
                                        'formed when the carboxy group is '
                                        'deprotonated.',
                          'parents': ['CHEBI:29067']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 4,
    'success': True,
    'best': True,
    'error': '[02:03:16] SMILES Parse Error: syntax error while parsing: '
             'C(CC/C=C\\C=C\x01/O[C@H]1C/C=C\\CC)CCCCC(=O)[O-]\n'
             '[02:03:16] SMILES Parse Error: Failed parsing SMILES '
             "'C(CC/C=C\\C=C\x01/O[C@H]1C/C=C\\CC)CCCCC(=O)[O-]' for input: "
             "'C(CC/C=C\\C=C\x01/O[C@H]1C/C=C\\CC)CCCCC(=O)[O-]'\n",
    'stdout': '',
    'num_true_positives': 27,
    'num_false_positives': 2,
    'num_true_negatives': 18,
    'num_false_negatives': 0,
    'precision': 0.9310344827586207,
    'recall': 1.0,
    'f1': 0.9642857142857143,
    'accuracy': None}