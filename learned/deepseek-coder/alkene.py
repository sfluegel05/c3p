"""
Classifies: CHEBI:32878 alkene
"""
"""
Classifies: CHEBI:32878 alkene
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkene(smiles: str):
    """
    Determines if a molecule is an alkene based on its SMILES string.
    An alkene is an acyclic hydrocarbon with exactly one carbon-carbon double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add hydrogens to the molecule
    mol = Chem.AddHs(mol)

    # Check if the molecule is acyclic
    if not rdMolDescriptors.CalcNumRings(mol) == 0:
        return False, "Molecule contains rings"

    # Count the number of carbon-carbon double bonds
    double_bond_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE and \
           bond.GetBeginAtom().GetAtomicNum() == 6 and \
           bond.GetEndAtom().GetAtomicNum() == 6:
            double_bond_count += 1

    if double_bond_count != 1:
        return False, f"Found {double_bond_count} carbon-carbon double bonds, need exactly 1"

    # Check if the molecule is a hydrocarbon (only contains C and H)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [1, 6]:
            return False, "Molecule contains non-carbon and non-hydrogen atoms"

    # Check the general formula CnH2n
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    if h_count != 2 * c_count:
        return False, f"Hydrogen count ({h_count}) does not match 2 * carbon count ({2 * c_count})"

    return True, "Acyclic hydrocarbon with exactly one carbon-carbon double bond"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:32878',
                          'name': 'alkene',
                          'definition': 'An acyclic branched or unbranched '
                                        'hydrocarbon having one carbon-carbon '
                                        'double bond and the general formula '
                                        'CnH2n. Acyclic branched or unbranched '
                                        'hydrocarbons having more than one '
                                        'double bond are alkadienes, '
                                        'alkatrienes, etc.',
                          'parents': ['CHEBI:32877', 'CHEBI:33281']},
    'config': {   'llm_model_name': 'lbl/claude-sonnet',
                  'f1_threshold': 0.8,
                  'max_attempts': 5,
                  'max_positive_instances': None,
                  'max_positive_to_test': None,
                  'max_negative_to_test': None,
                  'max_positive_in_prompt': 50,
                  'max_negative_in_prompt': 20,
                  'max_instances_in_prompt': 100,
                  'test_proportion': 0.1},
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199}