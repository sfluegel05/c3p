"""
Classifies: CHEBI:26199 polyprenol
"""
"""
Classifies: CHEBI:26195 polyprenol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    A polyprenol is a member of the class of prenols with the general formula H-[CH2C(Me)=CHCH2]nOH, 
    where the carbon skeleton is composed of more than one isoprene units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a hydroxyl group (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found"

    # Check for the presence of isoprene units (C5H8) with alternating double bonds
    isoprene_pattern = Chem.MolFromSmarts("[CH2][CH]=[CH][CH2]")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 2:
        return False, f"Found {len(isoprene_matches)} isoprene units, need at least 2"

    # Check if the molecule has more than one isoprene unit
    # Each isoprene unit has 5 carbons, so the total number of carbons should be at least 10
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
        return False, "Too few carbons for polyprenol (need at least 10)"

    # Check for alternating double bonds in the chain
    # This is a simplified check and may not cover all cases
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bond_count < 2:
        return False, "Not enough double bonds for polyprenol structure"

    return True, "Contains a hydroxyl group and a chain of more than one isoprene units with alternating double bonds"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26195',
                          'name': 'polyprenol',
                          'definition': 'Any member of the class of prenols '
                                        'possessing the general formula '
                                        'H-[CH2C(Me)=CHCH2]nOH in which the '
                                        'carbon skeleton is composed of more '
                                        'than one isoprene units.',
                          'parents': ['CHEBI:26195', 'CHEBI:26195']},
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