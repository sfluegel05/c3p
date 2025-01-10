"""
Classifies: CHEBI:26333 prostaglandin
"""
"""
Classifies: CHEBI:26347 prostaglandin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin based on its SMILES string.
    A prostaglandin is a naturally occurring compound derived from the parent C20 acid, prostanoic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prostaglandin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a cyclopentane ring
    cyclopentane_pattern = Chem.MolFromSmarts("[C]1[C][C][C][C]1")
    if not mol.HasSubstructMatch(cyclopentane_pattern):
        return False, "No cyclopentane ring found"

    # Check for the presence of at least one carboxyl group (COOH)
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group found"

    # Check for the presence of at least one hydroxyl group (OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found"

    # Check for the presence of at least one double bond (C=C)
    double_bond_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    if not mol.HasSubstructMatch(double_bond_pattern):
        return False, "No double bond found"

    # Check molecular weight - prostaglandins typically have a molecular weight around 350-400 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 500:
        return False, f"Molecular weight {mol_wt:.2f} is outside the typical range for prostaglandins"

    # Count carbons - prostaglandins typically have 20 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18 or c_count > 22:
        return False, f"Number of carbons {c_count} is outside the typical range for prostaglandins"

    return True, "Contains cyclopentane ring, carboxyl group, hydroxyl group, and double bond, with appropriate molecular weight and carbon count"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26347',
                          'name': 'prostaglandin',
                          'definition': 'Naturally occurring compounds derived from the parent C20 acid, prostanoic acid.',
                          'parents': ['CHEBI:26347', 'CHEBI:26347']},
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