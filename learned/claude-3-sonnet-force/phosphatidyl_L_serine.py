"""
Classifies: CHEBI:18303 phosphatidyl-L-serine
"""
"""
Classifies: CHEBI:26811 phosphatidyl-L-serine
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a phosphatidyl-L-serine based on its SMILES string.
    Phosphatidyl-L-serine is an aminophospholipid with a phosphatidyl group esterified
    to the hydroxy group of L-serine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphatidyl group
    phosphatidyl_pattern = Chem.MolFromSmarts("[P](=[O])([O-])O[C@@H]([C@@H](O)[C@@H](O)[C@@H](O)[C@@H](O)[C@@H](O)[C@@H](O)[C@@H](O)[C@@H](O)[C@@H](O)[C@@H](O)[C@@H](O)[C@@H](O)[C@@H](O)[C@@H](O)[C@@H](O)[C@@H](O)[C@@H](O)[C@@H](O)[C@@H](O)CO)")
    if not mol.HasSubstructMatch(phosphatidyl_pattern):
        return False, "Missing phosphatidyl group"

    # Check for serine residue
    serine_pattern = Chem.MolFromSmarts("C(N)C(=O)O")
    if not mol.HasSubstructMatch(serine_pattern):
        return False, "Missing serine residue"

    # Check for ester linkage between phosphatidyl group and serine
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "Missing ester linkage"

    # Check for fatty acid chains
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, "Missing fatty acid chains"

    # Check molecular weight (typical range: 700-900 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700 or mol_wt > 900:
        return False, "Molecular weight outside typical range"

    return True, "Contains phosphatidyl group esterified to serine, with fatty acid chains"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:26811',
        'name': 'phosphatidyl-L-serine',
        'definition': 'A class of aminophospholipids in which a phosphatidyl group is esterified to the hydroxy group of serine.',
        'parents': ['CHEBI:26865', 'CHEBI:18059']
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    },
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 30,
    'num_false_positives': 0,
    'num_true_negatives': 182418,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': 1.0
}