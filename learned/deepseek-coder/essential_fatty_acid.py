"""
Classifies: CHEBI:59549 essential fatty acid
"""
"""
Classifies: CHEBI:59554 essential fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is an essential fatty acid based on its SMILES string.
    Essential fatty acids are polyunsaturated fatty acids with multiple cis double bonds,
    a carboxyl group, and a chain length typically between 18 and 36 carbons.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an essential fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxyl group (-COOH)
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group found"

    # Check for polyunsaturated fatty acid structure (multiple cis double bonds)
    # Pattern for cis double bonds separated by a methylene group
    polyunsaturated_pattern = Chem.MolFromSmarts("[CX4][CX4]/[CX4]=[CX4]")
    polyunsaturated_matches = mol.GetSubstructMatches(polyunsaturated_pattern)
    if len(polyunsaturated_matches) < 2:
        return False, "Not enough cis double bonds for polyunsaturated fatty acid"

    # Check chain length (number of carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18 or c_count > 36:
        return False, f"Chain length {c_count} is outside the typical range for essential fatty acids (18-36 carbons)"

    # Check for at least 2 double bonds (polyunsaturated)
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bond_count < 2:
        return False, f"Only {double_bond_count} double bonds found, need at least 2 for polyunsaturated fatty acid"

    # Check that double bonds are in the cis configuration
    # RDKit doesn't directly support cis/trans checking, but we can infer based on the SMILES string
    # If the SMILES contains "/" or "\", it indicates cis/trans configuration
    if "/" not in smiles and "\\" not in smiles:
        return False, "Double bonds are not in cis configuration"

    return True, "Polyunsaturated fatty acid with multiple cis double bonds, carboxyl group, and appropriate chain length"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:59554',
                          'name': 'essential fatty acid',
                          'definition': 'Any member of the sub-set of polyunsaturated fatty acid for which there is an absolute dietary requirement.',
                          'parents': ['CHEBI:25681', 'CHEBI:26607']},
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