"""
Classifies: CHEBI:23053 catechin
"""
"""
Classifies: CHEBI:23053 catechin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.
    A catechin is a member of the class of hydroxyflavan that has a flavan-3-ol skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catechin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more general flavan-3-ol core structure pattern
    # This pattern matches the 2-phenylchromane skeleton with at least one hydroxyl group
    flavan_core = Chem.MolFromSmarts("[C@H]1[C@H](O)Cc2c(O)cc(O)cc2O1")
    if not mol.HasSubstructMatch(flavan_core):
        # Try a more relaxed pattern without stereochemistry
        flavan_core_relaxed = Chem.MolFromSmarts("C1C(O)Cc2c(O)cc(O)cc2O1")
        if not mol.HasSubstructMatch(flavan_core_relaxed):
            return False, "No flavan-3-ol core structure found"

    # Check for characteristic catechin features
    # 1. At least one aromatic ring with hydroxyl groups
    aromatic_oh_pattern = Chem.MolFromSmarts("c[OH]")
    if not mol.HasSubstructMatch(aromatic_oh_pattern):
        return False, "No aromatic hydroxyl groups found"

    # 2. Check molecular weight range (150-1000 Da typical for catechin derivatives)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150 or mol_wt > 1000:
        return False, f"Molecular weight {mol_wt:.1f} outside typical catechin range"

    # 3. Check for common catechin substitutions
    # Gallate ester pattern
    gallate_pattern = Chem.MolFromSmarts("[OH]c1cc(O)c(O)c(O)c1")
    # Sulfate pattern
    sulfate_pattern = Chem.MolFromSmarts("[OH]S(=O)(=O)[OH]")
    # Common substitutions
    substitution_patterns = [gallate_pattern, sulfate_pattern]
    
    has_substitution = any(mol.HasSubstructMatch(patt) for patt in substitution_patterns)
    
    # If no common substitutions, check for basic catechin structure
    if not has_substitution:
        # Basic catechin should have at least 3 hydroxyl groups
        hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
        hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
        if len(hydroxyl_matches) < 3:
            return False, f"Found {len(hydroxyl_matches)} hydroxyl groups, need at least 3 for basic catechin"

    return True, "Contains flavan-3-ol core structure with characteristic features"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23053',
                          'name': 'catechin',
                          'definition': 'Members of the class of hydroxyflavan that have a flavan-3-ol skeleton and its substituted derivatives.',
                          'parents': ['CHEBI:23053', 'CHEBI:23053']},
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