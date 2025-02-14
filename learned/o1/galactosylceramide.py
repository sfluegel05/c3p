"""
Classifies: CHEBI:36498 galactosylceramide
"""
"""
Classifies: CHEBI:34219 galactosylceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide based on its SMILES string.
    A galactosylceramide is a ceramide with a galactose monosaccharide head group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a galactosylceramide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ceramide backbone: sphingosine base linked via amide bond to fatty acid
    # Sphingosine backbone pattern: long-chain amino alcohol with trans double bond
    sphingosine_pattern = Chem.MolFromSmarts("C(O)[C@@H](NC(=O))C=C")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine base found"

    # Look for amide bond to fatty acid
    amide_pattern = Chem.MolFromSmarts("NC(=O)C")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond to fatty acid found"

    # Check for galactose head group attached via glycosidic bond
    # Galactose pattern (pyranose form)
    galactose_pattern = Chem.MolFromSmarts("CO[C@H]1O[C@H](O)[C@@H](O)[C@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(galactose_pattern):
        return False, "No galactose head group found"

    # Check for glycosidic bond between sphingosine and galactose
    glycosidic_bond_pattern = Chem.MolFromSmarts("C[NH]C(=O)[C@@H](CO[C@H]1O[C@H](O)[C@@H](O)[C@H](O)[C@H]1O)")
    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No glycosidic bond between sphingosine and galactose found"

    # Count carbon atoms to verify long-chain fatty acid
    c_chain = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(c_chain) < 20:
        return False, "Too few carbons for ceramide with long-chain fatty acid"

    # Verify the presence of hydroxyl groups
    oh_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalDegree() == 1)
    if oh_count < 4:
        return False, "Insufficient hydroxyl groups for galactosylceramide"

    return True, "Contains ceramide backbone with galactose head group attached via glycosidic bond"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:34219',
                              'name': 'galactosylceramide',
                              'definition': 'Any of the cerebrosides in which '
                                            'the monosaccharide head group is '
                                            'galactose.',
                              'parents': ['CHEBI:33563', 'CHEBI:75771']},
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