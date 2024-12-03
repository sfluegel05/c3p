"""
Classifies: CHEBI:64574 lysophosphatidylethanolamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_lysophosphatidylethanolamine(smiles: str):
    """
    Determines if a molecule is a lysophosphatidylethanolamine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lysophosphatidylethanolamine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the substructure patterns for lysophosphatidylethanolamine
    glycerol_pattern = Chem.MolFromSmarts("[C@H](CO)CO")
    phosphate_pattern = Chem.MolFromSmarts("COP(=O)(O)OCCN")
    fatty_acid_pattern = Chem.MolFromSmarts("C(=O)C")
    fatty_acid_pattern_2 = Chem.MolFromSmarts("C(=O)")

    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone not found"

    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate group with ethanolamine not found"

    # Check for the presence of only one fatty acid chain
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    fatty_acid_matches_2 = mol.GetSubstructMatches(fatty_acid_pattern_2)
    
    if len(fatty_acid_matches) != 1 and len(fatty_acid_matches_2) != 1:
        return False, "Number of fatty acid chains is not equal to one"

    return True, "Lysophosphatidylethanolamine structure found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:64574',
                          'name': 'lysophosphatidylethanolamine',
                          'definition': 'A glycerophosphoethanolamine '
                                        'resulting from partial hydrolysis of '
                                        'a phosphatidylethanolamine, which '
                                        'removes one of the fatty acid groups. '
                                        'The structure is depicted in the '
                                        'image where R(1) = acyl, R(2) = H or '
                                        'where R(1) = H, R(2) = acyl.',
                          'parents': ['CHEBI:36314']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 10,
    'num_false_positives': 0,
    'num_true_negatives': 14,
    'num_false_negatives': 4,
    'precision': 1.0,
    'recall': 0.7142857142857143,
    'f1': 0.8333333333333333,
    'accuracy': None}