"""
Classifies: CHEBI:33567 catecholamine
"""
"""
Classifies: CHEBI:33567 catecholamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    A catecholamine is defined as 4-(2-Aminoethyl)pyrocatechol [4-(2-aminoethyl)benzene-1,2-diol] 
    and derivatives formed by substitution.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catecholamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the core catechol structure (benzene-1,2-diol)
    catechol_pattern = Chem.MolFromSmarts("[c]1([OH])[c][c][c][c]1[OH]")
    if not mol.HasSubstructMatch(catechol_pattern):
        return False, "No catechol (benzene-1,2-diol) structure found"

    # Look for the 2-aminoethyl side chain (CCN)
    aminoethyl_pattern = Chem.MolFromSmarts("[CH2][CH2][NH2]")
    if not mol.HasSubstructMatch(aminoethyl_pattern):
        return False, "No 2-aminoethyl side chain found"

    # Check if the aminoethyl side chain is attached to the catechol ring
    combined_pattern = Chem.MolFromSmarts("[c]1([OH])[c][c]([CH2][CH2][NH2])[c][c]1[OH]")
    if not mol.HasSubstructMatch(combined_pattern):
        return False, "2-aminoethyl side chain not attached to catechol ring"

    # Check for additional substitutions on the catechol ring or side chain
    # This allows for derivatives formed by substitution
    # We can count the number of atoms to see if there are additional substituents
    num_atoms = mol.GetNumAtoms()
    if num_atoms > 14:  # Basic catecholamine has 14 atoms
        return True, "Catecholamine with additional substitutions"

    return True, "Basic catecholamine structure found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33567',
                          'name': 'catecholamine',
                          'definition': '4-(2-Aminoethyl)pyrocatechol [4-(2-aminoethyl)benzene-1,2-diol] and derivatives formed by substitution.',
                          'parents': ['CHEBI:33566', 'CHEBI:33568']},
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