"""
Classifies: CHEBI:23086 chalcones
"""
"""
Classifies: CHEBI:23003 chalcone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone based on its SMILES string.
    A chalcone is a ketone that is 1,3-diphenylpropenone (benzylideneacetophenone), ArCH=CH(=O)Ar, and its derivatives formed by substitution.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chalcone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the core chalcone structure: Ar-CH=CH-C(=O)-Ar
    chalcone_pattern = Chem.MolFromSmarts("[c;H1]:[c;H1]-[CH]=[CH]-[C](=O)-[c;H1]:[c;H1]")
    if not mol.HasSubstructMatch(chalcone_pattern):
        return False, "No chalcone core structure found (Ar-CH=CH-C(=O)-Ar)"

    # Check for the presence of two aromatic rings (Ar)
    aromatic_rings = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()]
    if len(aromatic_rings) < 2:
        return False, "Less than two aromatic rings found"

    # Check for the presence of the propenone chain (CH=CH-C=O)
    propenone_chain = Chem.MolFromSmarts("[CH]=[CH]-[C](=O)")
    if not mol.HasSubstructMatch(propenone_chain):
        return False, "No propenone chain (CH=CH-C=O) found"

    # Check for the presence of a ketone group (C=O)
    ketone_group = Chem.MolFromSmarts("[C](=O)")
    if not mol.HasSubstructMatch(ketone_group):
        return False, "No ketone group (C=O) found"

    # Check for the presence of at least one double bond in the propenone chain
    double_bond = Chem.MolFromSmarts("[CH]=[CH]")
    if not mol.HasSubstructMatch(double_bond):
        return False, "No double bond in the propenone chain"

    return True, "Contains the core chalcone structure (Ar-CH=CH-C(=O)-Ar) with possible substitutions on the aromatic rings"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23003',
                          'name': 'chalcone',
                          'definition': 'A ketone that is 1,3-diphenylpropenone (benzylideneacetophenone), ArCH=CH(=O)Ar, and its derivatives formed by substitution.',
                          'parents': ['CHEBI:23003', 'CHEBI:23003']},
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