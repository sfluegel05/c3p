"""
Classifies: CHEBI:18179 phosphoinositide
"""
"""
Classifies: CHEBI:18179 phosphoinositide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    A phosphoinositide is a phosphatidylinositol that is phosphorylated at one or more of the hydroxy groups of inositol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphoinositide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 3 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
        
    # Look for 2 ester groups (-O-C(=O)-) for fatty acid chains
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2 for fatty acid chains"

    # Look for phosphate group attached to glycerol (P-O-C)
    phosphate_glycerol_pattern = Chem.MolFromSmarts("[PX4](=[OX1])([OX2][CX4])[OX2]")
    if not mol.HasSubstructMatch(phosphate_glycerol_pattern):
        return False, "No phosphate group attached to glycerol backbone"

    # Look for inositol ring (6-membered ring with 6 oxygens)
    inositol_pattern = Chem.MolFromSmarts("[C@H]1([C@H]([C@H]([C@H]([C@H]([C@H]1O)O)O)O)O)O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"

    # Look for at least one phosphate group attached to inositol (P-O-C in inositol)
    phosphate_inositol_pattern = Chem.MolFromSmarts("[PX4](=[OX1])([OX2][C@H]1[C@H]([C@H]([C@H]([C@H]([C@H]1O)O)O)O)O)[OX2]")
    if not mol.HasSubstructMatch(phosphate_inositol_pattern):
        return False, "No phosphate group attached to inositol ring"

    # Count the number of phosphate groups attached to inositol
    phosphate_inositol_matches = mol.GetSubstructMatches(phosphate_inositol_pattern)
    if len(phosphate_inositol_matches) < 1:
        return False, "No phosphorylation on inositol ring"

    return True, "Contains glycerol backbone with 2 fatty acid chains, a phosphate group, and phosphorylated inositol ring"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:18179',
                          'name': 'phosphoinositide',
                          'definition': 'Any phosphatidylinositol that is '
                                        'phosphorylated at one or more of the '
                                        'hydroxy groups of inositol.',
                          'parents': ['CHEBI:26039', 'CHEBI:26040']},
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