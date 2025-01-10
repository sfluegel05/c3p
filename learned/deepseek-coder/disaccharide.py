"""
Classifies: CHEBI:36233 disaccharide
"""
"""
Classifies: CHEBI:36218 disaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide is a compound in which two monosaccharides are joined by a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a disaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of two sugar rings (monosaccharides)
    # A sugar ring typically has 5 or 6 atoms, with multiple hydroxyl groups
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 2:
        return False, "Less than 2 rings found, need exactly 2 for a disaccharide"

    # Check for glycosidic bond (oxygen linking two rings)
    glycosidic_bond_pattern = Chem.MolFromSmarts("[O;X2][C;X4][C;X4][O;X2]")
    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No glycosidic bond found between rings"

    # Count hydroxyl groups (OH) to ensure it's a sugar
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1)
    if hydroxyl_count < 4:
        return False, "Too few hydroxyl groups for a disaccharide"

    # Check molecular weight - disaccharides typically have a molecular weight between 300 and 500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 500:
        return False, "Molecular weight outside typical range for disaccharides"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 10 or c_count > 20:
        return False, "Carbon count outside typical range for disaccharides"
    if o_count < 8 or o_count > 12:
        return False, "Oxygen count outside typical range for disaccharides"

    return True, "Contains two monosaccharide units joined by a glycosidic bond"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36218',
                          'name': 'disaccharide',
                          'definition': 'A compound in which two monosaccharides are joined by a glycosidic bond.',
                          'parents': ['CHEBI:36218', 'CHEBI:47778']},
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