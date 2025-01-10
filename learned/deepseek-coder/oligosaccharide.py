"""
Classifies: CHEBI:50699 oligosaccharide
"""
"""
Classifies: CHEBI:50699 oligosaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_oligosaccharide(smiles: str):
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string.
    An oligosaccharide is a compound in which monosaccharide units are joined by glycosidic linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oligosaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycosidic linkage pattern (C-O-C between two sugar units)
    glycosidic_pattern = Chem.MolFromSmarts("[C][O][C]")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic linkage found"

    # Look for multiple hydroxyl groups (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 3:
        return False, f"Found {len(hydroxyl_matches)} hydroxyl groups, need at least 3"

    # Check for multiple sugar rings (pyranose or furanose)
    sugar_ring_pattern = Chem.MolFromSmarts("[C1][C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H]1O")
    sugar_ring_matches = mol.GetSubstructMatches(sugar_ring_pattern)
    if len(sugar_ring_matches) < 2:
        return False, f"Found {len(sugar_ring_matches)} sugar rings, need at least 2"

    # Check molecular weight - oligosaccharides typically >200 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for oligosaccharide"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 6:
        return False, "Too few carbons for oligosaccharide"
    if o_count < 5:
        return False, "Too few oxygens for oligosaccharide"

    return True, "Contains multiple sugar units joined by glycosidic linkages"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:50699',
                          'name': 'oligosaccharide',
                          'definition': 'A compound in which monosaccharide units are joined by glycosidic linkages. The term is commonly used to refer to a defined structure as opposed to a polymer of unspecified length or a homologous mixture. When the linkages are of other types the compounds are regarded as oligosaccharide analogues.',
                          'parents': ['CHEBI:16646', 'CHEBI:47778']},
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