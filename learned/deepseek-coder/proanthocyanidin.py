"""
Classifies: CHEBI:26267 proanthocyanidin
"""
"""
Classifies: CHEBI:26154 proanthocyanidin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_proanthocyanidin(smiles: str):
    """
    Determines if a molecule is a proanthocyanidin based on its SMILES string.
    A proanthocyanidin is a flavonoid oligomer obtained by the condensation of two or more units of hydroxyflavans.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a proanthocyanidin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more flexible hydroxyflavan pattern
    hydroxyflavan_pattern = Chem.MolFromSmarts("[cX3]1[cX3][cX3][cX3][cX3][cX3]1[C@@H]2[C@H](O)[C@H](O)[C@H](O)[C@H](O2)")
    if not mol.HasSubstructMatch(hydroxyflavan_pattern):
        return False, "No hydroxyflavan unit found"

    # Count the number of hydroxyflavan units
    hydroxyflavan_matches = mol.GetSubstructMatches(hydroxyflavan_pattern)
    if len(hydroxyflavan_matches) < 2:
        return False, f"Found {len(hydroxyflavan_matches)} hydroxyflavan units, need at least 2"

    # Look for multiple aromatic rings (at least two)
    aromatic_rings = Chem.GetSSSR(mol)
    if aromatic_rings < 2:
        return False, f"Found {aromatic_rings} aromatic rings, need at least 2"

    # Check for multiple hydroxyl groups (at least 4)
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1)
    if hydroxyl_count < 4:
        return False, f"Found {hydroxyl_count} hydroxyl groups, need at least 4"

    # Check for ether linkages (C-O-C)
    ether_pattern = Chem.MolFromSmarts("[CX4][OX2][CX4]")
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    if len(ether_matches) < 1:
        return False, "No ether linkages found"

    # Check molecular weight - proanthocyanidins typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for proanthocyanidin"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:
        return False, "Too few carbons for proanthocyanidin"
    if o_count < 6:
        return False, "Too few oxygens for proanthocyanidin"

    return True, "Contains multiple hydroxyflavan units with aromatic rings, hydroxyl groups, and ether linkages"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26154',
                          'name': 'proanthocyanidin',
                          'definition': 'A flavonoid oligomer obtained by the condensation of two or more units of hydroxyflavans.',
                          'parents': ['CHEBI:26154', 'CHEBI:26154']},
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