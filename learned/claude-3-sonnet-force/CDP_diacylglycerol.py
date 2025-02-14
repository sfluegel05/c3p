"""
Classifies: CHEBI:17962 CDP-diacylglycerol
"""
"""
Classifies: CHEBI:38836 CDP-diacylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_CDP_diacylglycerol(smiles: str):
    """
    Determines if a molecule is a CDP-diacylglycerol based on its SMILES string.
    A CDP-diacylglycerol is a CDP-glycerol having unspecified acyl groups
    (most commonly fatty acyl groups) at the 1- and 2-positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a CDP-diacylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for CDP-glycerol backbone pattern
    cdp_glycerol_pattern = Chem.MolFromSmarts("[C@H]1[C@@H]([C@H]([C@@H](O1)N1C=NC(=NC1=O)N)O)OP(=O)(O)OP(=O)(O)O")
    if not mol.HasSubstructMatch(cdp_glycerol_pattern):
        return False, "No CDP-glycerol backbone found"
        
    # Look for at least 2 ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"

    # Look for fatty acid chains (long aliphatic chains)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"

    # Optionally, check for long chains based on molecular weight or rotatable bond count
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if mol_wt < 500 or n_rotatable < 10:
        return False, "Chains too short to be fatty acids"

    return True, "Contains CDP-glycerol backbone with at least 2 fatty acid chains attached via ester bonds"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:38836',
        'name': 'CDP-diacylglycerol',
        'definition': 'A CDP-glycerol having unspecified acyl groups (most commonly fatty acyl groups) at the 1- and 2-positions.',
        'parents': ['CHEBI:38075', 'CHEBI:36342']
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    },
    'message': None,
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 450,
    'num_false_positives': 3,
    'num_true_negatives': 182398,
    'num_false_negatives': 39,
    'num_negatives': None,
    'precision': 0.9934640522875817,
    'recall': 0.920081967213115,
    'f1': 0.9560439560439561,
    'accuracy': 0.9997879011140626
}