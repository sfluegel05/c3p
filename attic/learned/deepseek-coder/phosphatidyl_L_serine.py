"""
Classifies: CHEBI:18303 phosphatidyl-L-serine
"""
"""
Classifies: CHEBI:18303 phosphatidyl-L-serine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a phosphatidyl-L-serine based on its SMILES string.
    A phosphatidyl-L-serine has a glycerol backbone with two fatty acid chains,
    a phosphate group, and a serine moiety attached to the phosphate.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the glycerol backbone pattern (C-C-C with two oxygens attached to fatty acids)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for two ester groups (-O-C(=O)-) attached to the glycerol backbone
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2 for fatty acids"

    # Look for the phosphate group (P with oxygens attached)
    phosphate_pattern = Chem.MolFromSmarts("[PX4](=[OX1])([OX2])[OX2]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Look for the serine moiety (C-CO-NH2 with a carboxyl group)
    serine_pattern = Chem.MolFromSmarts("[CX4][CX4]([NX3])([CX3](=[OX1])[OX2H])")
    if not mol.HasSubstructMatch(serine_pattern):
        return False, "No serine moiety found"

    # Check if the serine is attached to the phosphate group
    serine_phosphate_pattern = Chem.MolFromSmarts("[PX4](=[OX1])([OX2])[OX2][CX4][CX4]([NX3])([CX3](=[OX1])[OX2H])")
    if not mol.HasSubstructMatch(serine_phosphate_pattern):
        return False, "Serine moiety is not attached to the phosphate group"

    # Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short to be fatty acids"

    # Check molecular weight - phosphatidyl-L-serine typically >600 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, "Molecular weight too low for phosphatidyl-L-serine"

    return True, "Contains glycerol backbone with two fatty acid chains, a phosphate group, and a serine moiety attached to the phosphate"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:18303',
                          'name': 'phosphatidyl-L-serine',
                          'definition': 'A class of aminophospholipids in which a phosphatidyl group is esterified to the hydroxy group of serine.',
                          'parents': ['CHEBI:18303', 'CHEBI:47778']},
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