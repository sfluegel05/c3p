"""
Classifies: CHEBI:16749 1-phosphatidyl-1D-myo-inositol
"""
"""
Classifies: CHEBI:28874 1-phosphatidyl-1D-myo-inositol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_1_phosphatidyl_1D_myo_inositol(smiles: str):
    """
    Determines if a molecule is a 1-phosphatidyl-1D-myo-inositol based on its SMILES string.
    This involves checking for:
    - A glycerol backbone esterified at positions sn-1 and sn-2 with fatty acids
    - A phosphate group at position sn-3 of glycerol
    - An inositol ring attached to the phosphate group at position 1
    - The inositol ring should be the 1D-myo isomer (specific stereochemistry)
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a 1-phosphatidyl-1D-myo-inositol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for glycerol backbone esterified at sn-1 and sn-2 with fatty acids
    # Glycerol backbone with esters at sn-1 and sn-2 positions
    glycerol_pattern = Chem.MolFromSmarts('[C@@H](COC(=O)[#6])[C@H](COC(=O)[#6])[CH2O]')
    if glycerol_pattern is None:
        return False, "Invalid glycerol pattern"
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone esterified at sn-1 and sn-2 positions found"
    
    # Check for phosphate group at sn-3 position of glycerol
    phosphate_pattern = Chem.MolFromSmarts('O[P](=O)(O)O')
    if phosphate_pattern is None:
        return False, "Invalid phosphate pattern"
    # Attach phosphate group to sn-3 position
    glycerol_phosphate_pattern = Chem.MolFromSmarts('[CH2O][P](=O)(O)O')
    if glycerol_phosphate_pattern is None:
        return False, "Invalid glycerol phosphate pattern"
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "No phosphate group at sn-3 position of glycerol found"
    
    # Check for inositol ring attached to phosphate group
    # Inositol ring: six-membered ring with six hydroxyl groups
    # Connect phosphate to inositol at position 1
    inositol_pattern = Chem.MolFromSmarts('O[P](=O)(O)OC1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O')
    if inositol_pattern is None:
        return False, "Invalid inositol pattern"
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring attached to phosphate group found"
    
    # Check for correct stereochemistry of 1D-myo-inositol
    # Stereochemistry of inositol ring
    stereo_inositol_pattern = Chem.MolFromSmarts('C1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O')
    if stereo_inositol_pattern is None:
        return False, "Invalid stereochemistry pattern"
    matches = mol.GetSubstructMatches(stereo_inositol_pattern, useChirality=True)
    if not matches:
        return False, "Inositol ring does not have correct stereochemistry (1D-myo-inositol)"
    
    return True, "Molecule is a 1-phosphatidyl-1D-myo-inositol with correct structure"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:28874',
        'name': '1-phosphatidyl-1D-myo-inositol',
        'definition': 'A phosphatidylinositol in which the inositol moiety is the 1D-myo isomer and the phosphatidyl group is located at its position 1.',
        'parents': ['CHEBI:49117', 'CHEBI:24851']
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
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None
}