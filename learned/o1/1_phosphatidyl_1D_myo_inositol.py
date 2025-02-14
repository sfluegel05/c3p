"""
Classifies: CHEBI:16749 1-phosphatidyl-1D-myo-inositol
"""
"""
Classifies: CHEBI:28874 1-phosphatidyl-1D-myo-inositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_phosphatidyl_1D_myo_inositol(smiles: str):
    """
    Determines if a molecule is a 1-phosphatidyl-1D-myo-inositol based on its SMILES string.
    This involves checking for:
    - A glycerol backbone esterified at positions 1 and 2 with fatty acids
    - A phosphate group at position 3 of glycerol
    - An inositol ring attached to the phosphate group at position 1 of the inositol
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
    
    # Check for glycerol backbone esterified at positions 1 and 2 with fatty acids
    glycerol_pattern = Chem.MolFromSmarts("[C@@H]([CH2]OC(=O)*)([CH2]OC(=O)*)")  # Glycerol with esters at positions 1 and 2
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone esterified at positions 1 and 2 found"
    
    # Check for phosphate group at position 3 of glycerol
    phosphate_pattern = Chem.MolFromSmarts("COP(=O)(O)O")  # Phosphate group attached via oxygen
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group attached to glycerol found"
    
    # Check for inositol ring attached to phosphate group at position 1
    # Inositol ring: six-membered ring with six hydroxyl groups
    inositol_ring = Chem.MolFromSmarts("C1(CO)C(O)C(O)C(O)C(O)O1")
    if not mol.HasSubstructMatch(inositol_ring):
        return False, "No inositol ring attached to phosphate group found"
    
    # Check for correct stereochemistry in inositol ring (1D-myo-inositol)
    # 1D-myo-inositol has specific stereochemistry at the ring carbons
    stereo_inositol = Chem.MolFromSmiles("O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O")
    stereo_match = mol.HasSubstructMatch(stereo_inositol, useChirality=True)
    if not stereo_match:
        return False, "Inositol ring does not have correct stereochemistry (1D-myo-inositol)"
    
    # Confirm attachment point of inositol ring to phosphate group is at position 1
    phosphate_inositol_pattern = Chem.MolFromSmarts("OP(=O)(O)OC1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)CO1")  # Phosphate connected to O1 of inositol
    if not mol.HasSubstructMatch(phosphate_inositol_pattern):
        return False, "Phosphate group not attached at position 1 of inositol ring"
    
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
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None
}