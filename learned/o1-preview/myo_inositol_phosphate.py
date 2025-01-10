"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
"""
Classifies: CHEBI:12348 myo-inositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    A myo-inositol phosphate is an inositol phosphate in which the inositol component has myo-configuration.
    This involves a cyclohexane ring with six hydroxyl or phosphate ester groups in the myo configuration.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a myo-inositol phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for myo-inositol phosphate
    # Cyclohexane ring with specific stereochemistry and oxygen substituents
    pattern_smarts = """
    [C@H]1([O]) [C@@H]([O]) [C@H]([O]) [C@@H]([O]) [C@H]([O]) [C@@H]1([O])
    """
    # Clean up the SMARTS pattern
    pattern_smarts = ''.join(pattern_smarts.split())

    # Create the pattern molecule
    pattern = Chem.MolFromSmarts(pattern_smarts)
    if pattern is None:
        return False, "Invalid SMARTS pattern"

    # Check for substructure match with chirality
    matches = mol.GetSubstructMatches(pattern, useChirality=True)
    if not matches:
        return False, "Molecule does not match myo-inositol core with correct stereochemistry"

    # Check that each substituent is oxygen
    # and identify if at least one substituent is a phosphate group
    has_phosphate = False
    for match in matches:
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            # Get neighbors excluding ring bonds
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in match:
                    # Neighbor should be oxygen
                    if neighbor.GetAtomicNum() != 8:
                        return False, "Substituent is not oxygen"
                    # Check if oxygen is connected to phosphorus
                    is_phosphate = False
                    for nb in neighbor.GetNeighbors():
                        if nb.GetAtomicNum() == 15:
                            is_phosphate = True
                            has_phosphate = True
                            break
                    if not is_phosphate:
                        # Confirm it's a hydroxyl group (OH)
                        if neighbor.GetTotalNumHs() != 1:
                            return False, "Oxygen substituent is neither hydroxyl nor phosphate group"
    if not has_phosphate:
        return False, "No phosphate group found in molecule"

    return True, "Molecule matches myo-inositol phosphate pattern"

__metadata__ = {   
    'chemical_class': {   
        'id': 'CHEBI:12348',
        'name': 'myo-inositol phosphate',
        'definition': 'An inositol phosphate in which the inositol component has myo-configuration.',
        'parents': ['CHEBI:24848', 'CHEBI:24845']
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
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': None,
    'num_false_positives': None,
    'num_true_negatives': None,
    'num_false_negatives': None,
    'num_negatives': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None
}