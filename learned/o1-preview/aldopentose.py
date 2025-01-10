"""
Classifies: CHEBI:33916 aldopentose
"""
"""
Classifies: CHEBI:28053 aldopentose
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_aldopentose(smiles: str):
    """
    Determines if a molecule is an aldopentose based on its SMILES string.
    An aldopentose is a monosaccharide with five carbon atoms and an aldehyde group at one end,
    which can exist in open-chain or cyclic forms (furanose or pyranose rings).
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is an aldopentose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count number of carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons != 5:
        return False, f"Contains {num_carbons} carbon atoms, should be 5 for an aldopentose"
    
    # Check for open-chain form with aldehyde group at one end
    aldehyde_chain_pattern = Chem.MolFromSmarts("[O]=[C!-;!$(*=,#[!#6])]C([O,H])([O,H])[C@H]([O,H])[C@H]([O,H])CO")
    if mol.HasSubstructMatch(aldehyde_chain_pattern):
        return True, "Contains aldehyde group at end of five-carbon chain (open-chain aldopentose)"
    
    # Check for cyclic forms (furanose and pyranose rings)
    # Furanose ring pattern (5-membered ring with oxygen, derived from aldopentose)
    furanose_pattern = Chem.MolFromSmarts("C1OC([C@@H]([O,H])[C@H]([O,H])CO)C1[O,H]")
    if mol.HasSubstructMatch(furanose_pattern):
        return True, "Contains furanose ring derived from aldopentose"
    
    # Pyranose ring pattern (6-membered ring with oxygen, derived from aldopentose)
    pyranose_pattern = Chem.MolFromSmarts("C1OC([C@@H]([O,H])[C@H]([O,H])CO)C([O,H])C1[O,H]")
    if mol.HasSubstructMatch(pyranose_pattern):
        return True, "Contains pyranose ring derived from aldopentose"
    
    # Alternative patterns for cyclic forms (less strict stereochemistry)
    # General furanose ring from aldopentose
    general_furanose = Chem.MolFromSmarts("C1OC(C(O)C(O)CO)C1O")
    if mol.HasSubstructMatch(general_furanose):
        return True, "Contains furanose ring indicative of aldopentose"
    
    # General pyranose ring from aldopentose
    general_pyranose = Chem.MolFromSmarts("C1OC(C(O)C(O)C(O)CO)C1O")
    if mol.HasSubstructMatch(general_pyranose):
        return True, "Contains pyranose ring indicative of aldopentose"
    
    # If none of the patterns match, check for aldehyde at terminal carbon
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[CH2][CH]([O,H])[CH]([O,H])CO")
    if mol.HasSubstructMatch(aldehyde_pattern):
        return True, "Contains terminal aldehyde on five-carbon chain"
    
    # If still no match, molecule is not an aldopentose
    return False, "Does not match aldopentose structures (must have five carbons and aldehyde group at end or form specific cyclic structures)"
    

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:28053',
                              'name': 'aldopentose',
                              'definition': 'A pentose with a (potential) aldehyde group at one end.',
                              'parents': ['CHEBI:16842', 'CHEBI:24869']},
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
        'attempt': 1,
        'success': True,
        'best': True,
        'error': '',
        'stdout': None}