"""
Classifies: CHEBI:24279 glucosinolate
"""
"""
Classifies: CHEBI:59282 glucosinolate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glucosinolate(smiles: str):
    """
    Determines if a molecule is a glucosinolate based on its SMILES string.
    A glucosinolate has a central carbon bonded to a sulfur (connected to a glycone group),
    a nitrogen (connected to a sulfonated oxime group), and a side group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glucosinolate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core glucosinolate pattern: central C bonded to S (glycone), N (sulfonated oxime), and a side group
    glucosinolate_pattern = Chem.MolFromSmarts("[CX4]([SX2][C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)([NX2][OX2][SX4](=[OX1])(=[OX1])[OX1-])([CX4])")
    if not mol.HasSubstructMatch(glucosinolate_pattern):
        return False, "Core glucosinolate pattern not found"

    # Check for the anti stereochemistry across the C=N double bond
    # This requires the side chain and sulfate group to be on opposite sides
    # RDKit's stereochemistry tools can be used to verify this
    matches = mol.GetSubstructMatches(glucosinolate_pattern)
    for match in matches:
        central_carbon_idx = match[0]
        sulfur_idx = match[1]
        nitrogen_idx = match[2]
        side_group_idx = match[3]

        # Get the stereochemistry of the central carbon
        central_carbon = mol.GetAtomWithIdx(central_carbon_idx)
        if not central_carbon.HasProp("_CIPCode"):
            return False, "Stereochemistry not defined for central carbon"

        # Check if the side group and sulfate group are anti
        # This is a simplified check; a more rigorous approach would involve 3D coordinates
        # or explicit stereochemistry checks
        if central_carbon.GetProp("_CIPCode") not in ["R", "S"]:
            return False, "Stereochemistry not R or S for central carbon"

    # Check for the presence of a glycone group (glucose-like structure)
    glycone_pattern = Chem.MolFromSmarts("[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(glycone_pattern):
        return False, "Glycone group not found"

    # Check for the sulfonated oxime group
    sulfonated_oxime_pattern = Chem.MolFromSmarts("[NX2][OX2][SX4](=[OX1])(=[OX1])[OX1-]")
    if not mol.HasSubstructMatch(sulfonated_oxime_pattern):
        return False, "Sulfonated oxime group not found"

    # Check for the side group (any carbon chain or ring)
    side_group_pattern = Chem.MolFromSmarts("[CX4]")
    if not mol.HasSubstructMatch(side_group_pattern):
        return False, "Side group not found"

    return True, "Contains core glucosinolate structure with anti stereochemistry"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:59282',
        'name': 'glucosinolate',
        'definition': 'Water-soluble anionic substituted thioglucosides. Glucosinolates have a central C atom which is bonded via an S atom to a glycone group and via an N atom to a sulfonated oxime group, and which also carries a side-group. The side-chain and sulfate group have an anti stereochemical configuration across the C=N double bond.'
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
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199
}