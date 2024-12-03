"""
Classifies: CHEBI:167095 secondary fatty alcohol
"""
from rdkit import Chem

def is_secondary_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a secondary fatty alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary fatty alcohol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a hydroxyl group
    hydroxyl_groups = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'O' and atom.GetDegree() == 1]
    if not hydroxyl_groups:
        return False, "No hydroxyl group found"

    # Check if the hydroxyl group is attached to a non-terminal carbon
    for hydroxyl in hydroxyl_groups:
        attached_carbons = [neighbor for neighbor in hydroxyl.GetNeighbors() if neighbor.GetSymbol() == 'C']
        if len(attached_carbons) != 1:
            continue
        attached_carbon = attached_carbons[0]
        
        # Check if the attached carbon is non-terminal
        carbon_neighbors = [neighbor for neighbor in attached_carbon.GetNeighbors() if neighbor.GetSymbol() == 'C']
        if len(carbon_neighbors) < 2:
            continue
        
        # Check the length of the carbon chain
        carbon_chain_length = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
        if carbon_chain_length < 3:
            return False, f"Carbon chain length {carbon_chain_length} is less than 3"

        return True, "Valid secondary fatty alcohol"

    return False, "Hydroxyl group not attached to a non-terminal carbon"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:167095',
                          'name': 'secondary fatty alcohol',
                          'definition': 'A fatty alcohol consisting of a chain '
                                        'of 3 to greater than 27 carbon atoms '
                                        'in which a hydroxy group is attached '
                                        'to a saturated carbon atom different '
                                        'from the terminal carbons. Secondary '
                                        'fatty alcohols may be saturated or '
                                        'unsaturated and may be branched or '
                                        'unbranched.',
                          'parents': ['CHEBI:24026', 'CHEBI:35681']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 13,
    'num_false_positives': 4,
    'num_true_negatives': 9,
    'num_false_negatives': 0,
    'precision': 0.7647058823529411,
    'recall': 1.0,
    'f1': 0.8666666666666666,
    'accuracy': None}