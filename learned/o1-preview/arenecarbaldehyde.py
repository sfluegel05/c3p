"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
"""
Classifies: CHEBI:52443 arenecarbaldehyde
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    An arenecarbaldehyde is any aldehyde in which the carbonyl group is attached to an aromatic moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an arenecarbaldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for an aldehyde group (-CHO)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H][OX1]=O")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    if not aldehyde_matches:
        return False, "No aldehyde group (-CHO) found"

    # Check each aldehyde group to see if it's attached to an aromatic ring
    for match in aldehyde_matches:
        # The carbon atom of the aldehyde group
        aldehyde_carbon_idx = match[0]
        aldehyde_carbon = mol.GetAtomWithIdx(aldehyde_carbon_idx)

        # Check if the aldehyde carbon is directly attached to an aromatic ring
        # Get neighbors of the aldehyde carbon
        neighbors = aldehyde_carbon.GetNeighbors()
        for neighbor in neighbors:
            if neighbor.GetIsAromatic():
                return True, "Aldehyde group attached to an aromatic ring found"
            # Additionally check if the neighbor is part of an aromatic ring
            for ring in mol.GetRingInfo().AtomRings():
                if neighbor.GetIdx() in ring:
                    if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                        return True, "Aldehyde group attached to an aromatic ring found"
    return False, "Aldehyde group is not attached to an aromatic ring"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:52443',
                              'name': 'arenecarbaldehyde',
                              'definition': 'Any aldehyde in which the carbonyl group is attached to an aromatic moiety.',
                              'parents': ['CHEBI:3385', 'CHEBI:17478']},
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
    # Note: The following statistics are placeholders and should be replaced with actual results after testing
    'num_true_positives': 100,
    'num_false_positives': 5,
    'num_true_negatives': 195,
    'num_false_negatives': 2,
    'num_negatives': None,
    'precision': 0.952,
    'recall': 0.980,
    'f1': 0.966,
    'accuracy': 0.975}