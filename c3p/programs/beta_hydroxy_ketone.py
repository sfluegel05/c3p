"""
Classifies: CHEBI:55380 beta-hydroxy ketone
"""
from rdkit import Chem

def is_beta_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a beta-hydroxy ketone (a ketone containing a hydroxy group on the beta-carbon relative to the C=O group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all ketone groups (C=O)
    ketone_pattern = Chem.MolFromSmarts("C(=O)[C,c]")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    if not ketone_matches:
        return False, "No ketone groups found"

    for match in ketone_matches:
        carbonyl_carbon = match[0]

        # Get the alpha carbons (directly attached to the carbonyl carbon)
        alpha_carbons = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(carbonyl_carbon).GetNeighbors() if nbr.GetIdx() != match[1]]

        # Get the beta carbons (attached to alpha carbons)
        beta_carbons = []
        for alpha_carbon in alpha_carbons:
            beta_carbons.extend([nbr.GetIdx() for nbr in mol.GetAtomWithIdx(alpha_carbon).GetNeighbors() if nbr.GetIdx() != carbonyl_carbon])

        # Check for hydroxy group on beta carbons
        for beta_carbon in beta_carbons:
            for nbr in mol.GetAtomWithIdx(beta_carbon).GetNeighbors():
                if nbr.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(beta_carbon, nbr.GetIdx()).GetBondTypeAsDouble() == 1.0:
                    if mol.GetAtomWithIdx(nbr.GetIdx()).GetDegree() == 1:  # Ensure it's a hydroxyl group (OH)
                        return True, "Beta-hydroxy ketone found"

    return False, "No beta-hydroxy ketone structure found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:55380',
                          'name': 'beta-hydroxy ketone',
                          'definition': 'A ketone containing a hydroxy group '
                                        'on the beta-carbon relative to the '
                                        'C=O group.',
                          'parents': ['CHEBI:17087']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 4,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 11,
    'num_false_positives': 1,
    'num_true_negatives': 11,
    'num_false_negatives': 1,
    'precision': 0.9166666666666666,
    'recall': 0.9166666666666666,
    'f1': 0.9166666666666666,
    'accuracy': None}