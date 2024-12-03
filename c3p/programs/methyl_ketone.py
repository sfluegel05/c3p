"""
Classifies: CHEBI:51867 methyl ketone
"""
from rdkit import Chem

def is_methyl_ketone(smiles: str):
    """
    Determines if a molecule is a methyl ketone (RC(=O)CH3 where R is not H).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl ketone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a methyl ketone
    methyl_ketone_pattern = Chem.MolFromSmarts("C(=O)C")

    if mol.HasSubstructMatch(methyl_ketone_pattern):
        for match in mol.GetSubstructMatches(methyl_ketone_pattern):
            carbonyl_carbon = mol.GetAtomWithIdx(match[0])
            methyl_carbon = mol.GetAtomWithIdx(match[1])

            # Check if the carbonyl carbon is bonded to another carbon (R is not H)
            for neighbor in carbonyl_carbon.GetNeighbors():
                if neighbor.GetIdx() != methyl_carbon.GetIdx() and neighbor.GetSymbol() == 'C':
                    return True, "Molecule contains a methyl ketone group"
        return False, "Molecule has a C(=O)C group but R is H or non-carbon"
    else:
        return False, "Molecule does not contain a methyl ketone group"

# Example usage:
# print(is_methyl_ketone("CCC(C)=O"))  # Should return (True, "Molecule contains a methyl ketone group")
# print(is_methyl_ketone("CC(=O)C"))  # Should return (True, "Molecule contains a methyl ketone group")
# print(is_methyl_ketone("CCO"))  # Should return (False, "Molecule does not contain a methyl ketone group")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:51867',
                          'name': 'methyl ketone',
                          'definition': 'A ketone of formula RC(=O)CH3 (R =/= '
                                        'H).',
                          'parents': ['CHEBI:17087']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 27,
    'num_false_positives': 13,
    'num_true_negatives': 7,
    'num_false_negatives': 7,
    'precision': 0.675,
    'recall': 0.7941176470588235,
    'f1': 0.7297297297297296,
    'accuracy': None}