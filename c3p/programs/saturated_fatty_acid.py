"""
Classifies: CHEBI:26607 saturated fatty acid
"""
from rdkit import Chem

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saturated fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (COOH)
    carboxylic_acid = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetTotalDegree() == 3:
            neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
            if neighbors.count('O') == 2 and neighbors.count('C') == 1:
                carboxylic_acid = True
                break

    if not carboxylic_acid:
        return False, "No carboxylic acid group (COOH) found"

    # Check for absence of carbon to carbon multiple bonds
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.rdchem.BondType.SINGLE and bond.GetBeginAtom().GetSymbol() == 'C' and bond.GetEndAtom().GetSymbol() == 'C':
            return False, "Contains carbon to carbon multiple bonds"

    return True, "Saturated fatty acid"

# Example usage:
# print(is_saturated_fatty_acid('CCC(C)CCCCCCCCCCCCC(O)=O'))  # Should return (True, "Saturated fatty acid")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26607',
                          'name': 'saturated fatty acid',
                          'definition': 'Any fatty acid containing no carbon '
                                        'to carbon multiple bonds. Known to '
                                        'produce adverse biological effects '
                                        'when ingested to excess.',
                          'parents': ['CHEBI:35366']},
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
    'num_true_positives': 10,
    'num_false_positives': 2,
    'num_true_negatives': 8,
    'num_false_negatives': 0,
    'precision': 0.8333333333333334,
    'recall': 1.0,
    'f1': 0.9090909090909091,
    'accuracy': None}