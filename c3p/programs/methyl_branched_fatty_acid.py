"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
from rdkit import Chem

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl-branched fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl-branched fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the carboxylic acid group (COOH)
    carboxylic_acid = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetTotalDegree() == 3:
            neighbors = [neighbor.GetSymbol() for neighbor in atom.GetNeighbors()]
            if 'O' in neighbors and neighbors.count('O') == 2:
                carboxylic_acid = True
                break

    if not carboxylic_acid:
        return False, "No carboxylic acid group found"

    # Check for methyl branches
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetDegree() == 4:
            # Check if the carbon has exactly three hydrogen atoms (methyl group)
            hydrogens = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetSymbol() == 'H')
            if hydrogens == 3:
                # Check if the methyl carbon is connected to a carbon chain
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'C' and neighbor.GetDegree() > 1:
                        break
                else:
                    continue
                break
    else:
        return False, "No methyl branches found"

    # Check if there are only methyl branches
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetDegree() == 4:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C' and neighbor.GetDegree() > 1:
                    return True, "Methyl-branched fatty acid"

    return False, "Contains non-methyl branches"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:62499',
                          'name': 'methyl-branched fatty acid',
                          'definition': 'Any branched-chain fatty acid '
                                        'containing methyl branches only.',
                          'parents': ['CHEBI:35819']},
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
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 12,
    'num_false_negatives': 12,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}