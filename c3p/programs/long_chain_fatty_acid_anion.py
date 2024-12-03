"""
Classifies: CHEBI:57560 long-chain fatty acid anion
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acid anion (C13 to C22).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acid anion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylate group ([O-]C=O)
    carboxylate_group = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetFormalCharge() == 0:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'O' and neighbor.GetFormalCharge() == -1:
                    carboxylate_group = True
                    break
            if carboxylate_group:
                break
    
    if not carboxylate_group:
        return False, "No carboxylate group found"

    # Check chain length (C13 to C22)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')

    if not (13 <= carbon_count <= 22):
        return False, "Chain length is not between C13 and C22"

    return True, "Valid long-chain fatty acid anion"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:57560',
                          'name': 'long-chain fatty acid anion',
                          'definition': 'A fatty acid anion with a chain '
                                        'length of C13 to C22.',
                          'parents': ['CHEBI:28868']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 46,
    'num_false_positives': 5,
    'num_true_negatives': 15,
    'num_false_negatives': 8,
    'precision': 0.9019607843137255,
    'recall': 0.8518518518518519,
    'f1': 0.8761904761904761,
    'accuracy': None}