"""
Classifies: CHEBI:58950 very long-chain fatty acid anion
"""
from rdkit import Chem

def is_very_long_chain_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid anion (chain length greater than C22).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acid anion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylate anion group ([O-]C=O)
    carboxylate_anion = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetTotalNumHs() == 0:
            neighbors = [nbr for nbr in atom.GetNeighbors()]
            if len(neighbors) == 2:
                if neighbors[0].GetSymbol() == 'O' and neighbors[0].GetFormalCharge() == -1:
                    if neighbors[1].GetSymbol() == 'O' and neighbors[1].GetFormalCharge() == 0:
                        carboxylate_anion = True
                        break

    if not carboxylate_anion:
        return False, "No carboxylate anion group found"

    # Count the number of carbon atoms in the longest chain
    chains = mol.GetRingInfo().AtomRings()
    if chains:
        return False, "Molecule contains rings"

    carbon_count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            carbon_count += 1

    if carbon_count > 22:
        return True, f"Very long-chain fatty acid anion with {carbon_count} carbon atoms"
    else:
        return False, f"Chain length is {carbon_count}, which is not greater than C22"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:58950',
                          'name': 'very long-chain fatty acid anion',
                          'definition': 'Any fatty acid anion with a chain '
                                        'length greater than C22.',
                          'parents': ['CHEBI:28868']},
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
    'num_true_negatives': 10,
    'num_false_negatives': 10,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}