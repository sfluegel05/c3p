"""
Classifies: CHEBI:35903 oxo carboxylic acid anion
"""
from rdkit import Chem

def is_oxo_carboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is an oxo carboxylic acid anion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oxo carboxylic acid anion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    carboxylate_found = False
    oxo_group_found = False

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1:
            # Check for carboxylate group (COO-)
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 1 and neighbors[0].GetAtomicNum() == 6:
                carbon = neighbors[0]
                carbon_neighbors = carbon.GetNeighbors()
                if len(carbon_neighbors) == 3:
                    oxygens = [n for n in carbon_neighbors if n.GetAtomicNum() == 8 and n.GetIdx() != atom.GetIdx()]
                    if len(oxygens) == 2 and all(o.GetFormalCharge() == 0 for o in oxygens):
                        carboxylate_found = True

        if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == 0:
            # Check for oxo group (C=O)
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 1 and neighbors[0].GetAtomicNum() == 6:
                carbon = neighbors[0]
                if carbon.GetTotalValence() == 4:
                    oxo_group_found = True

    if carboxylate_found and oxo_group_found:
        return True, "Contains both carboxylate and oxo groups"
    elif not carboxylate_found:
        return False, "No carboxylate group found"
    elif not oxo_group_found:
        return False, "No oxo group found"
    else:
        return None, None


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35903',
                          'name': 'oxo carboxylic acid anion',
                          'definition': 'Any carboxylic acid anion containing '
                                        'at least one oxo group.',
                          'parents': ['CHEBI:29067']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 13-14: malformed \\N character escape (<string>, line '
             '1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}