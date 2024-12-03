"""
Classifies: CHEBI:134396 secondary allylic alcohol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_secondary_allylic_alcohol(smiles: str):
    """
    Determines if a molecule is a secondary allylic alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary allylic alcohol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate through all atoms to find hydroxyl groups
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and atom.GetDegree() == 1:
            # Check if the oxygen is bonded to a carbon
            neighbor = atom.GetNeighbors()[0]
            if neighbor.GetSymbol() == 'C':
                # Check if the carbon is an allylic carbon (bonded to a double bond)
                is_allylic = False
                for nbr in neighbor.GetNeighbors():
                    if nbr.GetSymbol() == 'C' and mol.GetBondBetweenAtoms(neighbor.GetIdx(), nbr.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        is_allylic = True
                        break
                if not is_allylic:
                    continue

                # Check if the carbon is secondary (bonded to one other carbon and one hydrogen)
                carbon_count = 0
                hydrogen_count = 0
                for nbr in neighbor.GetNeighbors():
                    if nbr.GetSymbol() == 'C':
                        carbon_count += 1
                    elif nbr.GetSymbol() == 'H':
                        hydrogen_count += 1

                if carbon_count == 1 and hydrogen_count == 1:
                    return True, "Secondary allylic alcohol found"

    return False, "No secondary allylic alcohol found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:134396',
                          'name': 'secondary allylic alcohol',
                          'definition': 'An allylic alcohol in which the '
                                        'carbon atom that links the double '
                                        'bond to the hydroxy group is also '
                                        'attached to one other carbon and one '
                                        'hydrogen.',
                          'parents': ['CHEBI:134361', 'CHEBI:35681']},
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
    'num_true_negatives': 15,
    'num_false_negatives': 15,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}