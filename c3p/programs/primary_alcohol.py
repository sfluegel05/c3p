"""
Classifies: CHEBI:15734 primary alcohol
"""
from rdkit import Chem

def is_primary_alcohol(smiles: str):
    """
    Determines if a molecule is a primary alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary alcohol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all hydroxyl groups
    hydroxyl_groups = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'O' and any(neigh.GetSymbol() == 'H' for neigh in atom.GetNeighbors())]
    if not hydroxyl_groups:
        return False, "No hydroxyl groups found"

    for hydroxyl in hydroxyl_groups:
        # Check if the oxygen is bonded to a carbon
        carbon_neighbors = [neigh for neigh in hydroxyl.GetNeighbors() if neigh.GetSymbol() == 'C']
        if not carbon_neighbors:
            continue

        carbon = carbon_neighbors[0]

        # Check if the carbon is saturated (sp3 hybridized)
        if carbon.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue

        # Check the number of hydrogen atoms attached to the carbon
        hydrogen_count = sum(1 for neigh in carbon.GetNeighbors() if neigh.GetSymbol() == 'H')

        # Check the number of carbon atoms attached to the carbon
        carbon_count = sum(1 for neigh in carbon.GetNeighbors() if neigh.GetSymbol() == 'C')

        if hydrogen_count == 2 and carbon_count == 1:
            return True, "Primary alcohol with one carbon and two hydrogen atoms attached to the carbon bearing the hydroxyl group"
        elif hydrogen_count == 3 and carbon_count == 0:
            return True, "Primary alcohol with three hydrogen atoms attached to the carbon bearing the hydroxyl group"

    return False, "No primary alcohol groups found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:15734',
                          'name': 'primary alcohol',
                          'definition': 'A primary alcohol is a compound in '
                                        'which a hydroxy group, -OH, is '
                                        'attached to a saturated carbon atom '
                                        'which has either three hydrogen atoms '
                                        'attached to it or only one other '
                                        'carbon atom and two hydrogen atoms '
                                        'attached to it.',
                          'parents': ['CHEBI:30879']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[19:30:15] SMILES Parse Error: syntax error while parsing: '
             'C(=C/C=C/CCCCC)\\C=C\\C(=C\x01/C(C(NC1=O)CO)=O)\\O\n'
             '[19:30:15] SMILES Parse Error: Failed parsing SMILES '
             "'C(=C/C=C/CCCCC)\\C=C\\C(=C\x01/C(C(NC1=O)CO)=O)\\O' for input: "
             "'C(=C/C=C/CCCCC)\\C=C\\C(=C\x01/C(C(NC1=O)CO)=O)\\O'\n",
    'stdout': '',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 102,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}