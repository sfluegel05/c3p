"""
Classifies: CHEBI:131863 docosanoid
"""
from rdkit import Chem

def is_docosanoid(smiles: str):
    """
    Determines if a molecule is a docosanoid (oxygenated derivative of C22 polyunsaturated fatty acids).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a docosanoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of 22 carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    if num_carbons != 22:
        return False, f"Number of carbon atoms is {num_carbons}, not 22"

    # Check for the presence of oxygen atoms
    num_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O')
    if num_oxygens == 0:
        return False, "No oxygen atoms found"

    # Check for the presence of double bonds
    num_double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondTypeAsDouble() == 2.0)
    if num_double_bonds < 3:
        return False, f"Number of double bonds is {num_double_bonds}, less than 3"

    # Check if molecule is a fatty acid (presence of carboxylic acid group)
    carboxylic_acid = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetTotalDegree() == 3:
            neighbors = [neighbor.GetSymbol() for neighbor in atom.GetNeighbors()]
            if neighbors.count('O') == 2:
                carboxylic_acid = True
                break

    if not carboxylic_acid:
        return False, "No carboxylic acid group found"

    return True, "Molecule is a docosanoid"

# Example usage
smiles = "CC/C=C\C/C=C\C=C\C(/C=C/C=C/C(C(CCCCCC(=O)O)O)O)O"  # resolvin T3
print(is_docosanoid(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:131863',
                          'name': 'docosanoid',
                          'definition': 'Any oxygenated derivative of C22 '
                                        'polyunsaturated fatty acids, such as '
                                        'docosapentaenoic acid (DPA) and '
                                        'docosahexaenoic acid (DHA).',
                          'parents': ['CHEBI:15904', 'CHEBI:26208']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(True, 'Molecule is a docosanoid')\n",
    'num_true_positives': 10,
    'num_false_positives': 0,
    'num_true_negatives': 10,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}