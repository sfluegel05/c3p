"""
Classifies: CHEBI:60004 mixture
"""
from rdkit import Chem

def is_mixture(smiles: str):
    """
    Determines if a molecule is a mixture (composed of multiple molecules, at least two of which are of a different kind).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mixture, False otherwise
        str: Reason for classification
    """
    mols = Chem.MolFromSmiles(smiles, sanitize=False)
    if mols is None:
        return False, "Invalid SMILES string"

    # Split the SMILES string into individual components
    components = smiles.split('.')
    
    if len(components) < 2:
        return False, "SMILES string does not represent a mixture"

    unique_molecules = set()
    
    for component in components:
        mol = Chem.MolFromSmiles(component)
        if mol is None:
            return False, f"Invalid component in SMILES string: {component}"
        unique_molecules.add(Chem.MolToSmiles(mol, isomericSmiles=True))

    if len(unique_molecules) < 2:
        return False, "Not enough unique molecules to be considered a mixture"

    return True, f"Mixture identified with {len(unique_molecules)} unique components"

# Example usage
smiles = "OC[C@H]([C@@H](O)[C@H](O)[C@@H](O)C([O-])=O)O.[Fe+3].OC[C@H]([C@@H](O)[C@H](O)[C@@H](O)C([O-])=O)O.OC[C@H]([C@@H](O)[C@H](O)[C@@H](O)C([O-])=O)O"
print(is_mixture(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:60004',
                          'name': 'mixture',
                          'definition': 'A mixture is a chemical substance '
                                        'composed of multiple molecules, at '
                                        'least two of which are of a different '
                                        'kind.',
                          'parents': ['CHEBI:59999']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(True, 'Mixture identified with 2 unique components')\n",
    'num_true_positives': 7,
    'num_false_positives': 3,
    'num_true_negatives': 9,
    'num_false_negatives': 5,
    'precision': 0.7,
    'recall': 0.5833333333333334,
    'f1': 0.6363636363636365,
    'accuracy': None}