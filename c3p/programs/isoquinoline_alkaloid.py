"""
Classifies: CHEBI:24921 isoquinoline alkaloid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_isoquinoline_alkaloid(smiles: str):
    """
    Determines if a molecule is an isoquinoline alkaloid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoquinoline alkaloid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the isoquinoline substructure
    isoquinoline_smarts = "c1ccc2ncccc2c1"
    isoquinoline = Chem.MolFromSmarts(isoquinoline_smarts)

    if mol.HasSubstructMatch(isoquinoline):
        # Check for the presence of nitrogen atoms as part of the alkaloid structure
        nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'N']
        if nitrogen_atoms:
            return True, "Isoquinoline alkaloid structure detected"
        else:
            return False, "Isoquinoline nucleus present but no nitrogen atoms detected, not an alkaloid"
    else:
        return False, "No isoquinoline nucleus detected"

# Example usage
smiles = "COc1ccc2C=CC3=CC(=O)[C@]4(C[C@@]3(CCN4C)c2c1O)OC"
print(is_isoquinoline_alkaloid(smiles))  # Example SMILES string for testing


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24921',
                          'name': 'isoquinoline alkaloid',
                          'definition': 'Any alkaloid that has a structure '
                                        'based on an isoquinoline nucleus. '
                                        'They are derived from the amino acids '
                                        'like tyrosine and phenylalanine.',
                          'parents': ['CHEBI:22315']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(False, 'No isoquinoline nucleus detected')\n",
    'num_true_positives': 1,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 37,
    'precision': 1.0,
    'recall': 0.02631578947368421,
    'f1': 0.05128205128205127,
    'accuracy': None}