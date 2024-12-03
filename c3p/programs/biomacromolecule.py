"""
Classifies: CHEBI:33694 biomacromolecule
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_biomacromolecule(smiles: str):
    """
    Determines if a molecule is a biomacromolecule (a macromolecule formed by a living organism).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a biomacromolecule, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate molecular weight
    mol_weight = Descriptors.MolWt(mol)

    # Check for high molecular weight (indicative of macromolecules)
    if mol_weight < 1000:
        return False, "Molecular weight is too low to be a macromolecule"

    # Check for presence of repeating units (common in biomacromolecules like proteins, polysaccharides, etc.)
    repeating_units = ["C(O)C(O)", "NC(=O)C", "C(O)C(O)C(O)"]  # Common repeating units in biomacromolecules
    mol_smiles = Chem.MolToSmiles(mol)
    for unit in repeating_units:
        if mol_smiles.count(unit) > 1:
            return True, f"Contains repeating unit: {unit}"

    return False, "Does not contain common repeating units of biomacromolecules"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33694',
                          'name': 'biomacromolecule',
                          'definition': 'A macromolecule formed by a living '
                                        'organism.',
                          'parents': ['CHEBI:33839', 'CHEBI:50860']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 118-119: malformed \\N character escape (<string>, line '
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