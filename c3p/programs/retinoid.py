"""
Classifies: CHEBI:26537 retinoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_retinoid(smiles: str):
    """
    Determines if a molecule is a retinoid (Oxygenated derivatives of 3,7-dimethyl-1-(2,6,6-trimethylcyclohex-1-enyl)nona-1,3,5,7-tetraene and derivatives thereof).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a retinoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core structure of a retinoid
    core_smarts = "C1=C(C)CCCC1(C)C=C(C)C=C(C)C=C(C)C=C(C)C"
    core_mol = Chem.MolFromSmarts(core_smarts)

    if core_mol is None:
        return False, "Invalid core structure"

    # Check if the molecule contains the core structure
    if not mol.HasSubstructMatch(core_mol):
        return False, "Molecule does not contain the retinoid core structure"

    # Check for oxygenated derivatives
    oxygen_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'O']
    if not oxygen_atoms:
        return False, "Molecule does not contain any oxygen atoms"

    return True, "Molecule is a retinoid"

# Example usage:
# smiles = "C\C(\C=C\C=C(\C=C\C1=C(C)CCCC1(C)C)/C1=CC=C(C)C=C1)=C/C(O)=O"  # SR11302
# print(is_retinoid(smiles))  # Should return (True, "Molecule is a retinoid")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26537',
                          'name': 'retinoid',
                          'definition': 'Oxygenated derivatives of '
                                        '3,7-dimethyl-1-(2,6,6-trimethylcyclohex-1-enyl)nona-1,3,5,7-tetraene '
                                        'and derivatives thereof.',
                          'parents': ['CHEBI:23849']},
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