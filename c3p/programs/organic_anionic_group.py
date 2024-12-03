"""
Classifies: CHEBI:64775 organic anionic group
"""
from rdkit import Chem

def is_organic_anionic_group(smiles: str):
    """
    Determines if a molecule is an organic anionic group (an anionic group that contains carbon).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organic anionic group, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of carbon atoms
    has_carbon = any(atom.GetSymbol() == 'C' for atom in mol.GetAtoms())
    if not has_carbon:
        return False, "No carbon atoms found"

    # Check for the presence of anionic groups (e.g., negative charges)
    has_anionic_group = any(atom.GetFormalCharge() < 0 for atom in mol.GetAtoms())
    if not has_anionic_group:
        return False, "No anionic groups found"

    return True, "Molecule is an organic anionic group"

# Example usage:
smiles = "C([C@H](C(NCCC(NCCSC(=O)C[C@@H](CCCCCCCCCCCCC)C)=O)=O)O)(COP(OC[C@@H](C(*)=O)N*)([O-])=O)(C)C"
print(is_organic_anionic_group(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:64775',
                          'name': 'organic anionic group',
                          'definition': 'An anionic group that contains '
                                        'carbon.',
                          'parents': ['CHEBI:64767']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(True, 'Molecule is an organic anionic group')\n",
    'num_true_positives': 74,
    'num_false_positives': 2,
    'num_true_negatives': 18,
    'num_false_negatives': 0,
    'precision': 0.9736842105263158,
    'recall': 1.0,
    'f1': 0.9866666666666666,
    'accuracy': None}