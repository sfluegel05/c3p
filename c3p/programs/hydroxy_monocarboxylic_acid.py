"""
Classifies: CHEBI:35868 hydroxy monocarboxylic acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_hydroxy_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy monocarboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxy monocarboxylic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (COOH)
    carboxylic_acid = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # carbon
            if any(neigh.GetAtomicNum() == 8 and neigh.GetTotalNumHs() == 1 for neigh in atom.GetNeighbors()):
                if any(neigh.GetAtomicNum() == 8 and neigh.GetTotalNumHs() == 0 for neigh in atom.GetNeighbors()):
                    carboxylic_acid = True
                    break

    if not carboxylic_acid:
        return False, "No carboxylic acid group found"

    # Check for hydroxy group (OH)
    hydroxy = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # oxygen
            if atom.GetTotalNumHs() == 1:
                hydroxy = True
                break

    if not hydroxy:
        return False, "No hydroxy group found"

    # Ensure the hydroxy group is not part of the carboxylic acid group
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # carbon
            if any(neigh.GetAtomicNum() == 8 and neigh.GetTotalNumHs() == 1 for neigh in atom.GetNeighbors()):
                if any(neigh.GetAtomicNum() == 8 and neigh.GetTotalNumHs() == 0 for neigh in atom.GetNeighbors()):
                    # This is the carboxylic acid carbon, check its neighbors
                    for neigh in atom.GetNeighbors():
                        if neigh.GetAtomicNum() == 8 and neigh.GetTotalNumHs() == 1:
                            hydroxy = False
                            break
                    if not hydroxy:
                        break

    if not hydroxy:
        return False, "Hydroxy group is part of the carboxylic acid group"

    return True, "Molecule is a hydroxy monocarboxylic acid"

# Example usage
smiles = "OC(=O)C(\O)=C\C=C"  # (2Z)-2-hydroxypenta-2,4-dienoic acid
print(is_hydroxy_monocarboxylic_acid(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35868',
                          'name': 'hydroxy monocarboxylic acid',
                          'definition': 'Any monocarboxylic acid which also '
                                        'contains a separate (alcoholic or '
                                        'phenolic) hydroxy substituent.',
                          'parents': ['CHEBI:24669', 'CHEBI:25384']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 8-9: malformed \\N character escape (<string>, line 1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}