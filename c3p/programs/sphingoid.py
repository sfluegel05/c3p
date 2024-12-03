"""
Classifies: CHEBI:35785 sphingoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid (Sphinganine, its homologs and stereoisomers, and the hydroxy and unsaturated derivatives of these compounds).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a long aliphatic chain (at least 12 carbons)
    chain_length = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # carbon atom
            if len([nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]) == 2:
                chain_length += 1

    if chain_length < 12:
        return False, "Aliphatic chain is too short"

    # Check for the presence of hydroxyl groups
    hydroxyl_groups = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1]
    if len(hydroxyl_groups) == 0:
        return False, "No hydroxyl groups found"

    # Check for the presence of an amino group
    amino_groups = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7 and atom.GetDegree() in [1, 2, 3]]
    if len(amino_groups) == 0:
        return False, "No amino groups found"

    return True, "Molecule is a sphingoid"

# Example usage:
# print(is_sphingoid("CCCCCCCCCCCCCCC[C@@H](O)CN"))  # 1-deoxymethylsphinganine
# print(is_sphingoid("CCCCCCCCCCCCCCC[C@@H](O)[C@H](C)N"))  # 1-deoxysphinganine


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35785',
                          'name': 'sphingoid',
                          'definition': 'Sphinganine, its homologs and '
                                        'stereoisomers, and the hydroxy and '
                                        'unsaturated derivatives of these '
                                        'compounds.',
                          'parents': ['CHEBI:26739']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 12,
    'num_false_positives': 10,
    'num_true_negatives': 2,
    'num_false_negatives': 0,
    'precision': 0.5454545454545454,
    'recall': 1.0,
    'f1': 0.7058823529411764,
    'accuracy': None}