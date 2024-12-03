"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify carbonyl groups (C=O)
    carbonyl_groups = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and any(neigh.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(atom.GetIdx(), neigh.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE for neigh in atom.GetNeighbors())]

    if not carbonyl_groups:
        return False, "No carbonyl groups found"

    for carbonyl in carbonyl_groups:
        carbonyl_idx = carbonyl.GetIdx()
        neighbors = carbonyl.GetNeighbors()

        # Check for alpha-hydroxy group
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() == 6:  # Check if neighbor is a carbon
                alpha_carbon_idx = neighbor.GetIdx()
                alpha_carbon_neighbors = neighbor.GetNeighbors()

                for alpha_carbon_neighbor in alpha_carbon_neighbors:
                    if alpha_carbon_neighbor.GetAtomicNum() == 8 and alpha_carbon_neighbor.GetTotalNumHs() > 0:
                        return True, "Alpha-hydroxy ketone found"

    return False, "No alpha-hydroxy ketone structure found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:139588',
                          'name': 'alpha-hydroxy ketone',
                          'definition': 'An alpha-oxyketone that has a hydroxy '
                                        'group as the alpha-oxy moiety.',
                          'parents': ['CHEBI:30879', 'CHEBI:52396']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[19:23:43] SMILES Parse Error: unclosed ring for input: '
             "'C[C@H]1C\\C=C\\[C@H]2[C@H](O)C(C)=C(C)[C@H]3[C@H](Cc4c[nH]c5ccccc45)NC(=O)[C@@]23C(=O)\\C=C\\C(=O)[C@H](O)\\C(C)=C\x01'\n",
    'stdout': '',
    'num_true_positives': 50,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 1,
    'precision': 0.9803921568627451,
    'recall': 0.9803921568627451,
    'f1': 0.9803921568627451,
    'accuracy': None}