"""
Classifies: CHEBI:36244 dicarboxylic acid monoester
"""
from rdkit import Chem

def is_dicarboxylic_acid_monoester(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid monoester.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dicarboxylic acid monoester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    carboxyl_groups = []
    ester_groups = []

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            neighbors = [neighbor.GetAtomicNum() for neighbor in atom.GetNeighbors()]
            if neighbors.count(8) == 2:  # Two oxygen neighbors
                oxygens = [neighbor for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 8]
                if any(o.GetTotalNumHs() == 1 for o in oxygens):  # Carboxyl group
                    carboxyl_groups.append(atom.GetIdx())
                if any(o.GetTotalNumHs() == 0 and o.GetDegree() == 2 for o in oxygens):  # Ester group
                    ester_groups.append(atom.GetIdx())

    if len(carboxyl_groups) < 2:
        return False, "Less than two carboxyl groups found"
    if len(ester_groups) < 1:
        return False, "No ester groups found"

    return True, "Dicarboxylic acid monoester identified"

# Example usage:
# smiles = "CC(=O)[C@H]1CCC2[C@@]1(CCC3C2CC=C4[C@@]3(CC[C@@H](C4)OC(=O)CCC(=O)O)C)C"
# print(is_dicarboxylic_acid_monoester(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36244',
                          'name': 'dicarboxylic acid monoester',
                          'definition': 'A monoester of a dicarboxylic acid.',
                          'parents': [   'CHEBI:131927',
                                         'CHEBI:33308',
                                         'CHEBI:33575']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 12-13: malformed \\N character escape (<string>, line '
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