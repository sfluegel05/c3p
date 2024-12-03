"""
Classifies: CHEBI:28965 dicarboxylic acid dianion
"""
from rdkit import Chem


def is_dicarboxylic_acid_dianion(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid dianion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dicarboxylic acid dianion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all carboxylate groups (COO-)
    carboxylate_groups = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetFormalCharge() == 0:
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 3:
                oxygens = [n for n in neighbors if n.GetSymbol() == 'O' and n.GetFormalCharge() == -1]
                if len(oxygens) == 2:
                    carboxylate_groups.append(atom.GetIdx())

    if len(carboxylate_groups) < 2:
        return False, "Less than two carboxylate groups found"

    # Ensure the carboxylate groups are connected through a carbon chain
    for i in range(len(carboxylate_groups)):
        for j in range(i + 1, len(carboxylate_groups)):
            path = Chem.rdmolops.GetShortestPath(mol, carboxylate_groups[i], carboxylate_groups[j])
            if all(mol.GetAtomWithIdx(idx).GetSymbol() == 'C' for idx in path):
                return True, "Dicarboxylic acid dianion found"

    return False, "Carboxylate groups are not connected through a carbon chain"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:28965',
                          'name': 'dicarboxylic acid dianion',
                          'definition': 'A carboxylic acid dianion obtained by '
                                        'deprotonation of both carboxy groups '
                                        'of any dicarboxylic acid.',
                          'parents': ['CHEBI:35693', 'CHEBI:38716']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 10-11: malformed \\N character escape (<string>, line '
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