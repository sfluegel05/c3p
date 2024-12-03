"""
Classifies: CHEBI:58953 saturated fatty acid anion
"""
from rdkit import Chem

def is_saturated_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid anion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saturated fatty acid anion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the carboxylate group (anion form of carboxylic acid)
    carboxylate = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'C' and atom.GetTotalDegree() == 3]
    if not carboxylate:
        return False, "No carboxylate group found"

    carboxylate = carboxylate[0]
    neighbors = carboxylate.GetNeighbors()
    if not any(neighbor.GetSymbol() == 'O' and neighbor.GetFormalCharge() == -1 for neighbor in neighbors):
        return False, "No carboxylate anion found"

    # Check for C-C unsaturation
    for bond in mol.GetBonds():
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        if bond.GetBondType() in (Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE):
            # Check if either atom is carbon
            if begin_atom.GetSymbol() == 'C' and end_atom.GetSymbol() == 'C':
                return False, "C-C unsaturation found"

    return True, "Saturated fatty acid anion"

# Example usage:
# print(is_saturated_fatty_acid_anion("[O-]C(=O)CCCCCCCCC(CCCCCCCC)O"))
# print(is_saturated_fatty_acid_anion("CCCCCCCCCCCCCCCC([O-])=O"))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:58953',
                          'name': 'saturated fatty acid anion',
                          'definition': 'Any fatty acid anion in which there '
                                        'is no C-C unsaturation.',
                          'parents': ['CHEBI:28868']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 11,
    'num_false_positives': 2,
    'num_true_negatives': 9,
    'num_false_negatives': 0,
    'precision': 0.8461538461538461,
    'recall': 1.0,
    'f1': 0.9166666666666666,
    'accuracy': None}