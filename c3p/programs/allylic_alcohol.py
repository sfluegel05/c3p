"""
Classifies: CHEBI:134361 allylic alcohol
"""
from rdkit import Chem

def is_allylic_alcohol(smiles: str):
    """
    Determines if a molecule is an allylic alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an allylic alcohol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate through all atoms to find hydroxyl groups
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and any(neigh.GetSymbol() == 'H' for neigh in atom.GetNeighbors()):
            # Check if the oxygen is part of a hydroxyl group
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C' and neighbor.GetDegree() == 3:
                    # Check if the carbon is saturated (single bonds only)
                    if all(bond.GetBondType() == Chem.rdchem.BondType.SINGLE for bond in neighbor.GetBonds()):
                        # Check if the carbon is adjacent to a double bond
                        for neigh_neigh in neighbor.GetNeighbors():
                            if neigh_neigh.GetSymbol() == 'C' and any(bond.GetBondType() == Chem.rdchem.BondType.DOUBLE for bond in neigh_neigh.GetBonds()):
                                return True, "Molecule is an allylic alcohol"

    return False, "Molecule is not an allylic alcohol"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:134361',
                          'name': 'allylic alcohol',
                          'definition': 'An alcohol where the hydroxy group is '
                                        'attached to a saturated carbon atom '
                                        'adjacent to a double bond (R groups '
                                        'may be H, organyl, etc.).',
                          'parents': ['CHEBI:30879', 'CHEBI:78840']},
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
    'num_true_negatives': 20,
    'num_false_negatives': 23,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}