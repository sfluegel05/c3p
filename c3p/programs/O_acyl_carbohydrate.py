"""
Classifies: CHEBI:52782 O-acyl carbohydrate
"""
from rdkit import Chem

def is_O_acyl_carbohydrate(smiles: str):
    """
    Determines if a molecule is an O-acyl carbohydrate.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an O-acyl carbohydrate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all ester groups (O=C-O)
    ester_groups = []
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            if atom1.GetAtomicNum() == 8 and atom2.GetAtomicNum() == 6:
                # Check if there is a single bond from the carbon to another oxygen
                for nbr in atom2.GetNeighbors():
                    if nbr.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(atom2.GetIdx(), nbr.GetIdx()).GetBondType() == Chem.rdchem.BondType.SINGLE:
                        ester_groups.append((atom1, atom2, nbr))
            elif atom2.GetAtomicNum() == 8 and atom1.GetAtomicNum() == 6:
                # Check if there is a single bond from the carbon to another oxygen
                for nbr in atom1.GetNeighbors():
                    if nbr.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(atom1.GetIdx(), nbr.GetIdx()).GetBondType() == Chem.rdchem.BondType.SINGLE:
                        ester_groups.append((atom2, atom1, nbr))

    if not ester_groups:
        return False, "No ester groups found"

    # Check if the ester group is attached to a carbohydrate
    for ester_group in ester_groups:
        carbon = ester_group[1]
        oxygen = ester_group[2]
        for nbr in oxygen.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != carbon.GetIdx():
                # Check if the neighbor carbon is part of a carbohydrate
                if any(atom.GetAtomicNum() == 6 and atom.GetDegree() == 3 for atom in nbr.GetNeighbors()):
                    return True, "O-acyl carbohydrate found"

    return False, "No O-acyl carbohydrate structure found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:52782',
                          'name': 'O-acyl carbohydrate',
                          'definition': 'A carbohydrate derivative in which '
                                        'the hydrogen atom of at least one '
                                        'alcoholic hydroxy group of a '
                                        'carbohydrate has been replaced by an '
                                        'acyl substituent.',
                          'parents': ['CHEBI:63299']},
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
    'num_true_positives': 11,
    'num_false_positives': 2,
    'num_true_negatives': 18,
    'num_false_negatives': 34,
    'precision': 0.8461538461538461,
    'recall': 0.24444444444444444,
    'f1': 0.37931034482758624,
    'accuracy': None}