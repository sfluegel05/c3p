"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
from rdkit import Chem

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid (an oxoacid containing three carboxy groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tricarboxylic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    carboxylic_acid_count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            # Check if the carbon is part of a carboxylic acid group (C(=O)O)
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 3:
                o_atoms = [n for n in neighbors if n.GetSymbol() == 'O']
                if len(o_atoms) == 2:
                    double_bonded_o = any(bond.GetBondType() == Chem.rdchem.BondType.DOUBLE for bond in atom.GetBonds() if bond.GetOtherAtom(atom).GetSymbol() == 'O')
                    single_bonded_o = any(bond.GetBondType() == Chem.rdchem.BondType.SINGLE for bond in atom.GetBonds() if bond.GetOtherAtom(atom).GetSymbol() == 'O')
                    if double_bonded_o and single_bonded_o:
                        carboxylic_acid_count += 1

    if carboxylic_acid_count == 3:
        return True, "Molecule contains three carboxy groups"
    else:
        return False, f"Molecule contains {carboxylic_acid_count} carboxy groups"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:27093',
                          'name': 'tricarboxylic acid',
                          'definition': 'An oxoacid containing three carboxy '
                                        'groups.',
                          'parents': ['CHEBI:33575']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 14-15: malformed \\N character escape (<string>, line '
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