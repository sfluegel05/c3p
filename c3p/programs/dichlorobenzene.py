"""
Classifies: CHEBI:23697 dichlorobenzene
"""
from rdkit import Chem

def is_dichlorobenzene(smiles: str):
    """
    Determines if a molecule is a dichlorobenzene.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dichlorobenzene, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for benzene ring
    benzene = Chem.MolFromSmarts('c1ccccc1')
    if not mol.HasSubstructMatch(benzene):
        return False, "No benzene ring found"

    # Check for exactly two chlorine atoms
    chlorine_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'Cl']
    if len(chlorine_atoms) != 2:
        return False, f"Incorrect number of chlorine atoms: found {len(chlorine_atoms)}, expected 2"

    # Check if both chlorine atoms are attached to the benzene ring
    chlorine_on_benzene = 0
    for atom in chlorine_atoms:
        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() == 'C' and neighbor.IsInRing():
                chlorine_on_benzene += 1
                break

    if chlorine_on_benzene != 2:
        return False, "Chlorine atoms are not both attached to the benzene ring"

    return True, "Molecule is a dichlorobenzene"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23697',
                          'name': 'dichlorobenzene',
                          'definition': 'Any member of the class of '
                                        'chlorobenzenes carrying two chloro '
                                        'groups at unspecified positions.',
                          'parents': ['CHEBI:23132']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 22-23: malformed \\N character escape (<string>, line '
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