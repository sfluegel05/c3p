"""
Classifies: CHEBI:26244 prenols
"""
"""
Classifies: CHEBI:26195 prenol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.
    A prenol is an alcohol with the general formula H-[CH2C(Me)=CHCH2]nOH, 
    where the carbon skeleton is composed of one or more isoprene units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a hydroxyl group (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found"

    # Check for the isoprene unit pattern (CH2C(Me)=CHCH2)
    isoprene_pattern = Chem.MolFromSmarts("[CH2][CH]([CH3])=[CH][CH2]")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 1:
        return False, "No isoprene units found"

    # Check if the molecule has a repeating isoprene pattern
    # We can count the number of isoprene units and ensure they are connected
    # in a way that forms a continuous chain.
    # This is a simplified check and may not cover all cases.
    # For a more robust check, we would need to verify the connectivity of the isoprene units.
    # However, this should work for most prenols.

    # Count the number of isoprene units
    n_isoprene = len(isoprene_matches)
    if n_isoprene < 1:
        return False, "No isoprene units found"

    # Check if the molecule has a long carbon chain (typical of prenols)
    # Prenols usually have a significant number of carbons due to the repeating isoprene units.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 5:
        return False, "Too few carbons for a prenol"

    # Check for the presence of double bonds (isoprene units have double bonds)
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bond_count < 1:
        return False, "No double bonds found (isoprene units have double bonds)"

    # Check if the hydroxyl group is at the end of the chain
    # This is a simplified check and may not cover all cases.
    # For a more robust check, we would need to verify the position of the hydroxyl group.
    # However, this should work for most prenols.
    hydroxyl_atom = mol.GetSubstructMatch(hydroxyl_pattern)[0]
    hydroxyl_neighbors = mol.GetAtomWithIdx(hydroxyl_atom).GetNeighbors()
    if len(hydroxyl_neighbors) != 1:
        return False, "Hydroxyl group is not at the end of the chain"

    return True, "Contains a hydroxyl group and one or more isoprene units"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26195',
                          'name': 'prenol',
                          'definition': 'Any alcohol possessing the general '
                                        'formula H-[CH2C(Me)=CHCH2]nOH in '
                                        'which the carbon skeleton is composed '
                                        'of one or more isoprene units '
                                        '(biogenetic precursors of the '
                                        'isoprenoids).',
                          'parents': ['CHEBI:23888', 'CHEBI:26194']},
    'config': {   'llm_model_name': 'lbl/claude-sonnet',
                  'f1_threshold': 0.8,
                  'max_attempts': 5,
                  'max_positive_instances': None,
                  'max_positive_to_test': None,
                  'max_negative_to_test': None,
                  'max_positive_in_prompt': 50,
                  'max_negative_in_prompt': 20,
                  'max_instances_in_prompt': 100,
                  'test_proportion': 0.1},
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199}