"""
Classifies: CHEBI:26199 polyprenol
"""
"""
Classifies: CHEBI:26195 polyprenol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    A polyprenol is a member of the class of prenols with the general formula H-[CH2C(Me)=CHCH2]nOH, 
    where the carbon skeleton is composed of more than one isoprene units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for terminal hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[CH2][OH]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No terminal hydroxyl group found"

    # Define isoprene unit pattern with methyl group
    isoprene_pattern = Chem.MolFromSmarts("[CH3][CH]=[CH][CH2]")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    
    # Need at least 2 isoprene units
    if len(isoprene_matches) < 2:
        return False, f"Found {len(isoprene_matches)} isoprene units, need at least 2"

    # Check connectivity of isoprene units
    # Get all carbon atoms
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    
    # Check if carbons form a continuous chain
    # This is a simplified check - a more rigorous approach would be needed
    # for complex structures
    chain_length = 0
    current_atom = None
    
    # Find the terminal carbon with hydroxyl group
    for atom in carbon_atoms:
        if atom.GetTotalNumHs() == 2 and len([n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]) == 1:
            current_atom = atom
            break
    
    if current_atom is None:
        return False, "Could not find terminal carbon with hydroxyl group"
    
    # Traverse the chain
    visited = set()
    while current_atom is not None:
        visited.add(current_atom.GetIdx())
        # Get next carbon in chain
        neighbors = [n for n in current_atom.GetNeighbors() 
                    if n.GetAtomicNum() == 6 and n.GetIdx() not in visited]
        if len(neighbors) == 0:
            break
        current_atom = neighbors[0]
        chain_length += 1
    
    # Minimum chain length for 2 isoprene units is 8 carbons
    if chain_length < 8:
        return False, "Chain too short for polyprenol structure"

    # Check double bond pattern
    double_bonds = [bond for bond in mol.GetBonds() 
                   if bond.GetBondType() == Chem.BondType.DOUBLE]
    if len(double_bonds) < 2:
        return False, "Not enough double bonds for polyprenol structure"

    # Check alternating double bond pattern
    # This is a simplified check
    prev_double = False
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            if prev_double:
                return False, "Double bonds not alternating properly"
            prev_double = True
        else:
            prev_double = False

    return True, "Contains terminal hydroxyl group and chain of more than one isoprene units with alternating double bonds"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26195',
                          'name': 'polyprenol',
                          'definition': 'Any member of the class of prenols '
                                        'possessing the general formula '
                                        'H-[CH2C(Me)=CHCH2]nOH in which the '
                                        'carbon skeleton is composed of more '
                                        'than one isoprene units.',
                          'parents': ['CHEBI:26195', 'CHEBI:26195']},
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