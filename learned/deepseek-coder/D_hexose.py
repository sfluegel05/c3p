"""
Classifies: CHEBI:4194 D-hexose
"""
"""
Classifies: CHEBI:15903 D-hexose
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    A D-hexose is a hexose that has D-configuration at position 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a D-hexose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule has a ring structure
    if not mol.GetRingInfo().NumRings():
        # Allow linear forms (aldehydo forms)
        pass

    # Check for at least 4 oxygen atoms (hydroxyl or ring oxygen)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 4:
        return False, "Not enough oxygen atoms for a hexose"

    # Find all potential position 5 carbons
    # Look for carbons with a hydroxyl group and chirality
    potential_position_5 = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            # Check for hydroxyl group
            has_hydroxyl = False
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:
                    has_hydroxyl = True
                    break
            if has_hydroxyl and atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                potential_position_5.append(atom)

    # Check if any potential position 5 carbon has D-configuration
    for carbon in potential_position_5:
        # Check stereochemistry - D-configuration is R in Cahn-Ingold-Prelog
        try:
            Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
            if carbon.GetProp('_CIPCode') == 'R':
                return True, "Hexose with D-configuration at position 5"
        except:
            # If stereochemistry assignment fails, try alternative method
            if carbon.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
                return True, "Hexose with D-configuration at position 5"

    return False, "No carbon with D-configuration found at position 5"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:15903',
                          'name': 'D-hexose',
                          'definition': 'A hexose that has D-configuration at position 5.',
                          'parents': ['CHEBI:15903', 'CHEBI:15903']},
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