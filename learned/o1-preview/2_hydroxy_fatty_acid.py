"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
"""
Classifies: CHEBI:15617 2-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    A 2-hydroxy fatty acid is any fatty acid with a hydroxy functional group in the alpha- or 2-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxylic_acid_matches:
        return False, "No carboxylic acid group found"
    
    # Assume first carboxylic acid group
    carboxyl_carbon_idx = carboxylic_acid_matches[0][0]  # Index of carbon in -C(=O)OH
    
    # Get the alpha carbon (adjacent to carboxyl carbon)
    carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_carbon_idx)
    neighbors = [atom for atom in carboxyl_carbon.GetNeighbors() if atom.GetIdx() != carboxyl_carbon_idx]
    alpha_carbons = [atom for atom in neighbors if atom.GetAtomicNum() == 6]  # Carbon atoms
    if not alpha_carbons:
        return False, "No alpha carbon adjacent to carboxyl carbon"
    alpha_carbon = alpha_carbons[0]
    
    # Check if alpha carbon has a hydroxy group attached
    alpha_carbon_idx = alpha_carbon.GetIdx()
    alpha_neighbors = alpha_carbon.GetNeighbors()
    has_hydroxy = False
    for neighbor in alpha_neighbors:
        if neighbor.GetAtomicNum() == 8:  # Oxygen atom
            # Check if oxygen is part of a hydroxy group (-OH)
            if mol.GetBondBetweenAtoms(alpha_carbon_idx, neighbor.GetIdx()).GetBondType() == Chem.rdchem.BondType.SINGLE:
                oxy_neighbors = neighbor.GetNeighbors()
                if len(oxy_neighbors) == 1:
                    has_hydroxy = True
                    break
    if not has_hydroxy:
        return False, "Alpha carbon does not have a hydroxy group attached"
    
    # Optionally, check if the molecule has a long aliphatic chain (typical of fatty acids)
    # Count the total number of carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
        return False, f"Molecule has {c_count} carbon atoms, less than typical fatty acids"
    
    # Check if molecule is mostly linear (fatty acids are usually unbranched)
    ring_info = mol.GetRingInfo()
    if ring_info.IsInitialized() and ring_info.NumRings() > 0:
        return False, "Molecule contains ring structures, which is atypical for fatty acids"
    
    return True, "Molecule is a 2-hydroxy fatty acid with hydroxy group at alpha position"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:15617',
        'name': '2-hydroxy fatty acid',
        'definition': 'Any fatty acid with a hydroxy functional group in the alpha- or 2-position.',
        'parents': ['CHEBI:25360', 'CHEBI:33853']
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    },
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': None,
    'num_false_positives': None,
    'num_true_negatives': None,
    'num_false_negatives': None,
    'num_negatives': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None
}