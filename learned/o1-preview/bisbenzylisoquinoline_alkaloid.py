"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
# bisbenzylisoquinoline_alkaloid.py
"""
Classifies: bisbenzylisoquinoline alkaloid
Definition: A type of benzylisoquinoline alkaloid whose structures are built up of two benzylisoquinoline units linked by ether bridges. Various structural patterns resulting from additional bridging between the two units by direct carbon-carbon bridging or by methylenedioxy groups are common.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_bisbenzylisoquinoline_alkaloid(smiles: str):
    """
    Determines if a molecule is a bisbenzylisoquinoline alkaloid based on its SMILES string.
    A bisbenzylisoquinoline alkaloid consists of two benzylisoquinoline units linked by ether bridges.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a bisbenzylisoquinoline alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define isoquinoline pattern
    isoquinoline_smarts = 'c1ccc2ncccc2c1'  # Isoquinoline ring
    isoquinoline_pattern = Chem.MolFromSmarts(isoquinoline_smarts)
    if isoquinoline_pattern is None:
        return False, "Invalid isoquinoline SMARTS pattern"
    
    # Find isoquinoline units
    isoquinoline_matches = mol.GetSubstructMatches(isoquinoline_pattern)
    num_isoquinoline_units = len(isoquinoline_matches)
    if num_isoquinoline_units < 2:
        return False, f"Found {num_isoquinoline_units} isoquinoline unit(s), need at least 2"
    
    # Collect isoquinoline units atom indices
    isoquinoline_units = [set(match) for match in isoquinoline_matches]
    
    # Identify ether bridges (C-O-C)
    ether_oxygen_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # Oxygen atom
            neighbors = atom.GetNeighbors()
            carbon_neighbors = [n for n in neighbors if n.GetAtomicNum() == 6]
            if len(carbon_neighbors) == 2:
                # Ether oxygen found
                ether_oxygen_atoms.append(atom.GetIdx())
    
    # Check if ether bridges connect isoquinoline units
    for oxygen_idx in ether_oxygen_atoms:
        oxygen_atom = mol.GetAtomWithIdx(oxygen_idx)
        neighbor_carbons = [n for n in oxygen_atom.GetNeighbors() if n.GetAtomicNum() == 6]
        units_connected = set()
        for carbon_atom in neighbor_carbons:
            carbon_idx = carbon_atom.GetIdx()
            for i, unit in enumerate(isoquinoline_units):
                if carbon_idx in unit:
                    units_connected.add(i)
        if len(units_connected) >= 2:
            return True, "Contains two isoquinoline units connected via ether bridge(s)"
    
    # Check for methylenedioxy bridges (O-CH2-O)
    methylenedioxy_pattern = Chem.MolFromSmarts('COC')  # Simplified pattern
    methylenedioxy_matches = mol.GetSubstructMatches(methylenedioxy_pattern)
    for match in methylenedioxy_matches:
        # Check if methylenedioxy group connects two isoquinoline units
        units_connected = set()
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6:  # Carbon atom in CH2
                neighbors = atom.GetNeighbors()
                oxygen_neighbors = [n for n in neighbors if n.GetAtomicNum() == 8]
                for oxygen_atom in oxygen_neighbors:
                    carbon_neighbors = [n for n in oxygen_atom.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIdx() != atom_idx]
                    for carbon_neighbor in carbon_neighbors:
                        carbon_idx = carbon_neighbor.GetIdx()
                        for i, unit in enumerate(isoquinoline_units):
                            if carbon_idx in unit:
                                units_connected.add(i)
        if len(units_connected) >= 2:
            return True, "Contains two isoquinoline units connected via methylenedioxy bridge(s)"
    
    # Check for direct carbon-carbon bridges between isoquinoline units
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6:
            units_connected = set()
            for i, unit in enumerate(isoquinoline_units):
                if atom1.GetIdx() in unit or atom2.GetIdx() in unit:
                    units_connected.add(i)
            if len(units_connected) >= 2:
                return True, "Contains two isoquinoline units connected via direct carbon-carbon bridge(s)"
    
    return False, "No bridges connecting isoquinoline units found"


__metadata__ = {   
    'chemical_class': {   
        'id': None,
        'name': 'bisbenzylisoquinoline alkaloid',
        'definition': 'A type of benzylisoquinoline alkaloid whose structures are built up of two benzylisoquinoline units linked by ether bridges. Various structural patterns resulting from additional bridging between the two units by direct carbon-carbon bridging or by methylenedioxy groups are common.',
        'parents': ['benzylisoquinoline alkaloid']},
    'config': {   
        'llm_model_name': None,
        'f1_threshold': None,
        'max_attempts': None,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': None,
        'max_negative_in_prompt': None,
        'max_instances_in_prompt': None,
        'test_proportion': None},
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
    'accuracy': None}