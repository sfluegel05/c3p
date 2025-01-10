"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
# bisbenzylisoquinoline_alkaloid.py
"""
Classifies: bisbenzylisoquinoline alkaloid
Definition: A type of benzylisoquinoline alkaloid whose structures are built up of two benzylisoquinoline units linked by ether bridges. Various structural patterns resulting from additional bridging between the two units by direct carbon-carbon bridging or by methylenedioxy groups are common.
"""

from rdkit import Chem

def is_bisbenzylisoquinoline_alkaloid(smiles: str):
    """
    Determines if a molecule is a bisbenzylisoquinoline alkaloid based on its SMILES string.
    A bisbenzylisoquinoline alkaloid consists of two benzylisoquinoline units linked by ether bridges,
    methylenedioxy groups, or direct carbon-carbon bonds.
    
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
    
    # Define benzylisoquinoline unit pattern (allows for substitutions and saturation)
    benzylisoquinoline_smarts = '''
    [
        # Isoquinoline or tetrahydroisoquinoline core
        $([
            $([nH]1ccccc1),                # Tetrahydroisoquinoline
            $([n]1ccccc1)                  # Isoquinoline
        ]),
        # Connected to benzyl group
        $(C[cH1,cH0][cR]),                # Benzyl group connected to aromatic ring
        # Allow for substitutions on the rings
        $(*)
    ]
    '''
    benzylisoquinoline_pattern = Chem.MolFromSmarts(benzylisoquinoline_smarts)
    if benzylisoquinoline_pattern is None:
        return False, "Invalid benzylisoquinoline SMARTS pattern"
    
    # Find benzylisoquinoline units
    benzylisoquinoline_matches = mol.GetSubstructMatches(benzylisoquinoline_pattern)
    num_units = len(benzylisoquinoline_matches)
    if num_units < 2:
        return False, f"Found {num_units} benzylisoquinoline unit(s), need at least 2"
    
    # Get atom indices of benzylisoquinoline units
    unit_atoms = []
    for match in benzylisoquinoline_matches:
        unit_atoms.append(set(match))
    
    # Identify bridges between units
    units_connected = set()
    for bond in mol.GetBonds():
        atom1_idx = bond.GetBeginAtomIdx()
        atom2_idx = bond.GetEndAtomIdx()
        units_involved = []
        for i, unit in enumerate(unit_atoms):
            if atom1_idx in unit or atom2_idx in unit:
                units_involved.append(i)
        if len(units_involved) == 2:
            # Check if bridge is ether (C-O-C), methylenedioxy (O-CH2-O), or C-C bond
            atom1 = mol.GetAtomWithIdx(atom1_idx)
            atom2 = mol.GetAtomWithIdx(atom2_idx)
            if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                if atom1.GetAtomicNum() == 8 or atom2.GetAtomicNum() == 8:
                    # Oxygen bridge
                    units_connected.update(units_involved)
                elif atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6:
                    # Carbon bridge
                    units_connected.update(units_involved)
            elif bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                # Possible methylenedioxy bridge (but simplified here)
                units_connected.update(units_involved)
    
    if len(units_connected) < 2:
        return False, "No bridges connecting benzylisoquinoline units found"
    
    return True, "Contains two benzylisoquinoline units connected via bridge(s)"

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