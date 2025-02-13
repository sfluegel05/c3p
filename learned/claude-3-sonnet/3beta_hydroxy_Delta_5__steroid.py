"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
"""
Classifies: 3beta-hydroxy-Delta(5)-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy-Delta(5)-steroid based on its SMILES string.
    These compounds have a steroid core with a 3-beta hydroxyl group and a double bond between C5-C6.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3beta-hydroxy-Delta(5)-steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic steroid core pattern - four fused rings
    # More flexible pattern that matches the cyclopentanoperhydrophenanthrene core
    steroid_core = Chem.MolFromSmarts("C1C2CCC3CCCC4CCCC(C4)C3C2C1")
    
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Look for 3β-hydroxyl group with correct stereochemistry
    # Multiple patterns to catch different SMILES representations
    beta_hydroxyl_patterns = [
        Chem.MolFromSmarts("[C@@H](O)CC"), # Direct pattern
        Chem.MolFromSmarts("[C@H](CC)O"),  # Alternative representation
        Chem.MolFromSmarts("C[C@@H](O)C"),  # Another common form
    ]
    
    has_beta_hydroxyl = any(mol.HasSubstructMatch(pattern) 
                           for pattern in beta_hydroxyl_patterns 
                           if pattern is not None)
    
    if not has_beta_hydroxyl:
        return False, "No 3-beta-hydroxyl group found"

    # Check for double bond between C5-C6
    # Multiple patterns to catch different representations of the Δ5 double bond
    delta_5_patterns = [
        Chem.MolFromSmarts("C=CC1CCC"), # Basic pattern
        Chem.MolFromSmarts("C1CC=CC(C1)"), # Ring-based pattern
        Chem.MolFromSmarts("C=C1CCCC1"), # Alternative ring pattern
    ]
    
    has_delta_5 = any(mol.HasSubstructMatch(pattern) 
                     for pattern in delta_5_patterns 
                     if pattern is not None)
    
    if not has_delta_5:
        return False, "No double bond between positions 5 and 6"

    # Additional validation checks
    
    # Count rings (steroids should have 4 or more)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings for steroid structure"

    # Count carbons (steroids typically have 17+ carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 17:
        return False, "Too few carbons for a steroid structure"

    # Count oxygens (should have at least one for the hydroxyl group)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 1:
        return False, "No oxygen atoms found"

    # Verify sp2 carbons (should have at least 2 for the double bond)
    sp2_carbons = sum(1 for atom in mol.GetAtoms() 
                     if atom.GetAtomicNum() == 6 
                     and atom.GetHybridization() == Chem.HybridizationType.SP2)
    if sp2_carbons < 2:
        return False, "No carbon-carbon double bond found"

    return True, "Contains steroid core with 3-beta-hydroxyl group and double bond between C5-C6"