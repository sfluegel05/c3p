"""
Classifies: CHEBI:15734 primary alcohol
"""
"""
Classifies: primary alcohol
A primary alcohol has a hydroxyl group attached to a saturated carbon that has either:
- three hydrogens and one other atom, or
- two hydrogens and one carbon atom
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_primary_alcohol(smiles: str):
    """
    Determines if a molecule contains a primary alcohol group.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains a primary alcohol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Pattern 1: CH2OH connected to carbon (most common case)
    # [CD1,CD2,CD3,CD4] means a carbon with 1-4 connections total
    # [$([CH2][OH])] means CH2 group connected to OH
    pattern1 = Chem.MolFromSmarts('[$([CH2][OH])]([CD1,CD2,CD3,CD4])')
    
    # Pattern 2: Isolated CH2OH group (terminal)
    pattern2 = Chem.MolFromSmarts('[CH2][OH]')
    
    # Pattern 3: CH3OH type (methanol type)
    pattern3 = Chem.MolFromSmarts('[CH3][OH]')
    
    # Check each pattern
    matches1 = mol.GetSubstructMatches(pattern1)
    matches2 = mol.GetSubstructMatches(pattern2)
    matches3 = mol.GetSubstructMatches(pattern3)
    
    # For each match, verify it's a true primary alcohol
    for match in matches1 + matches2 + matches3:
        # Get the carbon atom
        carbon = mol.GetAtomWithIdx(match[0])
        
        # Skip if carbon is not sp3 (must be saturated)
        if carbon.GetHybridization() != Chem.HybridizationType.SP3:
            continue
            
        # Count non-hydrogen connections
        heavy_neighbors = [n for n in carbon.GetNeighbors() 
                         if n.GetAtomicNum() != 1]
        
        # For CH3OH type
        if len(heavy_neighbors) == 1 and carbon.GetTotalNumHs() == 3:
            if any(n.GetAtomicNum() == 8 for n in heavy_neighbors):
                return True, "Found methanol-type primary alcohol (CH3-OH)"
                
        # For RCH2OH type
        elif len(heavy_neighbors) == 2 and carbon.GetTotalNumHs() == 2:
            # One must be oxygen, other must be carbon
            if (sum(1 for n in heavy_neighbors if n.GetAtomicNum() == 8) == 1 and
                sum(1 for n in heavy_neighbors if n.GetAtomicNum() == 6) == 1):
                return True, "Found RCH2-OH type primary alcohol"
    
    return False, "No primary alcohol groups found"