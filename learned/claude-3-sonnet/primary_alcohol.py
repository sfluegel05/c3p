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
    
    # Add explicit hydrogens
    mol = Chem.AddHs(mol)
    
    # Find all OH groups
    oh_pattern = Chem.MolFromSmarts("[OH]")
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    
    if not oh_matches:
        return False, "No hydroxyl groups found"
    
    # For each OH group, check if it's a primary alcohol
    for match in oh_matches:
        oh_oxygen = mol.GetAtomWithIdx(match[0])
        
        # Get the carbon attached to the OH
        for neighbor in oh_oxygen.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon
                # Check if carbon is sp3 (saturated)
                if neighbor.GetHybridization() != Chem.HybridizationType.SP3:
                    continue
                
                # Count hydrogens and non-hydrogen neighbors on the carbon
                h_count = neighbor.GetTotalNumHs()
                heavy_neighbors = [n for n in neighbor.GetNeighbors() 
                                 if n.GetAtomicNum() != 1]  # Exclude hydrogens
                
                # Case 1: CH3-OH type (3 H, 1 O)
                if h_count == 3 and len(heavy_neighbors) == 1:
                    return True, "Found CH3-OH type primary alcohol"
                
                # Case 2: RCH2-OH type (2 H, 1 C, 1 O)
                if h_count == 2 and len(heavy_neighbors) == 2:
                    # One neighbor must be oxygen, other must be carbon
                    non_oh_neighbors = [n for n in heavy_neighbors if n.GetIdx() != oh_oxygen.GetIdx()]
                    if len(non_oh_neighbors) == 1 and non_oh_neighbors[0].GetAtomicNum() == 6:
                        return True, "Found RCH2-OH type primary alcohol"
    
    return False, "No primary alcohol groups found"