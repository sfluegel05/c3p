"""
Classifies: CHEBI:35681 secondary alcohol
"""
"""
Classifies: secondary alcohol
A compound with a hydroxy group attached to a saturated carbon atom 
which has two other carbon atoms attached to it.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_secondary_alcohol(smiles: str):
    """
    Determines if a molecule contains a secondary alcohol group.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains secondary alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Add explicit hydrogens to properly count connections
    mol = Chem.AddHs(mol)
    
    # Find all hydroxyl groups
    oh_pattern = Chem.MolFromSmarts("[OH1]")
    if not mol.HasSubstructMatch(oh_pattern):
        return False, "No hydroxyl group found"
    
    # Look for secondary alcohol pattern:
    # [C;X4] means sp3 carbon (exactly 4 connections)
    # [#6] means any carbon
    # The carbon bearing OH must have exactly 2 other carbons attached
    sec_alcohol_pattern = Chem.MolFromSmarts("[#6]-[C;X4]([OH1])-[#6]")
    matches = mol.GetSubstructMatches(sec_alcohol_pattern)
    
    if not matches:
        return False, "No secondary alcohol pattern found"
        
    # For each potential secondary alcohol, verify it's really a secondary alcohol
    for match in matches:
        central_carbon_idx = match[1]  # Index of carbon bearing OH
        central_carbon = mol.GetAtomWithIdx(central_carbon_idx)
        
        # Count number of carbon neighbors
        carbon_neighbors = sum(1 for neighbor in central_carbon.GetNeighbors() 
                             if neighbor.GetAtomicNum() == 6)
        
        if carbon_neighbors == 2:
            # Verify the carbon is sp3 (saturated)
            if central_carbon.GetHybridization() == Chem.HybridizationType.SP3:
                return True, "Contains hydroxyl group attached to saturated carbon with exactly two carbon neighbors"
    
    return False, "No valid secondary alcohol found - hydroxyl groups present but not meeting secondary alcohol criteria"