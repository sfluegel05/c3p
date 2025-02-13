"""
Classifies: CHEBI:22307 aldoxime
"""
"""
Classifies: CHEBI:33566 aldoxime
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime based on its SMILES string.
    Aldoximes are compounds with the general structure H-CR=N-OH, where R is any carbon group.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an aldoxime, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Main pattern for aldoxime group: CH=N-OH where C has exactly one H
    # Use recursive SMARTS to exclude C=N-OH attached to C=O (glyoximes)
    pattern = "[CH1]([!$(C=O)])(=N[OH1])"
    
    smarts = Chem.MolFromSmarts(pattern)
    if not mol.HasSubstructMatch(smarts):
        return False, "No aldoxime group (H-CR=N-OH) found"
        
    matches = mol.GetSubstructMatches(smarts)
    
    for match in matches:
        carbon_idx = match[0]  # First atom in pattern is the carbon
        carbon = mol.GetAtomWithIdx(carbon_idx)
        
        # Verify carbon has exactly one hydrogen
        if carbon.GetTotalNumHs() != 1:
            continue
            
        # Get all neighbors of the carbon
        neighbors = carbon.GetNeighbors()
        
        # Count non-nitrogen neighbors (should be only one)
        non_n_neighbors = [n for n in neighbors if n.GetAtomicNum() != 7]
        if len(non_n_neighbors) != 1:
            continue
            
        # Verify the nitrogen neighbor is indeed part of the oxime
        n_neighbors = [n for n in neighbors if n.GetAtomicNum() == 7]
        if not n_neighbors:
            continue
            
        nitrogen = n_neighbors[0]
        o_neighbors = [n for n in nitrogen.GetNeighbors() if n.GetAtomicNum() == 8]
        if not o_neighbors:
            continue
            
        # Check that the oxygen is -OH (has one H)
        oxygen = o_neighbors[0]
        if oxygen.GetTotalNumHs() != 1:
            continue
            
        # If we get here, we have a valid aldoxime group
        return True, "Contains aldoxime group (H-CR=N-OH)"
        
    return False, "No aldoxime group (H-CR=N-OH) found"