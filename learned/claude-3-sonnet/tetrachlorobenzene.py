"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
"""
Classifies: tetrachlorobenzene
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule contains a tetrachlorobenzene substructure.
    Tetrachlorobenzene is a benzene ring with exactly 4 chlorine atoms attached.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains tetrachlorobenzene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Create a SMARTS pattern that matches a benzene ring with any substituents
    # and then we'll check for exactly 4 chlorines
    # The pattern uses recursive SMARTS to ensure we match a benzene ring 
    # where exactly 4 carbons have chlorine substituents
    pattern = """
        c1c(*)c(*)c(*)c(*)c(*)1    # benzene ring with any substituents
    """
    
    substructure = Chem.MolFromSmarts(pattern)
    if substructure is None:
        return False, "Invalid SMARTS pattern"
            
    matches = mol.GetSubstructMatches(substructure)
    
    for match in matches:
        # Get the benzene carbons from the match
        benzene_carbons = set(match[:6])  # First 6 atoms are the benzene ring carbons
        
        # Count chlorines directly attached to these carbons
        chlorine_count = 0
        for carbon_idx in benzene_carbons:
            atom = mol.GetAtomWithIdx(carbon_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 17:  # Chlorine
                    chlorine_count += 1
        
        # If we find exactly 4 chlorines on this benzene ring, it's a match
        if chlorine_count == 4:
            return True, "Found benzene ring with exactly 4 chlorine substituents"
    
    return False, "No benzene ring with exactly 4 chlorine substituents found"