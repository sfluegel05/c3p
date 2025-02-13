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

    # Create SMARTS patterns for different tetrachlorobenzene arrangements
    # Match a benzene ring where exactly 4 carbons have chlorine substituents
    patterns = [
        # Generic pattern for any tetrachlorobenzene
        "c1(-[Cl])c(-[Cl])c(-[Cl])c(-[Cl])cc1",
        "c1(-[Cl])c(-[Cl])c(-[Cl])cc(-[Cl])1",
        "c1(-[Cl])c(-[Cl])cc(-[Cl])c(-[Cl])1",
    ]

    for pattern in patterns:
        substructure = Chem.MolFromSmarts(pattern)
        if substructure is None:
            continue
            
        matches = mol.GetSubstructMatches(substructure)
        if matches:
            # For each match, verify it's a valid tetrachlorobenzene
            for match in matches:
                # Get the benzene carbons from the match
                benzene_carbons = set(match[:6])  # First 6 atoms are the benzene carbons
                
                # Count chlorines attached to these specific carbons
                chlorine_count = 0
                for carbon_idx in benzene_carbons:
                    atom = mol.GetAtomWithIdx(carbon_idx)
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetAtomicNum() == 17:  # Chlorine
                            chlorine_count += 1
                
                # Verify exactly 4 chlorines are attached to the benzene ring
                if chlorine_count == 4:
                    return True, "Found valid tetrachlorobenzene structure"
    
    return False, "No tetrachlorobenzene substructure found"