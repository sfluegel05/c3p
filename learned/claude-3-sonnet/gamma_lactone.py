"""
Classifies: CHEBI:37581 gamma-lactone
"""
"""
Classifies: gamma-lactone (5-membered lactone ring)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule contains a gamma-lactone (5-membered lactone ring).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains a gamma-lactone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic gamma-lactone pattern:
    # - [#8] represents oxygen
    # - [#6] represents carbon
    # - ;!@1 ensures it's part of a ring
    # The pattern matches any 5-membered ring with an ester group
    patterns = [
        # Saturated gamma-lactone
        "[#8;R1]1-[#6;R1]-[#6;R1]-[#6;R1]-[#6;R1](=[O;R1])1",
        
        # Unsaturated gamma-lactone (one double bond)
        "[#8;R1]1-[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1](=[O;R1])1",
        
        # Alternative unsaturated pattern
        "[#8;R1]1-[#6;R1]=[#6;R1]-[#6;R1]-[#6;R1](=[O;R1])1",
        
        # More general pattern that can catch other variants
        "[#8;R1]1[#6;R1][#6;R1][#6;R1][#6;R1](=[O;R1])1"
    ]
    
    all_matches = 0
    for pattern in patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt is not None:
            matches = mol.GetSubstructMatches(patt)
            all_matches += len(matches)
    
    if all_matches > 0:
        # Remove potential duplicate matches
        # Use a more general pattern to count unique lactone rings
        basic_pattern = Chem.MolFromSmarts("[#8;R1]1[#6][#6][#6][#6](=[O;R1])1")
        unique_matches = len(mol.GetSubstructMatches(basic_pattern))
        
        if unique_matches == 1:
            return True, "Contains a gamma-lactone (5-membered lactone ring)"
        else:
            return True, f"Contains {unique_matches} gamma-lactone rings"
            
    return False, "No gamma-lactone substructure found"