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

    # Patterns for gamma-lactone detection
    patterns = [
        # Basic gamma-lactone pattern - most general
        # Matches any 5-membered ring with O-C(=O) motif
        "[#8X2r5]1-[#6r5]-[#6r5]-[#6r5]-[#6X3r5](=[O])1",
        
        # Alternative pattern with any bond types
        "[#8X2r5]1~[#6r5]~[#6r5]~[#6r5]~[#6X3r5](=[O])1",
        
        # Pattern for fused systems
        "[#8X2r5](@[#6r5])(@[#6X3r5](=[O]))@[#6r5]@[#6r5]",
        
        # Pattern for substituted variants
        "[#8X2r5]1[#6r5][#6r5][#6r5][#6X3r5](=O)1"
    ]
    
    for pattern in patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt is not None:
            if mol.HasSubstructMatch(patt):
                # Count unique matches using the most general pattern
                matches = len(mol.GetSubstructMatches(patt))
                if matches == 1:
                    return True, "Contains a gamma-lactone (5-membered lactone ring)"
                else:
                    return True, f"Contains {matches} gamma-lactone rings"
    
    # Additional check for special cases using a more relaxed pattern
    backup_pattern = "[#8r5]1~[#6]~[#6]~[#6]~[#6](=[O])1"
    patt = Chem.MolFromSmarts(backup_pattern)
    if patt is not None and mol.HasSubstructMatch(patt):
        return True, "Contains a gamma-lactone variant"
            
    return False, "No gamma-lactone substructure found"