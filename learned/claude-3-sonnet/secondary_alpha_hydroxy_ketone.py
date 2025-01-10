"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
"""
Classifies: secondary alpha-hydroxy ketone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone based on its SMILES string.
    A secondary alpha-hydroxy ketone has a carbonyl group adjacent to a carbon bearing
    one hydroxy group, one hydrogen, and one organyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens for proper matching
    mol = Chem.AddHs(mol)
    
    # SMARTS pattern for secondary alpha-hydroxy ketone:
    # [CH1] - carbon with exactly one H
    # [CX3](=O) - carbonyl group
    # [OX2H1] - hydroxy group
    # The central carbon must have exactly one H and be connected to a carbon 
    # (representing the organyl group)
    pattern = Chem.MolFromSmarts("[CX3](=[OX1])[CH1]([OX2H1])[#6]")
    
    if not mol.HasSubstructMatch(pattern):
        return False, "No secondary alpha-hydroxy ketone pattern found"
    
    matches = mol.GetSubstructMatches(pattern)
    
    # Additional validation to ensure correct pattern
    for match in matches:
        carbonyl_c = mol.GetAtomWithIdx(match[0])
        alpha_c = mol.GetAtomWithIdx(match[1])
        hydroxy_o = mol.GetAtomWithIdx(match[2])
        
        # Verify the alpha carbon has exactly one hydrogen
        if alpha_c.GetTotalNumHs() != 1:
            continue
            
        # Verify the carbonyl carbon is sp2 hybridized
        if carbonyl_c.GetHybridization() != Chem.HybridizationType.SP2:
            continue
            
        # Verify hydroxy oxygen has one hydrogen
        if hydroxy_o.GetTotalNumHs() != 1:
            continue
            
        return True, "Contains carbonyl group adjacent to carbon with OH, H, and organyl group"
        
    return False, "Structure does not match secondary alpha-hydroxy ketone pattern"