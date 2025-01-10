"""
Classifies: CHEBI:35359 carboxamidine
"""
"""
Classifies: CHEBI:51381 carboxamidine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule contains a carboxamidine group (RC(=NR)NR2).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains carboxamidine group, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Basic carboxamidine pattern: carbon with double bond to N and single bond to N
    # [C;!$(C=O)] ensures we don't match amides
    # [$([NX2]=[CX3]);!$(N=C=N)] ensures we match C=N but not N=C=N
    # [NX3;!$(NC=O)] ensures we match C-N but not amide N
    pattern = Chem.MolFromSmarts("[C;!$(C=O)]([#6,#1,N])([NX3;!$(NC=O)])=[$([NX2]=[CX3]);!$(N=C=N)]")
    
    # Alternative pattern for cases where N order might be swapped
    pattern2 = Chem.MolFromSmarts("[C;!$(C=O)]([#6,#1,N])(=[$([NX2]=[CX3]);!$(N=C=N)])[NX3;!$(NC=O)]")
    
    # Check for matches
    matches = mol.GetSubstructMatches(pattern)
    matches2 = mol.GetSubstructMatches(pattern2)
    all_matches = matches + matches2
    
    if not all_matches:
        return False, "No carboxamidine group found"
    
    # Additional check to exclude guanidine (N=C(N)N) pattern
    guanidine_pattern = Chem.MolFromSmarts("[NX3][CX3](=[NX2])[NX3]")
    guanidine_matches = mol.GetSubstructMatches(guanidine_pattern)
    
    # If all matches are part of guanidine groups, it's not a carboxamidine
    non_guanidine_matches = False
    for match in all_matches:
        carbon_idx = match[0]
        is_guanidine = False
        for g_match in guanidine_matches:
            if carbon_idx in g_match:
                is_guanidine = True
                break
        if not is_guanidine:
            non_guanidine_matches = True
            break
    
    if not non_guanidine_matches and guanidine_matches:
        return False, "Contains guanidine but not carboxamidine group"
        
    # Count number of unique carboxamidine groups
    unique_carbons = set(match[0] for match in all_matches)
    
    return True, f"Contains {len(unique_carbons)} carboxamidine group(s)"