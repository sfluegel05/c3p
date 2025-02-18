"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
"""
Classifies: secondary alpha-hydroxy ketone (acyloins)
"""
from rdkit import Chem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone (acyloin) based on its SMILES string.
    A secondary alpha-hydroxy ketone has a carbonyl group (C=O) and a hydroxyl group (-OH) on adjacent carbons,
    with the alpha carbon (bearing the OH) having exactly one hydrogen and one organyl group (non-O, non-H).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Primary SMARTS pattern: [CH](O)(R)C=O where R is organyl (non-O, non-H)
    # This matches alpha carbon with exactly one H, one OH, one R group, and adjacent to ketone
    acyloin_pattern = Chem.MolFromSmarts('[CH]([OH])([!O;!H])C(=O)')
    if mol.HasSubstructMatch(acytoin_pattern):
        return True, "Alpha-hydroxy ketone with required substituents"
    
    # Secondary check for cases where SMARTS might miss due to stereochemistry or explicit H handling
    # Match any carbon with OH and adjacent to ketone, then verify substituents
    alt_pattern = Chem.MolFromSmarts('[C]([OH])C(=O)')
    matches = mol.GetSubstructMatches(alt_pattern)
    for match in matches:
        alpha_idx = match[0]
        alpha_atom = mol.GetAtomWithIdx(alpha_idx)
        
        # Count substituents: must have exactly one H and one R group (non-O, non-H)
        h_count = 0
        r_groups = 0
        for neighbor in alpha_atom.GetNeighbors():
            if neighbor.GetIdx() == match[1]:  # Skip ketone carbon
                continue
            if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() == 1:  # OH group
                continue
            if neighbor.GetAtomicNum() == 1:
                h_count += 1
            else:
                r_groups += 1
        
        if h_count >= 1 and r_groups >= 1:
            return True, "Verified alpha-hydroxy ketone with required substituents"
    
    return False, "No secondary alpha-hydroxy ketone structure detected"