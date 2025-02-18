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
    with the alpha carbon (bearing the OH) having one hydrogen and one organyl group.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Define SMARTS pattern for the alpha-hydroxy ketone structure
    # The pattern matches a carbon adjacent to a ketone, bearing an OH group,
    # one hydrogen, and one organyl group (non-O, non-H)
    acyloin_pattern = Chem.MolFromSmarts('[CH]([OH])([!O;!H])C(=O)')
    if mol.HasSubstructMatch(acyloin_pattern):
        return True, "Contains alpha-hydroxy ketone with required substituents"
    
    # Check for cases where the alpha carbon might have multiple organyl groups (like benzoin)
    # Alternative pattern without hydrogen check
    alt_pattern = Chem.MolFromSmarts('[C]([OH])([!O;!H])C(=O)')
    if mol.HasSubstructMatch(alt_pattern):
        # Verify the alpha carbon has at least one organyl group and the hydroxyl
        # Manually check each match to ensure at least one R group
        matches = mol.GetSubstructMatches(alt_pattern)
        for match in matches:
            alpha_carbon = match[0]
            # Count non-O, non-H substituents (excluding the ketone carbon)
            substituents = []
            for neighbor in alpha_carbon.GetNeighbors():
                if neighbor.GetIdx() == match[1]:  # ketone carbon is part of the match
                    continue
                if neighbor.GetSymbol() != 'O' and neighbor.GetSymbol() != 'H':
                    substituents.append(neighbor)
            if len(substituents) >= 1:
                return True, "Contains alpha-hydroxy ketone with organyl group"
    
    return False, "No secondary alpha-hydroxy ketone structure found"