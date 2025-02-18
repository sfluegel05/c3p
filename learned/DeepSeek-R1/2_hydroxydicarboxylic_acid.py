"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
"""
Classifies: 2-hydroxydicarboxylic acid (CHEBI:144810)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid based on its SMILES string.
    A 2-hydroxydicarboxylic acid has two carboxy groups (-COOH) and a hydroxyl group (-OH)
    on the alpha carbon (adjacent to at least one carboxy group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule matches the criteria, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Find carboxylic acid groups (-COOH)
    carboxy_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    carboxy_matches = mol.GetSubstructMatches(carboxy_pattern)
    if len(carboxy_matches) != 2:
        return False, f"Found {len(carboxy_matches)} carboxy groups (needs exactly 2)"

    # Check for hydroxyl group on alpha carbon (adjacent to any carboxy group)
    alpha_hydroxyl = False
    for carboxy_match in carboxy_matches:
        carboxy_carbon = carboxy_match[0]  # Carbon in the carboxy group (C=O)
        # Get neighboring atoms to the carboxy carbon
        neighbors = mol.GetAtomWithIdx(carboxy_carbon).GetNeighbors()
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() == 6:  # Adjacent carbon (alpha position)
                # Check if this alpha carbon has a hydroxyl group
                hydroxyls = [a for a in neighbor.GetNeighbors() if a.GetAtomicNum() == 8 and a.GetTotalNumHs() >= 1]
                if hydroxyls:
                    alpha_hydroxyl = True
                    break
        if alpha_hydroxyl:
            break

    if not alpha_hydroxyl:
        return False, "No hydroxyl group on alpha carbon adjacent to carboxy group"

    return True, "Two carboxy groups with hydroxyl on alpha carbon"