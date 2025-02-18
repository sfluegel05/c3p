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
    on an alpha carbon (adjacent to a carboxy group) that connects both carboxy groups.

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

    # Get carboxy carbons from matches
    carboxy_carbons = [match[0] for match in carboxy_matches]
    carb1, carb2 = carboxy_carbons

    # Check connectivity between carboxy groups through the hydroxyl-bearing carbon
    hydroxyl_found = False
    for carboxy_carbon in [carb1, carb2]:
        # Get alpha carbons (direct neighbors to carboxy carbon)
        alpha_carbons = [n.GetIdx() for n in mol.GetAtomWithIdx(carboxy_carbon).GetNeighbors() if n.GetAtomicNum() == 6]
        
        for alpha in alpha_carbons:
            # Check if alpha carbon has hydroxyl (-OH)
            alpha_atom = mol.GetAtomWithIdx(alpha)
            hydroxyls = [n for n in alpha_atom.GetNeighbors() 
                        if n.GetAtomicNum() == 8 and n.GetTotalNumHs() >= 1]
            
            if hydroxyls:
                # Check if this alpha carbon is on the path between both carboxy groups
                path = Chem.GetShortestPath(mol, alpha, carb2 if carboxy_carbon == carb1 else carb1)
                if path:
                    hydroxyl_found = True
                    break
        if hydroxyl_found:
            break

    if not hydroxyl_found:
        return False, "No hydroxyl on alpha carbon connecting both carboxy groups"

    return True, "Two carboxy groups with hydroxyl on connecting alpha carbon"