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
    # Split into components to handle salts/mixtures
    mols = Chem.SmilesMolSupplierFromText(smiles, sanitize=False, delimiter='.')
    mol = None
    for m in mols:
        if m is not None:
            if mol is None:
                mol = m
            else:
                return False, "Multiple molecules present"
    if mol is None:
        return False, "Invalid SMILES"
    
    try:
        Chem.SanitizeMol(mol)
    except:
        return False, "Sanitization failed"

    # Find carboxylic acid groups (-COOH)
    carboxy_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    carboxy_matches = mol.GetSubstructMatches(carboxy_pattern)
    if len(carboxy_matches) != 2:
        return False, f"Found {len(carboxy_matches)} carboxy groups (needs exactly 2)"

    # Get carboxy carbons from matches
    carb1 = carboxy_matches[0][0]
    carb2 = carboxy_matches[1][0]

    # Find all alpha carbons (adjacent to either carboxy carbon)
    alpha_carbons = set()
    for carb in [carb1, carb2]:
        neighbors = mol.GetAtomWithIdx(carb).GetNeighbors()
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() == 6:
                alpha_carbons.add(neighbor.GetIdx())

    # Check for hydroxyl groups on alpha carbons
    hydroxyl_found = False
    for alpha in alpha_carbons:
        alpha_atom = mol.GetAtomWithIdx(alpha)
        # Check for -OH attached to alpha carbon
        for bond in alpha_atom.GetBonds():
            if bond.GetBondType() == Chem.BondType.SINGLE:
                neighbor = bond.GetOtherAtom(alpha_atom)
                if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() >= 1:
                    hydroxyl_found = True
                    break
        if hydroxyl_found:
            break

    if not hydroxyl_found:
        return False, "No hydroxyl on alpha carbon"

    # Check if the hydroxyl-bearing alpha connects both carboxy groups
    try:
        path = Chem.GetShortestPath(mol, carb1, carb2)
    except:
        return False, "No path between carboxy groups"
    
    if alpha not in path:
        return False, "Hydroxyl not on connecting path"

    return True, "Two carboxy groups with hydroxyl on connecting alpha carbon"