"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
"""
Classifies: CHEBI:39424 fatty acid methyl ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester based on its SMILES string.
    A fatty acid methyl ester is the carboxylic ester formed by the condensation of a fatty acid with methanol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid methyl ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for exactly one methyl ester group (O=C-O-CH3)
    methyl_ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CH3]")
    methyl_ester_matches = mol.GetSubstructMatches(methyl_ester_pattern)
    if len(methyl_ester_matches) != 1:
        return False, f"Found {len(methyl_ester_matches)} methyl ester groups, need exactly 1"

    # Ensure no other ester groups (any other O-C=O-O- groups)
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Check for absence of carboxylic acid groups (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Carboxylic acid group present"

    # Extract the R group (the part attached to the carbonyl carbon)
    # Get the carbonyl carbon from the methyl ester match
    carbonyl_carbon = methyl_ester_matches[0][0]
    # Find the R group atom (connected to carbonyl carbon, not part of ester oxygen or carbonyl oxygen)
    neighbors = mol.GetAtomWithIdx(carbonyl_carbon).GetNeighbors()
    r_group_atom = None
    for neighbor in neighbors:
        if neighbor.GetIdx() not in {methyl_ester_matches[0][1], methyl_ester_matches[0][2]}:
            r_group_atom = neighbor
            break
    if r_group_atom is None:
        return False, "No R group attached to carbonyl carbon"

    # Traverse R group to count carbons and check for other heteroatoms (allowed in substituents)
    # Count carbons in the R group (excluding the ester part)
    visited = set()
    stack = [r_group_atom]
    r_group_carbons = 0

    while stack:
        atom = stack.pop()
        if atom.GetIdx() in visited:
            continue
        visited.add(atom.GetIdx())
        if atom.GetAtomicNum() == 6:
            r_group_carbons += 1
        # Add all neighbors except those in the ester group
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in {carbonyl_carbon, methyl_ester_matches[0][1], methyl_ester_matches[0][2], methyl_ester_matches[0][3]}:
                stack.append(neighbor)

    # Check minimum carbon count in R group (at least 4 carbons)
    if r_group_carbons < 4:
        return False, f"R group has only {r_group_carbons} carbons, need at least 4"

    return True, "Contains a methyl ester group with a hydrocarbon chain of sufficient length"