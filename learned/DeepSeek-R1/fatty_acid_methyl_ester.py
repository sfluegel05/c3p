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

    # Check for at least one methyl ester group (O=C-O-CH3)
    methyl_ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CH3]")
    methyl_ester_matches = mol.GetSubstructMatches(methyl_ester_pattern)
    if not methyl_ester_matches:
        return False, "No methyl ester groups found"

    # Ensure all ester groups are methyl esters
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    for ester in ester_matches:
        # Check if the ester is a methyl ester
        carbonyl_carbon = ester[0]
        oxygen = ester[1]
        # Find the atom attached to the ester oxygen (should be CH3)
        oxygen_atom = mol.GetAtomWithIdx(oxygen)
        neighbors = oxygen_atom.GetNeighbors()
        # The neighbor should be a carbon with three hydrogens (CH3)
        is_methyl = False
        for neighbor in neighbors:
            if neighbor.GetIdx() == carbonyl_carbon:
                continue
            if neighbor.GetAtomicNum() == 6 and neighbor.GetDegree() == 1:  # CH3
                is_methyl = True
        if not is_methyl:
            return False, "Non-methyl ester group present"

    # Check for absence of carboxylic acid groups (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Carboxylic acid group present"

    # Check each methyl ester's R group (hydrocarbon chain with >=4 carbons)
    for methyl_ester in methyl_ester_matches:
        carbonyl_carbon = methyl_ester[0]
        # Get the R group atom attached to the carbonyl carbon (not part of the ester)
        r_group_atom = None
        for neighbor in mol.GetAtomWithIdx(carbonyl_carbon).GetNeighbors():
            if neighbor.GetIdx() not in {methyl_ester[1], methyl_ester[2], methyl_ester[3]}:
                r_group_atom = neighbor
                break
        if r_group_atom is None:
            return False, "No R group attached to carbonyl carbon"

        # Traverse R group to check for hydrocarbons and count carbons
        visited = set()
        stack = [r_group_atom]
        r_group_carbons = 0
        valid_hydrocarbon = True

        while stack:
            atom = stack.pop()
            if atom.GetIdx() in visited:
                continue
            visited.add(atom.GetIdx())
            # Check if atom is carbon
            if atom.GetAtomicNum() != 6:
                valid_hydrocarbon = False
                break
            r_group_carbons += 1
            # Add all neighbors except those in the ester group
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in visited and neighbor.GetIdx() != carbonyl_carbon:
                    stack.append(neighbor)

        if not valid_hydrocarbon:
            return False, "R group contains non-carbon atoms"
        if r_group_carbons < 4:
            return False, f"R group has only {r_group_carbons} carbons, need at least 4"

    return True, "Contains methyl ester(s) with hydrocarbon chain(s) of sufficient length"