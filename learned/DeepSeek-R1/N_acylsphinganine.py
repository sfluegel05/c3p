"""
Classifies: CHEBI:31488 N-acylsphinganine
"""
"""
Classifies: CHEBI:89998 N-acylsphinganine
"""
from rdkit import Chem

def is_N_acylsphinganine(smiles: str):
    """
    Determines if a molecule is an N-acylsphinganine based on its SMILES string.
    An N-acylsphinganine consists of sphinganine with a fatty acyl group attached to the amino group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylsphinganine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define core structure pattern using SMARTS:
    # [NH] connected to carbonyl, followed by two adjacent carbons with hydroxyl groups
    # Stereochemistry-aware pattern for sphinganine backbone:
    # N-C(=O)-C(CO)-C(OH)- (with specific stereochemistry)
    pattern = Chem.MolFromSmarts('[NH1]C(=O)[C@@H](CO)[C@H](O)')
    if not mol.HasSubstructMatch(pattern):
        return False, "Missing core sphinganine backbone with acyl group"

    # Verify the fatty acyl chain has at least 2 carbons (acetyl is minimal example)
    # Find the amide bond atom indices
    amide_matches = mol.GetSubstructMatches(Chem.MolFromSmarts('[NH1]C(=O)'))
    if not amide_matches:
        return False, "No amide group found"
    
    # Get the carbonyl carbon (index 1 in the match [N, C=O])
    carbonyl_carbon = amide_matches[0][1]
    neighbor_atoms = mol.GetAtomWithIdx(carbonyl_carbon).GetNeighbors()
    
    # The acyl chain is the neighbor of carbonyl carbon not part of the amide
    acyl_chain_start = [n for n in neighbor_atoms if n.GetIdx() != amide_matches[0][0]][0]
    
    # Traverse the acyl chain to count carbons
    chain_length = 0
    visited = set()
    stack = [(acyl_chain_start, 0)]
    while stack:
        atom, depth = stack.pop()
        if atom.GetAtomicNum() == 6 and atom.GetIdx() not in visited:
            visited.add(atom.GetIdx())
            chain_length += 1
            # Follow non-carbonyl, non-amide connections
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in visited and neighbor.GetIdx() != carbonyl_carbon:
                    stack.append((neighbor, depth + 1))
    
    if chain_length < 2:  # At least acetyl (2 carbons including carbonyl)
        return False, f"Acyl chain too short (length {chain_length})"

    return True, "Contains N-acylsphinganine structure with proper acyl chain"