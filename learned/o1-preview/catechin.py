"""
Classifies: CHEBI:23053 catechin
"""
"""
Classifies: catechin
"""
from rdkit import Chem

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.
    A catechin is a hydroxyflavan that has a flavan-3-ol skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catechin, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the catechin core structure
    # The pattern includes:
    # - A flavan-3-ol skeleton
    # - Optional hydroxyl groups at positions 5 and 7 on ring A
    # - Optional hydroxyl groups at positions 3', 4', and 5' on ring B
    # - Hydroxyl group at position 3 on ring C
    # - Ring connectivity matching the catechin structure
    catechin_smarts = """
    [#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1
    -[#6]2(-[#8])-[#6](-[#8])-[#6]3:[#6]:[#6]:[#6]:[#6]:[#6]:3
    -[#8]-[#6]2-[#6](-[#8])-[#6](-[#1])-[#6]1
    """

    # Remove whitespace and newlines from SMARTS pattern
    catechin_smarts = ''.join(catechin_smarts.split())
    
    pattern = Chem.MolFromSmarts(catechin_smarts)
    if pattern is None:
        return False, "Invalid SMARTS pattern for catechin"

    # Check if the molecule matches the catechin pattern
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "Catechin core structure not found"

    # Now check for hydroxyl groups at specific positions
    # We will check for at least 4 phenolic hydroxyl groups
    num_pheno_oh = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # Oxygen atom
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 1:
                neighbor = neighbors[0]
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIsAromatic():
                    num_pheno_oh += 1
    if num_pheno_oh < 4:
        return False, f"Not enough phenolic hydroxyl groups ({num_pheno_oh} found), catechins typically have multiple phenolic hydroxyl groups"

    # Check for hydroxyl group at position 3 (attached to saturated carbon)
    has_3_oh = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # Oxygen atom
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 1:
                neighbor = neighbors[0]
                if neighbor.GetAtomicNum() == 6 and neighbor.GetDegree() == 4 and not neighbor.GetIsAromatic():
                    has_3_oh = True
                    break
    if not has_3_oh:
        return False, "Hydroxyl group at position 3 not found"

    # Check for absence of additional rings or large substitutions
    ri = mol.GetRingInfo()
    num_rings = ri.NumRings()
    if num_rings > 3:
        return False, f"Molecule has additional rings ({num_rings} rings found), which is unusual for catechins"

    return True, "Molecule matches catechin structure with necessary hydroxyl groups"