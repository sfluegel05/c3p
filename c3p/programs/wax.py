"""
Classifies: CHEBI:73702 wax
"""
"""
Classifies: Wax (long-chain ester)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_wax(smiles: str):
    """
    Determines if a molecule is a wax based on its SMILES string.
    Waxes are esters of long-chain alcohols and long-chain fatty acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for exactly one ester group
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Get ester oxygen and carbonyl carbon indices
    oxygen_idx = ester_matches[0][0]
    carbonyl_c_idx = ester_matches[0][1]

    # Identify alcohol component (connected to ester oxygen)
    oxygen = mol.GetAtomWithIdx(oxygen_idx)
    alcohol_start = None
    for neighbor in oxygen.GetNeighbors():
        if neighbor.GetIdx() != carbonyl_c_idx:
            alcohol_start = neighbor.GetIdx()
            break
    if alcohol_start is None:
        return False, "No alcohol component found"

    # Identify acid component (connected to carbonyl carbon)
    carbonyl_c = mol.GetAtomWithIdx(carbonyl_c_idx)
    acid_start = None
    for neighbor in carbonyl_c.GetNeighbors():
        if neighbor.GetIdx() == oxygen_idx:
            continue
        bond = mol.GetBondBetweenAtoms(carbonyl_c_idx, neighbor.GetIdx())
        if bond.GetBondType() == Chem.BondType.DOUBLE and neighbor.GetAtomicNum() == 8:
            continue  # Skip carbonyl oxygen
        acid_start = neighbor.GetIdx()
        break
    if acid_start is None:
        return False, "No acid component found"

    # Helper function to count carbons in a component
    def count_component_carbons(start_idx, exclude):
        visited = set(exclude)
        stack = [start_idx]
        carbon_count = 0
        while stack:
            atom_idx = stack.pop()
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6:
                carbon_count += 1
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in visited:
                    stack.append(neighbor.GetIdx())
        return carbon_count

    # Count carbons in both components
    exclude = {oxygen_idx, carbonyl_c_idx}
    alcohol_carbons = count_component_carbons(alcohol_start, exclude)
    acid_carbons = count_component_carbons(acid_start, exclude)

    # Check chain length requirements
    if alcohol_carbons < 8 or acid_carbons < 8:
        return False, f"Chain lengths too short (alcohol: {alcohol_carbons}, acid: {acid_carbons})"

    # Check for carboxylic acid groups (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OX2H1]")
    if mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Contains carboxylic acid group"

    return True, "Long-chain ester meeting wax criteria"