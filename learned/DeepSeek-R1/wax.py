"""
Classifies: CHEBI:73702 wax
"""
"""
Classifies: Wax (long-chain ester)
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

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

    # Helper function to analyze component
    def analyze_component(start_idx, exclude):
        visited = set(exclude)
        stack = [start_idx]
        carbon_count = 0
        in_ring = False
        while stack:
            atom_idx = stack.pop()
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            
            # Check for rings
            if atom.IsInRing():
                in_ring = True
                
            # Count carbons
            if atom.GetAtomicNum() == 6:
                carbon_count += 1
                
            # Check for forbidden groups (amide, nitrile, etc.)
            for bond in atom.GetBonds():
                neighbor_idx = bond.GetOtherAtomIdx(atom_idx)
                if neighbor_idx in visited:
                    continue
                stack.append(neighbor_idx)
                
        return carbon_count, in_ring

    # Analyze both components
    exclude = {oxygen_idx, carbonyl_c_idx}
    alcohol_carbons, alcohol_in_ring = analyze_component(alcohol_start, exclude)
    acid_carbons, acid_in_ring = analyze_component(acid_start, exclude)

    # Check for rings in components
    if alcohol_in_ring or acid_in_ring:
        return False, "Component in ring structure"

    # Check chain lengths (minimum 8 carbons each)
    if alcohol_carbons < 8 or acid_carbons < 8:
        return False, f"Chain lengths too short (alcohol: {alcohol_carbons}, acid: {acid_carbons})"

    # Check for other functional groups
    forbidden_patterns = {
        "carboxylic_acid": Chem.MolFromSmarts("C(=O)[OH]"),
        "amide": Chem.MolFromSmarts("[NX3][CX3](=[OX1])"),
        "nitrile": Chem.MolFromSmarts("C#N"),
        "sulfonic_acid": Chem.MolFromSmarts("S(=O)(=O)O"),
        "phosphate": Chem.MolFromSmarts("OP(=O)(O)O"),
    }
    for name, pattern in forbidden_patterns.items():
        if mol.HasSubstructMatch(pattern):
            return False, f"Contains {name.replace('_', ' ')} group"

    # Check molecular weight (typical waxes > 300 Da)
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da)"

    return True, "Long-chain ester meeting wax criteria"