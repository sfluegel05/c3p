"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
"""
Classifies: CHEBI:35189 sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    Sesquiterpenoids are terpenoids derived from sesquiterpenes (C15 skeleton),
    possibly modified by rearrangement or removal of some carbons.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Sesquiterpenoids should have around 15 carbons (allowing variation due to modifications)
    if c_count < 11 or c_count > 25:
        return False, f"Carbon count ({c_count}) outside typical range for sesquiterpenoids (11-25)"

    # Calculate rings
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    
    # Most sesquiterpenoids are cyclic, but some can be linear
    if ring_count > 5:
        return False, f"Too many rings ({ring_count}) for a typical sesquiterpenoid"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 180 or mol_wt > 600:
        return False, f"Molecular weight ({mol_wt:.1f}) outside typical range for sesquiterpenoids"

    # Count double bonds using SMARTS pattern
    double_bond_pattern = Chem.MolFromSmarts('C=C')
    if double_bond_pattern:
        double_bond_count = len(mol.GetSubstructMatches(double_bond_pattern))
    else:
        double_bond_count = 0

    # Look for typical functional groups found in sesquiterpenoids
    functional_groups = {
        'alcohol': '[OH1]',
        'ketone': '[CX3](=[OX1])[CX4]',
        'carboxylic_acid': '[CX3](=O)[OX2H1]',
        'ester': '[#6][CX3](=O)[OX2H0][#6]',
        'epoxide': '[C]1O[C]1',
        'ether': '[OX2]([CX4])[CX4]'
    }
    
    found_groups = []
    for group_name, smarts in functional_groups.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            found_groups.append(group_name)

    # Look for branching methyl groups - common in terpenoids
    methyl_pattern = Chem.MolFromSmarts('[CH3]')
    if methyl_pattern:
        methyl_count = len(mol.GetSubstructMatches(methyl_pattern))
    else:
        methyl_count = 0

    # Common sesquiterpenoid ring patterns
    ring_patterns = [
        ('[C]1[C][C][C][C][C]1', 'cyclohexane'),  # 6-membered ring
        ('[C]1[C][C][C][C]1', 'cyclopentane'),    # 5-membered ring
        ('[C]1[C][C][C]1', 'cyclobutane'),        # 4-membered ring
    ]
    
    found_rings = []
    for smarts, ring_name in ring_patterns:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            found_rings.append(ring_name)

    # Scoring system
    score = 0
    reasons = []
    
    # Base score from carbon count
    if 13 <= c_count <= 17:
        score += 2
        reasons.append(f"Typical sesquiterpenoid carbon count ({c_count})")
    elif 11 <= c_count <= 25:
        score += 1
        reasons.append(f"Acceptable carbon count ({c_count})")

    # Score from rings
    if 1 <= ring_count <= 4:
        score += 2
        reasons.append(f"Typical ring count ({ring_count})")
    elif ring_count == 0 and double_bond_count >= 2:
        score += 1
        reasons.append("Linear with multiple double bonds")

    # Score from functional groups
    if found_groups:
        score += len(found_groups)
        reasons.append(f"Contains {', '.join(found_groups)}")

    # Score from methyl groups
    if 3 <= methyl_count <= 6:
        score += 2
        reasons.append(f"Contains {methyl_count} methyl groups")
    elif methyl_count > 0:
        score += 1
        reasons.append(f"Contains {methyl_count} methyl groups")

    # Final decision
    if score >= 4:
        return True, "; ".join(reasons)
    else:
        return False, "Insufficient sesquiterpenoid characteristics"