"""
Classifies: CHEBI:27325 xanthophyll
"""
"""
Classifies: xanthophylls (oxygenated carotenoids)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    Xanthophylls are oxygenated carotenoids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a xanthophyll, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Allow C20-C50 to include apo-carotenoids
    if c_count < 20 or c_count > 50:
        return False, f"Carbon count {c_count} outside typical range for xanthophylls (20-50)"
    
    # Must have oxygen atoms
    if o_count == 0:
        return False, "No oxygen atoms found - xanthophylls must be oxygenated"

    # Look for characteristic carotenoid backbone patterns
    carotenoid_patterns = [
        "C=CC=CC=CC=CC=C",  # Extended conjugation
        "CC(C)=CC=CC=C",    # Typical end portion
        "C=CC(C)=CC=CC=C",  # Methylated conjugation
        "C=CC=C(C)C=CC=C"   # Alternative methylation
    ]
    
    backbone_matches = 0
    for pattern in carotenoid_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat and mol.HasSubstructMatch(pat):
            backbone_matches += len(mol.GetSubstructMatches(pat))
    
    if backbone_matches < 2:
        return False, "Missing characteristic carotenoid conjugated system"

    # Look for oxygen-containing functional groups in typical xanthophyll positions
    oxygen_patterns = {
        "hydroxy": "[OX2H1]C",  # Hydroxyl group
        "keto": "[CX3](=O)[CX4]",  # Ketone
        "epoxy": "[OX2r3]1[CX4r3][CX4r3]1",  # Epoxide
        "ether": "[OX2]([CX4])[CX4]"  # Ether
    }
    
    oxygen_group_counts = {}
    total_o_groups = 0
    
    for name, pattern in oxygen_patterns.items():
        pat = Chem.MolFromSmarts(pattern)
        if pat:
            matches = len(mol.GetSubstructMatches(pat))
            oxygen_group_counts[name] = matches
            total_o_groups += matches
    
    if total_o_groups == 0:
        return False, "No characteristic oxygen-containing functional groups found"

    # Count double bonds
    double_bond_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("C=C")))
    if double_bond_count < 5:
        return False, f"Insufficient conjugated system (found {double_bond_count} double bonds, need at least 5)"

    # Look for alternating methyl branching pattern characteristic of carotenoids
    methyl_patterns = [
        "CC(C)=CC=C",  # Terminal methyl pattern
        "C=C(C)C=CC",  # Internal methyl pattern
        "C=CC(C)=CC"   # Alternative internal pattern
    ]
    
    methyl_matches = 0
    for pattern in methyl_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat:
            methyl_matches += len(mol.GetSubstructMatches(pat))
    
    if methyl_matches < 2:
        return False, "Missing characteristic carotenoid methyl branching pattern"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350 or mol_wt > 850:
        return False, f"Molecular weight {mol_wt:.1f} outside typical range for xanthophylls"

    # If all checks pass, classify as xanthophyll
    features = []
    features.append(f"{double_bond_count} double bonds")
    features.append(f"{o_count} oxygen atoms")
    oxygen_groups_str = ", ".join(f"{count} {name}" for name, count in oxygen_group_counts.items() if count > 0)
    features.append(f"oxygen groups: {oxygen_groups_str}")
    
    return True, "Oxygenated carotenoid with " + "; ".join(features)