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
    
    # Xanthophylls are typically C40 terpenoids
    if c_count < 30 or c_count > 50:
        return False, f"Carbon count {c_count} outside typical range for xanthophylls (30-50)"
    
    # Must have oxygen atoms
    if o_count == 0:
        return False, "No oxygen atoms found - xanthophylls must be oxygenated"

    # Look for conjugated polyene patterns (characteristic of carotenoids)
    conjugated_patterns = [
        "C=CC=CC=CC=C",  # Basic conjugated system
        "C=CC=CC=CC=CC=C",  # Extended conjugation
        r"[CH2,CH3]/C=C/C=C/C=C"  # Trans conjugated system
    ]
    
    has_conjugation = False
    for pattern in conjugated_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat and mol.HasSubstructMatch(pat):
            has_conjugation = True
            break
            
    if not has_conjugation:
        return False, "Missing characteristic conjugated polyene system"

    # Look for oxygen-containing functional groups
    oxygen_patterns = {
        "hydroxy": "[OX2H1]",  # Hydroxyl group
        "keto": "[CX3](=O)[CX4]",  # Ketone
        "epoxy": "[OX2r3]1[CX4r3][CX4r3]1",  # Epoxide
        "ester": "[CX3](=O)[OX2]",  # Ester
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
        return False, "No oxygen-containing functional groups found"

    # Look for characteristic end groups (cyclohexene rings common in xanthophylls)
    end_group_patterns = [
        "C1=C(C)CCCC1(C)C",  # Common end group
        "C1C(C)=CCCC1(C)C",  # Alternative end group
    ]
    
    has_end_group = False
    for pattern in end_group_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat and mol.HasSubstructMatch(pat):
            has_end_group = True
            break

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 450 or mol_wt > 850:
        return False, f"Molecular weight {mol_wt:.1f} outside typical range for xanthophylls"

    # Count double bonds by looking at bond types
    double_bond_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("C=C")))
    if double_bond_count < 7:
        return False, f"Insufficient conjugated system (found {double_bond_count} double bonds, need at least 7)"

    # Check for characteristic branching methyl groups
    methyl_pattern = Chem.MolFromSmarts("C[C](C)C")
    methyl_groups = len(mol.GetSubstructMatches(methyl_pattern))
    if methyl_groups < 2:
        return False, "Missing characteristic methyl branching"

    # If all checks pass, classify as xanthophyll
    features = []
    if has_end_group:
        features.append("cyclic end groups")
    features.append(f"{double_bond_count} double bonds")
    features.append(f"{o_count} oxygen atoms")
    oxygen_groups_str = ", ".join(f"{count} {name}" for name, count in oxygen_group_counts.items() if count > 0)
    features.append(f"oxygen groups: {oxygen_groups_str}")
    
    return True, "Oxygenated carotenoid with " + "; ".join(features)