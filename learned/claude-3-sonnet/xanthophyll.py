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
    if c_count < 35 or c_count > 45:
        return False, f"Carbon count {c_count} outside typical range for xanthophylls (35-45)"
    
    # Must have oxygen atoms
    if o_count == 0:
        return False, "No oxygen atoms found - xanthophylls must be oxygenated"

    # Look for conjugated polyene patterns
    conjugated_patterns = [
        "C=CC=CC=CC=C",  # Basic conjugated system
        "C=CC=CC=CC=CC=C"  # Extended conjugation
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
        "hydroxy": "[OH]",
        "keto": "[CH0](=O)",
        "epoxy": "[O][CH1][CH1]",
        "ester": "C(=O)[O]"
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

    # Look for cyclic end groups (common in xanthophylls)
    ring_pattern = Chem.MolFromSmarts("C1CCCCC1")
    if ring_pattern and not mol.HasSubstructMatch(ring_pattern):
        return False, "Missing characteristic cyclic end groups"

    # Count rings
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count == 0:
        return False, "No rings found - most xanthophylls have cyclic end groups"

    # Calculate double bonds
    double_bonds = rdMolDescriptors.CalcNumDoubleBonds(mol)
    if double_bonds < 8:
        return False, "Insufficient double bonds for xanthophyll conjugated system"

    # Check for characteristic methyl branching
    methyl_pattern = Chem.MolFromSmarts("CC(C)C")
    if methyl_pattern and not mol.HasSubstructMatch(methyl_pattern):
        return False, "Missing characteristic methyl branching"

    # Additional check for molecular weight range (typical for xanthophylls)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500 or mol_wt > 800:
        return False, f"Molecular weight {mol_wt:.1f} outside typical range for xanthophylls"

    # If all checks pass, classify as xanthophyll
    oxygen_groups_str = ", ".join(f"{count} {name}" for name, count in oxygen_group_counts.items() if count > 0)
    return True, (f"Oxygenated carotenoid with {ring_count} rings, {double_bonds} double bonds, "
                 f"and oxygen groups: {oxygen_groups_str}")