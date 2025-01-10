"""
Classifies: CHEBI:25409 monoterpenoid
"""
"""
Classifies: CHEBI:35189 monoterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.
    Monoterpenoids are derived from monoterpenes (C10 skeleton) but may have some
    carbons removed or rearranged, and often contain oxygen-containing functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get basic molecular properties
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Most monoterpenoids have 7-10 carbons (some may have lost methyl groups)
    if c_count < 7 or c_count > 12:
        return False, f"Carbon count ({c_count}) outside typical range for monoterpenoids (7-12)"
    
    # Check molecular weight - should typically be between 120-250 Da
    # (allowing for addition of oxygen-containing groups)
    if mol_wt < 100 or mol_wt > 300:
        return False, f"Molecular weight ({mol_wt:.1f}) outside typical range for monoterpenoids"

    # Look for common structural features
    
    # Isopropyl group pattern
    isopropyl_pattern = Chem.MolFromSmarts("CC(C)[CH1,CH2]")
    
    # Gem-dimethyl pattern
    gem_dimethyl_pattern = Chem.MolFromSmarts("C[C](C)")
    
    # Common cyclic patterns in monoterpenoids
    cyclohexane_pattern = Chem.MolFromSmarts("C1CCCCC1")
    cyclopentane_pattern = Chem.MolFromSmarts("C1CCCC1")
    
    # Check for presence of at least one of these characteristic patterns
    has_isopropyl = mol.HasSubstructMatch(isopropyl_pattern)
    has_gem_dimethyl = mol.HasSubstructMatch(gem_dimethyl_pattern)
    has_cyclohexane = mol.HasSubstructMatch(cyclohexane_pattern)
    has_cyclopentane = mol.HasSubstructMatch(cyclopentane_pattern)
    
    structural_features = []
    if has_isopropyl:
        structural_features.append("isopropyl group")
    if has_gem_dimethyl:
        structural_features.append("gem-dimethyl group")
    if has_cyclohexane:
        structural_features.append("cyclohexane ring")
    if has_cyclopentane:
        structural_features.append("cyclopentane ring")
        
    if not any([has_isopropyl, has_gem_dimethyl, has_cyclohexane, has_cyclopentane]):
        return False, "Missing characteristic monoterpenoid structural features"
    
    # Check for presence of oxygen-containing functional groups
    if o_count > 0:
        # Look for common oxygen-containing functional groups
        alcohol_pattern = Chem.MolFromSmarts("[OH1]")
        ether_pattern = Chem.MolFromSmarts("[OR0]")
        ketone_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[#6]")
        ester_pattern = Chem.MolFromSmarts("[#6][CX3](=[OX1])[OX2][#6]")
        
        has_alcohol = mol.HasSubstructMatch(alcohol_pattern)
        has_ether = mol.HasSubstructMatch(ether_pattern)
        has_ketone = mol.HasSubstructMatch(ketone_pattern)
        has_ester = mol.HasSubstructMatch(ester_pattern)
        
        if any([has_alcohol, has_ether, has_ketone, has_ester]):
            structural_features.append("oxygen-containing functional groups")

    # Calculate ring count
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count > 0:
        structural_features.append(f"{ring_count} ring(s)")

    # Final classification
    if len(structural_features) >= 2:  # Require at least two characteristic features
        return True, f"Monoterpenoid features found: {', '.join(structural_features)}"
    else:
        return False, "Insufficient characteristic features for monoterpenoid classification"