"""
Classifies: CHEBI:26935 tetraterpenoid
"""
"""
Classifies: CHEBI:26764 tetraterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    Tetraterpenoids are derived from tetraterpenes (C40 skeleton) and may have 
    modifications like rearrangements or removal of methyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetraterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbons - expanded range to include degraded and modified tetraterpenoids
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20 or c_count > 50:
        return False, f"Carbon count ({c_count}) outside typical range for tetraterpenoids (20-50)"

    # Count methyl groups - tetraterpenoids typically have multiple
    methyl_pattern = Chem.MolFromSmarts("[CH3]")
    methyl_count = len(mol.GetSubstructMatches(methyl_pattern))
    if methyl_count < 4:
        return False, f"Too few methyl groups ({methyl_count}) for tetraterpenoid"

    # Check for conjugated double bond system
    conjugated_pattern = Chem.MolFromSmarts("C=CC=CC=C")
    if not mol.HasSubstructMatch(conjugated_pattern):
        return False, "Missing conjugated double bond system"

    # Count double bonds - tetraterpenoids typically have many
    double_bond_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("C=C")))
    if double_bond_count < 5:
        return False, f"Too few double bonds ({double_bond_count}) for tetraterpenoid"

    # Check molecular weight - adjusted range
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 900:
        return False, f"Molecular weight ({mol_wt}) outside typical range for tetraterpenoids"

    # Count rings - but don't require them
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    
    # Check for excessive ring count (typical for non-tetraterpenoids)
    if ring_count > 4:
        return False, f"Too many rings ({ring_count}) for typical tetraterpenoid"
    
    # Count nitrogens - tetraterpenoids rarely contain nitrogen
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count > 1:
        return False, f"Too many nitrogens ({n_count}) for typical tetraterpenoid"
    
    # Count oxygens - check for excessive oxygen content
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count > 8:
        return False, f"Too many oxygens ({o_count}) for typical tetraterpenoid"

    # Check for branching - tetraterpenoids are typically branched
    branching_points = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[*]([*])([*])[*]")))
    if branching_points < 1:
        return False, "Insufficient branching for tetraterpenoid"

    # Calculate degree of unsaturation
    du = rdMolDescriptors.CalcNumRotatableBonds(mol) + ring_count + double_bond_count
    if du < 10:
        return False, f"Insufficient degree of unsaturation ({du}) for tetraterpenoid"

    ring_info = f" and {ring_count} rings" if ring_count > 0 else ""
    oxygen_info = f" Contains {o_count} oxygen atoms." if o_count > 0 else ""

    return True, (f"Matches tetraterpenoid pattern with {c_count} carbons, "
                 f"{methyl_count} methyl groups, {double_bond_count} double bonds{ring_info}.{oxygen_info}")