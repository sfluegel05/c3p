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

    # Count carbons - should be around 40 (allowing for some loss of methyls)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 35 or c_count > 45:
        return False, f"Carbon count ({c_count}) outside typical range for tetraterpenoids (35-45)"

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
    if double_bond_count < 6:
        return False, f"Too few double bonds ({double_bond_count}) for tetraterpenoid"

    # Check molecular weight - should be in typical range for tetraterpenoids
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500 or mol_wt > 800:
        return False, f"Molecular weight ({mol_wt}) outside typical range for tetraterpenoids"

    # Count rings - most tetraterpenoids have at least one
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count == 0:
        return False, "No rings found - unusual for tetraterpenoid"
    
    # Check for branching - tetraterpenoids are typically branched
    branching_points = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[*]([*])([*])[*]")))
    if branching_points < 2:
        return False, "Insufficient branching for tetraterpenoid"

    # Calculate degree of unsaturation
    du = rdMolDescriptors.CalcNumRotatableBonds(mol) + ring_count + double_bond_count
    if du < 15:
        return False, f"Insufficient degree of unsaturation ({du}) for tetraterpenoid"

    # Most tetraterpenoids contain oxygen (though not all)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    oxygen_info = f" Contains {o_count} oxygen atoms." if o_count > 0 else ""

    return True, (f"Matches tetraterpenoid pattern with {c_count} carbons, "
                 f"{methyl_count} methyl groups, {double_bond_count} double bonds, "
                 f"and {ring_count} rings.{oxygen_info}")