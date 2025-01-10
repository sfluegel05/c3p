"""
Classifies: CHEBI:26199 polyprenol
"""
"""
Classifies: polyprenol compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    Polyprenols are compounds with multiple isoprene units and a terminal hydroxyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for primary alcohol group (CH2-OH)
    primary_alcohol = Chem.MolFromSmarts("[CH2][OH1]")
    if not mol.HasSubstructMatch(primary_alcohol):
        return False, "No terminal primary alcohol group found"

    # Look for isoprene units with different possible configurations
    # Pattern matches both cis and trans configurations
    isoprene_pattern1 = Chem.MolFromSmarts("[CH2,CH3]-[C]=[C]-[C]")
    isoprene_pattern2 = Chem.MolFromSmarts("[CH2,CH3]-[C](-[CH3])=[C]-[CH2]")
    
    isoprene_matches1 = len(mol.GetSubstructMatches(isoprene_pattern1))
    isoprene_matches2 = len(mol.GetSubstructMatches(isoprene_pattern2))
    total_isoprene_matches = max(isoprene_matches1, isoprene_matches2)
    
    if total_isoprene_matches < 2:  # Need at least 2 isoprene units
        return False, "Less than 2 isoprene units found"

    # Count carbons and check if consistent with polyprenol structure
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:  # Minimum 2 isoprene units (C10)
        return False, f"Too few carbons ({c_count}) for a polyprenol"
    
    # Count methyls (branching points characteristic of isoprene units)
    methyl_pattern = Chem.MolFromSmarts("[CH3]-[C]")
    methyl_matches = len(mol.GetSubstructMatches(methyl_pattern))
    expected_methyls = total_isoprene_matches  # One methyl per isoprene unit
    if methyl_matches < expected_methyls - 1:  # Allow some flexibility
        return False, "Not enough methyl branches for polyprenol structure"

    # Count double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = len(mol.GetSubstructMatches(double_bond_pattern))
    if double_bond_matches < total_isoprene_matches - 1:  # One double bond per isoprene unit (except possibly terminal)
        return False, "Not enough double bonds for polyprenol structure"

    # Check for rings - polyprenols should be mostly linear
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings > 1:  # Allow 1 ring as some examples might have a single ring
        return False, f"Too many rings ({num_rings}) for a polyprenol"

    # Check for continuous carbon chain using rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < (total_isoprene_matches * 2 - 2):  # Expect ~2 rotatable bonds per isoprene unit
        return False, "Carbon chain not long enough for claimed number of isoprene units"

    # Check oxygen count - should only have 1 oxygen (the terminal OH)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1:
        return False, "No oxygen found"
    if o_count > 5:  # Allow a few extra oxygens for modified polyprenols
        return False, "Too many oxygens for a typical polyprenol"

    # Success case - molecule meets all criteria
    return True, f"Contains approximately {total_isoprene_matches} isoprene units with terminal OH group"