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
    Polyprenols have the general formula H-[CH2C(Me)=CHCH2]nOH where n>1.

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

    # Count total oxygens - polyprenols typically have just one OH group
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count > 10:  # Allow some additional O atoms for modified polyprenols
        return False, f"Too many oxygens ({o_count}) for a polyprenol"

    # Look for the characteristic repeating isoprene unit pattern
    # This pattern matches the basic isoprene unit structure: -CH2-C(CH3)=CH-CH2-
    isoprene_pattern = Chem.MolFromSmarts("[CH2][C]([CH3])=[CH][CH2]")
    isoprene_matches = len(mol.GetSubstructMatches(isoprene_pattern))
    
    # Alternative pattern to catch some variations
    alt_pattern = Chem.MolFromSmarts("[CH2,CH3][C]([CH3])=[CH][CH2]")
    alt_matches = len(mol.GetSubstructMatches(alt_pattern))
    
    total_isoprene_units = max(isoprene_matches, alt_matches)
    
    if total_isoprene_units < 1:
        return False, "No isoprene units found"

    # Count carbons - should be 5n+1 where n is number of isoprene units
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    expected_c_min = total_isoprene_units * 5  # Minimum expected carbons
    if c_count < expected_c_min:
        return False, f"Too few carbons ({c_count}) for claimed number of isoprene units"

    # Check for double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bonds = len(mol.GetSubstructMatches(double_bond_pattern))
    if double_bonds < total_isoprene_units:
        return False, "Insufficient number of double bonds for polyprenol structure"

    # Verify linear chain structure - polyprenols should be mostly linear
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings > 1:  # Allow at most one ring as some natural polyprenols might have modifications
        return False, f"Too many rings ({num_rings}) for a polyprenol"

    # Check for branching pattern characteristic of polyprenols
    methyl_branch_pattern = Chem.MolFromSmarts("[CH3][C]=[C]")
    methyl_branches = len(mol.GetSubstructMatches(methyl_branch_pattern))
    if methyl_branches < total_isoprene_units - 1:  # Allow some flexibility
        return False, "Insufficient methyl branching pattern for polyprenol"

    # Check for continuous carbon chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    min_rotatable = total_isoprene_units * 2 - 1  # Expect at least 2 rotatable bonds per unit
    if n_rotatable < min_rotatable:
        return False, "Carbon chain not long enough for polyprenol structure"

    # Success case - molecule meets all criteria
    return True, f"Contains {total_isoprene_units} isoprene units in polyprenol arrangement"