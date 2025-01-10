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

    # Look for the characteristic repeating isoprene unit pattern
    # This matches the basic isoprene unit structure with proper connectivity
    isoprene_pattern = Chem.MolFromSmarts("[CH2,CH3][C]([CH3])=[CH][CH2]")
    isoprene_matches = len(mol.GetSubstructMatches(isoprene_pattern))
    
    if isoprene_matches < 2:
        return False, f"Insufficient isoprene units ({isoprene_matches}), need at least 2"

    # Check for continuous carbon chain characteristic of polyprenols
    # This pattern looks for connected isoprene units
    chain_pattern = Chem.MolFromSmarts("[CH2,CH3][C]([CH3])=[CH][CH2][C]([CH3])=[CH][CH2]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No continuous isoprene chain found"

    # Count carbons and check if consistent with polyprenol structure
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    expected_c_min = isoprene_matches * 5  # Minimum expected carbons
    if c_count < expected_c_min:
        return False, f"Too few carbons ({c_count}) for claimed number of isoprene units"

    # Verify linear structure - polyprenols should be mostly linear
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings > 0:  # Polyprenols should be linear
        return False, f"Contains rings ({num_rings}), should be linear"

    # Check for proper double bond pattern
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bonds = len(mol.GetSubstructMatches(double_bond_pattern))
    if double_bonds < isoprene_matches:
        return False, "Insufficient number of double bonds for polyprenol structure"

    # Check for proper methyl branching pattern
    methyl_pattern = Chem.MolFromSmarts("[CH3][C]=[CH][CH2]")
    methyl_branches = len(mol.GetSubstructMatches(methyl_pattern))
    if methyl_branches < isoprene_matches - 1:
        return False, "Incorrect methyl branching pattern"

    # Check for sugar moieties or other complex substituents
    sugar_pattern = Chem.MolFromSmarts("[OH1][CH1]([OH1])[CH1]([OH1])")
    if mol.HasSubstructMatch(sugar_pattern):
        return False, "Contains sugar moiety"

    # Success case - molecule meets all criteria
    return True, f"Contains {isoprene_matches} isoprene units in polyprenol arrangement"