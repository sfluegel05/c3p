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

    # Look for isoprene units (C=C-C(C)-C)
    isoprene_pattern = Chem.MolFromSmarts("[CH2,CH3]-[C]=[C]-[C]")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    
    if len(isoprene_matches) < 2:  # Need at least 2 isoprene units
        return False, "Less than 2 isoprene units found"

    # Count carbons and check if multiple of 5 (±1 for variations)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:  # Minimum 2 isoprene units (C10)
        return False, f"Too few carbons ({c_count}) for a polyprenol"
    
    # Count methyls (branching points characteristic of isoprene units)
    methyl_pattern = Chem.MolFromSmarts("[CH3]-[C]")
    methyl_matches = len(mol.GetSubstructMatches(methyl_pattern))
    if methyl_matches < 2:
        return False, "Not enough methyl branches for polyprenol structure"

    # Count double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = len(mol.GetSubstructMatches(double_bond_pattern))
    if double_bond_matches < 2:
        return False, "Not enough double bonds for polyprenol structure"

    # Calculate degree of unsaturation
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings > 1:  # Allow 1 ring as some examples might have a single ring
        return False, f"Too many rings ({num_rings}) for a polyprenol"

    # Additional check for carbon chain length
    longest_chain = max(len(path) for path in rdMolDescriptors.FindAllPathsOfLengthN(mol, 8))
    if longest_chain < 8:
        return False, "Carbon chain too short for polyprenol"

    # Success case - molecule meets all criteria
    return True, f"Contains {len(isoprene_matches)} isoprene units with terminal OH group"


__metadata__ = {
    'chemical_class': {
        'name': 'polyprenol',
        'definition': 'Any member of the class of prenols possessing the general formula '
                     'H-[CH2C(Me)=CHCH2]nOH in which the carbon skeleton is composed of '
                     'more than one isoprene units.'
    }
}