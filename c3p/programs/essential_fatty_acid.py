"""
Classifies: CHEBI:59549 essential fatty acid
"""
"""
Classifies: essential fatty acid (CHEBI:67011)
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is an essential fatty acid based on its SMILES string.
    Essential fatty acids are polyunsaturated with a carboxylic acid group, all cis double bonds,
    and no other functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an essential fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES"

    # Check for carboxylic acid group (including deprotonated form)
    carboxyl = Chem.MolFromSmarts("[CX3](=O)[OX2H1,O-]")
    if not mol.HasSubstructMatch(carboxyl):
        return False, "No carboxylic acid group"

    # Check oxygen count is exactly 2 (from COOH)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count != 2:
        return False, f"Expected 2 oxygen atoms, found {o_count}"

    # Check all other atoms are carbon or hydrogen
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [6, 1, 8]:
            return False, "Contains non-carbon/hydrogen/oxygen atoms"

    # Identify the carboxylic acid double bond to exclude
    carboxyl_matches = mol.GetSubstructMatches(carboxyl)
    if not carboxyl_matches:
        return False, "No carboxylic acid group found"
    # Get the carbonyl oxygen and adjacent atoms
    carboxyl_atoms = set()
    for match in carboxyl_matches:
        carboxyl_atoms.update(match)

    # Check all non-carboxyl double bonds are cis (Z)
    non_carboxyl_double_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            # Skip if part of carboxylic acid group
            atoms_in_bond = {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()}
            if atoms_in_bond.intersection(carboxyl_atoms):
                continue
            non_carboxyl_double_bonds += 1
            # Check stereo is cis (Z) or unspecified
            stereo = bond.GetStereo()
            if stereo not in [Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREONONE]:
                return False, f"Double bond {bond.GetIdx()} is trans or unspecified"
            # For STEREONONE, we can't confirm cis, so reject
            if stereo == Chem.BondStereo.STEREONONE:
                return False, "Unspecified double bond configuration"

    if non_carboxyl_double_bonds < 2:
        return False, f"Found {non_carboxyl_double_bonds} non-carboxyl double bonds, need at least 2"

    return True, "Polyunsaturated fatty acid with all cis double bonds and no other functional groups"