"""
Classifies: CHEBI:16219 cucurbitacin
"""
from rdkit import Chem

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin based on its SMILES string.
    Cucurbitacins are tetracyclic triterpenoids derived from cucurbitane.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool or None: True if molecule is a cucurbitacin, False otherwise, or None if unable to determine
        str or None: Reason for classification or None
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the cucurbitane core skeleton SMARTS (tetracyclic triterpenoid structure)
    cucurbitane_smarts = '[#6]12CC3CC(C1)CC4CCC(C2)C34'  # Corrected pattern without line breaks
    cucurbitane_pattern = Chem.MolFromSmarts(cucurbitane_smarts)
    if cucurbitane_pattern is None:
        return None, "Unable to create cucurbitane SMARTS pattern"

    # Check for cucurbitane skeleton
    if not mol.HasSubstructMatch(cucurbitane_pattern):
        return False, "Cucurbitane skeleton not found"

    # Check for ketone groups (C=O)
    ketone_pattern = Chem.MolFromSmarts("C=O")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    if len(ketone_matches) == 0:
        return False, "No ketone groups found"

    # Check for hydroxyl groups (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) == 0:
        return False, "No hydroxyl groups found"

    # Check for double bonds in rings (unsaturation)
    double_bond_in_ring = False
    for bond in mol.GetBonds():
        if bond.IsInRing() and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            double_bond_in_ring = True
            break
    if not double_bond_in_ring:
        return False, "No double bonds found in rings"

    # If all checks pass, classify as cucurbitacin
    return True, "Contains cucurbitane skeleton with characteristic functional groups"