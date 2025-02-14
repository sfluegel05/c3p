"""
Classifies: CHEBI:16219 cucurbitacin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin based on its SMILES string.
    Cucurbitacins are tetracyclic triterpenoids derived from cucurbitane.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cucurbitacin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the cucurbitane core skeleton SMARTS (tetracyclic triterpenoid structure)
    # The pattern represents four fused cyclohexane rings with the appropriate connectivity
    cucurbitane_smarts = """
    [#6]1([#6])[#6][#6][#6]2[#6]1[#6][#6][#6]3[#6]2[#6][#6][#6]4[#6]3[#6][#6][#6][#6]4
    """
    cucurbitane_pattern = Chem.MolFromSmarts(cucurbitane_smarts)
    if cucurbitane_pattern is None:
        return False, "Unable to create cucurbitane SMARTS pattern"

    # Check for cucurbitane skeleton
    if not mol.HasSubstructMatch(cucurbitane_pattern):
        return False, "Cucurbitane skeleton not found"

    # Check for characteristic functional groups (e.g., ketones, hydroxyls)
    ketone_pattern = Chem.MolFromSmarts("C(=O)[#6]")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    if len(ketone_matches) == 0:
        return False, "No ketone groups found"

    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) == 0:
        return False, "No hydroxyl groups found"

    # Check for acylated hydroxyl groups (e.g., acetyl esters)
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    # Count the number of rings
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 4:
        return False, f"Insufficient number of rings ({num_rings})"

    # Check for double bonds in specific positions if necessary
    # For cucurbitacins, double bonds may be present in the ring system or side chains

    # Optionally, check molecular weight range typical for cucurbitacins
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500 or mol_wt > 1200:
        return False, "Molecular weight not in typical range for cucurbitacins"

    return True, "Contains cucurbitane skeleton with characteristic functional groups"