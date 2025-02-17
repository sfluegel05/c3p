"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
"""
Classifies: triterpenoid saponin
Definition: A terpene glycoside in which the terpene moiety is a triterpenoid.
That is, the molecule must have at least one sugar (glycoside) group and a terpenoid
aglycone with roughly 30 carbon atoms and several fused rings.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin based on its SMILES string.
    
    Heuristic criteria:
      1. The molecule must be a valid structure.
      2. At least one sugar moiety must be present.
         (Here we search for a typical pyranose sugar pattern.)
      3. The aglycone part (i.e. total carbons minus those in sugar fragments) should
         be near 30 carbons (allowing a margin of ±3).
      4. The non-sugar (aglycone) part should have at least 4 rings (typical for
         triterpenoid cores).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a triterpenoid saponin, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Define a SMARTS pattern for a typical pyranose sugar unit.
    # This pattern is not universal but should pick up many sugar rings.
    sugar_smarts = "[C@H]1OC(CO)C(O)C(O)C1O"
    sugar_pattern = Chem.MolFromSmarts(sugar_smarts)
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No sugar (glycoside) moiety detected."
    
    # Gather unique atom indices that belong to any sugar fragment found.
    sugar_atom_indices = set()
    for match in sugar_matches:
        for idx in match:
            sugar_atom_indices.add(idx)
    
    # Count total carbons in the molecule.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    # Count carbons that are in sugar fragment(s).
    sugar_carbons = 0
    for idx in sugar_atom_indices:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 6:
            sugar_carbons += 1
    # Estimate aglycone (terpenoid) carbon count.
    aglycone_carbons = total_carbons - sugar_carbons

    # For a triterpenoid the aglycone should be ~30 carbons (allow a small margin).
    if not (27 <= aglycone_carbons <= 33):
        return False, f"Aglycone carbons not around 30 (found {aglycone_carbons})."
    
    # Check ring systems. A triterpenoid core usually exhibits multiple fused rings.
    # Get all rings in the molecule.
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    # Count rings that are not exclusively part of a sugar ring.
    non_sugar_ring_count = 0
    for ring in all_rings:
        # A sugar ring typically is 6-membered and all atoms are in one of sugar_matches.
        # Check if at least one atom in the ring is not in any sugar fragment.
        if not all(idx in sugar_atom_indices for idx in ring):
            non_sugar_ring_count += 1

    if non_sugar_ring_count < 4:
        return False, f"Insufficient number of rings in aglycone (found {non_sugar_ring_count})."
    
    # Optionally, one can check molecular weight or other descriptors:
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    # Many triterpenoid saponins have significant molecular weight.
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a triterpenoid saponin."
    
    return True, "Molecule contains a sugar moiety plus a triterpenoid-like aglycone (approx. 30 carbons and multiple rings)."