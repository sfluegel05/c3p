"""
Classifies: CHEBI:50128 biflavonoid
"""
"""
Classifies: Biflavonoid
Definition: A flavonoid oligomer that is obtained by the oxidative coupling of at least two units of aryl‐substituted benzopyran rings or its substituted derivatives, 
resulting in the two ring systems being joined together by a single atom or bond.

This improved heuristic:
  - Requires a moderate number of rings (≥6) and a molecular weight > 400 Da.
  - Uses two more permissive SMARTS patterns for flavonoid substructures:
       * A "flavone‐like" unit: roughly "c1ccc2c(c1)OC(=O)c(c2)"
       * A "flavan‐like" unit: roughly "c1ccc2c(c1)OC(c2)"
  - Groups overlapping matches (if they share >1 atom) into candidate flavonoid units.
  - Checks that at least two units are found and that a linking bond (or a shared atom) is present.
Note: This is only a heuristic and may mis‐classify some edge cases.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_biflavonoid(smiles: str):
    """
    Determines if a molecule is a biflavonoid based on its SMILES string.
    The improved heuristic uses two more permissive flavonoid SMARTS patterns.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is considered a biflavonoid, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Sanity checks: check for a reasonable ring count and molecular weight.
    ring_info = mol.GetRingInfo()
    num_rings = len(ring_info.AtomRings())
    if num_rings < 6:
        return False, f"Too few rings ({num_rings}) to be a biflavonoid"
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low for a biflavonoid"

    # Define two more permissive SMARTS for flavonoid-like units.
    # Pattern 1: Flavone-like. (contains an aromatic fused system with a carbonyl)
    flavone_smarts = "c1ccc2c(c1)OC(=O)c(c2)"
    flavone_pattern = Chem.MolFromSmarts(flavone_smarts)
    if flavone_pattern is None:
        return False, "Error generating flavone pattern"

    # Pattern 2: Flavan-like. (a similar fused aromatic system without an explicit carbonyl)
    flavan_smarts = "c1ccc2c(c1)OC(c2)"
    flavan_pattern = Chem.MolFromSmarts(flavan_smarts)
    if flavan_pattern is None:
        return False, "Error generating flavan pattern"

    # Collect all substructure matches (each as a frozenset of atom indices) from both patterns.
    matches = []
    for pattern in (flavone_pattern, flavan_pattern):
        for match in mol.GetSubstructMatches(pattern):
            # Save each match as a frozenset for grouping purposes.
            matches.append(frozenset(match))

    if not matches:
        return False, "No flavonoid-like substructure found"

    # Group matches that likely represent the same flavonoid unit.
    groups = []
    for match in matches:
        added = False
        for group in groups:
            # If this match overlaps with any match in the group by more than one atom,
            # then we consider it to belong to the same candidate unit.
            if any(len(match & other) > 1 for other in group):
                group.append(match)
                added = True
                break
        if not added:
            groups.append([match])
    num_units = len(groups)
    if num_units < 2:
        return False, f"Found only {num_units} flavonoid unit(s); at least two are needed for a biflavonoid"

    # Check that at least one linking connection exists between the two candidate units.
    # We require that the union of atoms of one unit and the other are connected by a bond
    # or share a very small overlap (up to 2 atoms).
    connection_found = False
    for i in range(len(groups)):
        group_i_atoms = set.union(*groups[i])
        for j in range(i+1, len(groups)):
            group_j_atoms = set.union(*groups[j])
            # Check direct overlap first.
            common = group_i_atoms & group_j_atoms
            if common and len(common) <= 2:
                connection_found = True
                break
            # If not overlapping, check if any atom in group_i is bonded to any in group_j.
            for a in group_i_atoms:
                for b in group_j_atoms:
                    if mol.GetBondBetweenAtoms(a, b) is not None:
                        connection_found = True
                        break
                if connection_found:
                    break
            if connection_found:
                break
        if connection_found:
            break
    if not connection_found:
        return False, "Did not detect a linking atom or bond between candidate flavonoid units"

    return True, f"Found {num_units} flavonoid unit(s) joined by a short linker; molecular weight = {mol_wt:.1f} Da"

# For sanity testing, you may call the function on one of the examples.
if __name__ == "__main__":
    test_smiles = "Oc1ccc(cc1)[C@@H]1Oc2cc(O)cc(O)c2C(=O)[C@H]1c1c(O)cc(O)c2c1oc(cc2=O)-c1ccc(O)c(O)c1"  # (+)-morelloflavone
    result, reason = is_biflavonoid(test_smiles)
    print("Biflavonoid:", result)
    print("Reason:", reason)