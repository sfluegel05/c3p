"""
Classifies: CHEBI:50128 biflavonoid
"""
"""
Classifies: Biflavonoid
Definition: A flavonoid oligomer that is obtained by the oxidative coupling of at least two units of aryl‐substituted benzopyran rings or its substituted derivatives, 
resulting in the two ring systems being joined together by a single atom or bond.
This program heuristically searches for two flavonoid‐like substructures in the input SMILES.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_biflavonoid(smiles: str):
    """
    Determines if a molecule is a biflavonoid based on its SMILES string.
    The heuristic uses substructure matching: it looks for at least two distinct flavonoid-like units.
    
    We define two simple SMARTS for flavonoid substructures:
      - A "flavan" type: a benzopyran where the heterocycle is saturated.
      - A "flavone" type: a benzopyran with a carbonyl group (chromone type).
    
    We then group overlapping matches (allowing up to one common atom as a joining point)
    to represent one flavonoid unit and require at least two such units.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a biflavonoid, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # As a sanity check, biflavonoids generally have a moderate to high number of rings.
    ring_info = mol.GetRingInfo()
    num_rings = len(ring_info.AtomRings())
    if num_rings < 6:
        return False, f"Too few rings ({num_rings}) to be a biflavonoid"
    
    # Define two heuristic substructure patterns for a flavonoid unit.
    # Pattern 1 (flavan-like): fused benzene and saturated oxygen-containing heterocycle.
    # e.g. matches substructure like: c1ccc2OCc2c1 
    flavan_pattern = Chem.MolFromSmarts("c1ccc2OCc2c1")
    if flavan_pattern is None:
        return False, "Error in generating flavan pattern"
    
    # Pattern 2 (flavone-like): fused benzene and unsaturated heterocycle with a carbonyl.
    # e.g. matches substructure like: c1ccc2c(c1)OC(=O)c2
    flavone_pattern = Chem.MolFromSmarts("c1ccc2c(c1)OC(=O)c2")
    if flavone_pattern is None:
        return False, "Error in generating flavone pattern"

    # Find all matches for either pattern
    matches = []
    for pattern in (flavan_pattern, flavone_pattern):
        sub_matches = mol.GetSubstructMatches(pattern)
        for m in sub_matches:
            # Store each match as a frozenset of atom indices for easier grouping
            matches.append(frozenset(m))
    
    if not matches:
        return False, "No flavonoid-like substructure found"
    
    # Group matches that represent the same flavonoid unit.
    # Two matches will be grouped together if they share more than 1 atom (allowing a connection by one atom).
    # We use a union-find like approach.
    groups = []
    
    for match in matches:
        added = False
        for group in groups:
            # if this match overlaps with any match already in the group by more than one atom, consider it part of the same unit.
            if any(len(match & other) > 1 for other in group):
                group.append(match)
                added = True
                break
        if not added:
            groups.append([match])
    
    num_units = len(groups)
    if num_units < 2:
        return False, f"Found only {num_units} flavonoid unit(s); at least two are needed for a biflavonoid"

    # Optional additional checks:
    # Check that the two flavonoid units are joined by a short linker.
    # We find the intersection (shared atoms) between any two groups.
    connection_found = False
    for i in range(len(groups)):
        # Get the union of atoms for the entire group i.
        group_i_atoms = set.union(*groups[i])
        for j in range(i+1, len(groups)):
            group_j_atoms = set.union(*groups[j])
            # Intersection might represent the linking atom(s)
            if group_i_atoms & group_j_atoms:
                # Since a biflavonoid is connected by a single atom or bond, we allow a small overlap.
                if len(group_i_atoms & group_j_atoms) <= 2:
                    connection_found = True
                    break
        if connection_found:
            break
    
    if not connection_found:
        # It is possible the two units are entirely separate which would be unusual for a biflavonoid.
        return False, "Did not detect a linking atom or bond between flavonoid units"
    
    # As an extra check, we verify that the molecular weight is in an expected range (typically >400 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low for a biflavonoid"
    
    return True, f"Found {num_units} flavonoid unit(s) joined by a short linker; molecular weight = {mol_wt:.1f} Da"

# For testing purposes, you could call the function with one of the examples:
if __name__ == "__main__":
    test_smiles = "Oc1ccc(cc1)[C@@H]1Oc2cc(O)cc(O)c2C(=O)[C@H]1c1c(O)cc(O)c2c1oc(cc2=O)-c1ccc(O)c(O)c1"  # (+)-morelloflavone
    result, reason = is_biflavonoid(test_smiles)
    print("Biflavonoid:", result)
    print("Reason:", reason)