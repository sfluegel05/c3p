"""
Classifies: CHEBI:20156 3-oxo-Delta(1) steroid
"""
#!/usr/bin/env python
"""
Classifies: 3-oxo-Delta(1) steroid, defined as any 3-oxo steroid that contains a
double bond between positions 1 and 2.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_Delta_1__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(1) steroid.
    This function uses a simplified algorithm that:
      1. Checks if the molecule has four fused rings (steroid nucleus).
      2. Searches for a ring carbonyl ([R]C(=O)[R]) as the 3-oxo feature.
      3. Checks if the carbonyl carbon has an adjacent double bond (C=C) in the same ring,
         as an approximation of the Delta(1) (1,2-double bond) feature.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a 3-oxo-Delta(1) steroid, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information. Steroid nucleus (cyclopenta[a]phenanthrene) consists of 4 fused rings.
    ring_info = mol.GetRingInfo().AtomRings()
    if len(ring_info) != 4:
        return False, f"Expected 4 fused rings for a steroid nucleus, found {len(ring_info)} rings"
    
    # Find a ring carbonyl group by matching the SMARTS pattern that enforces the carbonyl to be in a ring.
    # [R] means the atom is in a ring.
    carbonyl_pattern = Chem.MolFromSmarts("[R]C(=O)[R]")
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    if not carbonyl_matches:
        return False, "No ring carbonyl (3-oxo) group found"
    
    # Convert each ring (tuple of atom indices) into a set for membership testing.
    rings = [set(r) for r in ring_info]
    
    # Check each matched ring carbonyl group and see if the carbonyl carbon has an adjacent C=C
    # bond where both atoms are contained in at least one common ring.
    has_delta1 = False
    for match in carbonyl_matches:
        # match tuple: (carbonyl carbon index, oxygen index, adjacent carbon index)
        carbonyl_idx = match[0]
        # Get all rings that contain the carbonyl atom:
        rings_with_carbonyl = [r for r in rings if carbonyl_idx in r]
        if not rings_with_carbonyl:
            continue  # should not happen because of the [R] in the SMARTS
        
        # Inspect all bonds from the carbonyl carbon
        atom = mol.GetAtomWithIdx(carbonyl_idx)
        for bond in atom.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                neighbor_idx = bond.GetOtherAtomIdx(carbonyl_idx)
                # Check if the neighbor also lies in one of the same rings as the carbonyl.
                for ring in rings_with_carbonyl:
                    if neighbor_idx in ring:
                        has_delta1 = True
                        break
            if has_delta1:
                break
        if has_delta1:
            break
    
    if not has_delta1:
        return False, "No double bond adjacent to the ring carbonyl found (Delta(1) feature missing)"
    
    return True, "Molecule has a steroid nucleus with a 3-oxo group and a Delta(1) double bond"

# Example usage (uncomment the following lines to test):
# smiles_example = "C[C@]12CC[C@H]3[C@@H](CC[C@@H]4CC(=O)C=C[C@]34C)[C@@H]1CC[C@@H]2O"  # example steroid SMILES
# print(is_3_oxo_Delta_1__steroid(smiles_example))