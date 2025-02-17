"""
Classifies: CHEBI:46633 carbapenems
"""
"""
Classifies: Carbapenems (beta-lactam antibiotics having a carbapenem skeleton,
which is a four-membered beta-lactam ring fused with a five-membered ring that does not contain sulfur).

The approach:
1. Parse the SMILES.
2. Scan the ring systems:
   a. Look for a four-membered ring containing a nitrogen and a carbon with a double bond to oxygen (i.e. a carbonyl)
      to represent the beta-lactam.
   b. For each candidate beta-lactam, check for a five-membered ring that is fused to it (shares at least 2 atoms)
      and that does not contain any sulfur (differentiating it from e.g. penicillins).
3. If such a fused bicyclic system is found, return True along with a reason.
Otherwise, return False with a reason.
"""

from rdkit import Chem

def is_carbapenems(smiles: str):
    """
    Determines if a molecule has the carbapenem bicyclic skeleton based on its SMILES string.
    This searches for the presence of a beta-lactam ring (4-membered ring with one nitrogen and one carbonyl)
    fused to a 5-membered ring that lacks sulfur (differentiating it from e.g. penicillins).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule appears to be a carbapenem, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information as tuples of atom indices.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()

    beta_lactam_candidate = None

    # First loop: identify candidate 4-membered beta-lactam ring(s)
    for ring in rings:
        if len(ring) == 4:
            # Check that the ring contains exactly one nitrogen
            n_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            if n_count != 1:
                continue

            # Check that at least one carbon in the ring is part of a carbonyl (C=O)
            carbonyl_found = False
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 6:  # carbon
                    # Look for a double bond to oxygen
                    for bond in atom.GetBonds():
                        # Check that the bond is double and the neighbor is oxygen.
                        if bond.GetBondType() == Chem.BondType.DOUBLE:
                            nbr = bond.GetOtherAtom(atom)
                            if nbr.GetAtomicNum() == 8:
                                carbonyl_found = True
                                break
                    if carbonyl_found:
                        break
            if not carbonyl_found:
                continue

            # We have a candidate beta-lactam ring; save its indices.
            beta_lactam_candidate = set(ring)
            break  # stop after finding one candidate

    if beta_lactam_candidate is None:
        return False, "No four-membered beta-lactam ring (with N and carbonyl) found"

    # Now look for a five-membered ring that is fused (sharing at least 2 atoms) with the beta-lactam candidate.
    fused_ring_found = False
    for ring in rings:
        if len(ring) == 5:
            ring_set = set(ring)
            # Check for ring fusion: at least 2 atoms are common.
            if len(beta_lactam_candidate.intersection(ring_set)) < 2:
                continue
            # For carbapenems, the 5-membered ring should not contain sulfur (atomic num 16),
            # helping to differentiate from e.g. penicillins.
            if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 16 for idx in ring):
                continue
            # Found a fused 5-membered ring without sulfur.
            fused_ring_found = True
            break

    if not fused_ring_found:
        return False, "No fused five-membered ring (free of sulfur) found with the beta-lactam ring; not a carbapenem skeleton"

    return True, "Found carbapenem bicyclic skeleton (4-membered beta-lactam fused to a sulfur-free 5-membered ring)"

# Example usage (for quick ad-hoc testing):
if __name__ == "__main__":
    test_smiles = "[H][C@@]12CC(=O)N1C=C[C@@H]2C"  # 1beta-methylcarbapenem example
    result, reason = is_carbapenems(test_smiles)
    print(result, reason)