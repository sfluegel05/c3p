"""
Classifies: CHEBI:46633 carbapenems
"""
"""
Classifies: Carbapenems (beta-lactam antibiotics that have a carbapenem skeleton,
which is a four-membered beta-lactam ring fused to a five-membered ring lacking sulfur
and containing a C=C bond in the 5-membered ring).

The approach:
1. Parse the SMILES.
2. Identify a candidate beta-lactam ring (a ring of size 4 that contains exactly one nitrogen
   and at least one carbon (C) that is double-bonded to an oxygen (C=O)).
3. Look for a fused five-membered ring (ring of size 5) sharing at least 2 atoms with the beta-lactam;
   the five-membered ring must not contain any sulfur and must contain at least one C=C double bond.
4. If both are found, return True and a reason; otherwise, return False.
"""

from rdkit import Chem

def is_carbapenems(smiles: str):
    """
    Determines if a molecule has the carbapenem bicyclic skeleton based on its SMILES string.
    It looks for a candidate beta-lactam ring (4-membered with one nitrogen and a carbonyl) fused
    to a 5-membered ring that (a) is free of sulfur and (b) contains at least one carbon-carbon double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule appears to be a carbapenem, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES string to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()

    beta_lactam_candidate = None

    # Identify a candidate beta-lactam ring: 4-membered, exactly one nitrogen,
    # and at least one carbon (C) in the ring which is double-bonded to an oxygen.
    for ring in rings:
        if len(ring) != 4:
            continue
        # Check for exactly one nitrogen in the ring.
        n_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
        if n_count != 1:
            continue

        carbonyl_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Check only carbon atoms
            if atom.GetAtomicNum() == 6:
                # Look through bonds for a double bond to oxygen.
                for bond in atom.GetBonds():
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        nbr = bond.GetOtherAtom(atom)
                        if nbr.GetAtomicNum() == 8:
                            carbonyl_found = True
                            break
                if carbonyl_found:
                    break
        if not carbonyl_found:
            continue

        # Candidate beta-lactam ring found.
        beta_lactam_candidate = set(ring)
        break

    if beta_lactam_candidate is None:
        return False, "No four-membered beta-lactam ring (with one nitrogen and carbonyl) found"

    # Look for a fused 5-membered ring that meets additional criteria.
    fused_ring_found = False
    for ring in rings:
        if len(ring) != 5:
            continue
        ring_set = set(ring)
        # Check that the ring and the beta-lactam candidate share at least 2 atoms (fused rings)
        if len(beta_lactam_candidate.intersection(ring_set)) < 2:
            continue
        # Exclude rings containing sulfur (atomic number 16).
        if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 16 for idx in ring):
            continue

        # Additionally, require that within this 5-membered ring there is a carbon-carbon double bond.
        # Iterate over bonds whose both atoms are in the ring.
        has_CC_double = False
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetIdx() in ring_set and a2.GetIdx() in ring_set:
                # Check if both are carbons and bond is double.
                if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6 and bond.GetBondType() == Chem.BondType.DOUBLE:
                    has_CC_double = True
                    break
        if not has_CC_double:
            continue

        # Found a fused 5-membered ring fulfilling all conditions.
        fused_ring_found = True
        break

    if not fused_ring_found:
        return False, "No fused five-membered ring (free of sulfur and containing a C=C bond) found with the beta-lactam ring; not a carbapenem skeleton"

    return True, "Found carbapenem bicyclic skeleton (4-membered beta-lactam fused to a sulfur-free 5-membered ring with C=C bond)"

# Example usage (for quick ad-hoc testing)
if __name__ == "__main__":
    # 1beta-methylcarbapenem example
    test_smiles = "[H][C@@]12CC(=O)N1C=C[C@@H]2C"
    result, reason = is_carbapenems(test_smiles)
    print(result, reason)