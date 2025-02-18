"""
Classifies: CHEBI:17334 penicillin
"""
#!/usr/bin/env python
"""
Classifies: Penicillins

Definition:
  “Any member of the group of substituted penams containing two methyl substituents at position 2,
   a carboxylate substituent at position 3 and a carboxamido group at position 6.”

This implementation first removes common counterions (and small uninteresting fragments)
using RDKit’s SaltRemover, then detects a fused bicyclic penicillin core. Finally, it
checks that key substituents (two methyl groups on the carbon adjacent to the core sulfur,
a carboxylate group and a carboxamido group directly attached to the core) are present.

Note: This is a heuristic based on SMARTS and substructure patterns and may not be perfect.
"""

from rdkit import Chem
from rdkit.Chem import SaltRemover

def is_penicillin(smiles: str):
    """
    Determines if a molecule is a penicillin based on its SMILES string.

    Criteria:
      (a) Remove salts/counterions by using RDKit’s SaltRemover and choosing the largest fragment.
      (b) The major fragment must be neutral.
      (c) It must contain a fused bicyclic core:
             – A 4-membered β-lactam ring (with at least one nitrogen and one carbonyl)
             – Fused via exactly two atoms with a 5-membered thiazolidine ring (that contains at least one S)
      (d) On that core, one of the carbons adjacent to the sulfur (position 2) must carry two methyl substituents.
      (e) A carboxylate group (C(=O)[O-] or C(=O)O) must be attached directly to the core.
      (f) A carboxamido fragment (N–C(=O)) must be attached directly to the core.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule appears to be a penicillin; False otherwise.
        str: Explanation of the decision.
    """
    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Remove common salts and small fragments using the default salt remover
    remover = SaltRemover.SaltRemover()
    mol = remover.StripMol(mol)

    # Get separate fragments and select the largest (by heavy atoms)
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return False, "No fragments found in molecule"
    frag = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    if frag.GetNumHeavyAtoms() < 10:
        return False, "Major fragment too small to be a penicillin core"
    if Chem.GetFormalCharge(frag) != 0:
        return False, "Major fragment is charged (even after salt removal)"

    # --- Step 1: Identify the fused penicillin core ---
    # The penicillin core comprises a 4-membered β‑lactam ring (with an N and at least one carbonyl)
    # fused to a 5-membered thiazolidine ring (containing at least one S) sharing exactly 2 atoms.
    def get_penicillin_core(mol):
        ring_info = mol.GetRingInfo()
        rings = ring_info.AtomRings()  # list of tuples of atom indices
        # Search for a candidate 4-membered ring (β-lactam)
        for ring4 in rings:
            if len(ring4) != 4:
                continue
            ring4_set = set(ring4)
            atoms4 = [mol.GetAtomWithIdx(i) for i in ring4_set]
            # Must contain at least one nitrogen
            if not any(atom.GetSymbol() == "N" for atom in atoms4):
                continue
            # At least one carbon should have a double bond to oxygen (carbonyl)
            has_carbonyl = False
            for atom in atoms4:
                if atom.GetSymbol() == "C":
                    for bond in atom.GetBonds():
                        if bond.GetBondTypeAsDouble() == 2.0:
                            other = bond.GetOtherAtom(atom)
                            if other.GetSymbol() == "O":
                                has_carbonyl = True
                                break
                    if has_carbonyl:
                        break
            if not has_carbonyl:
                continue
            # Now find a 5-membered ring (thiazolidine) sharing exactly 2 atoms with the 4-membered ring
            for ring5 in rings:
                if len(ring5) != 5:
                    continue
                common = ring4_set.intersection(ring5)
                if len(common) == 2:
                    atoms5 = [mol.GetAtomWithIdx(i) for i in ring5]
                    if any(atom.GetSymbol() == "S" for atom in atoms5):
                        return ring4_set.union(ring5)
        return None

    core_atoms = get_penicillin_core(frag)
    if core_atoms is None:
        return False, "Molecule does not contain the required fused penicillin core (β‐lactam fused with thiazolidine containing sulfur)"

    # --- Step 2: Check for two methyl substituents on the carbon adjacent to the core sulfur (position 2) ---
    # Loop over S atoms in the core and for each, check its core neighbor carbons.
    found_dimethyl = False
    for atom in frag.GetAtoms():
        if atom.GetIdx() not in core_atoms:
            continue
        if atom.GetSymbol() != "S":
            continue
        # For each S, look at connected carbons that are also in the core.
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() not in core_atoms or nbr.GetSymbol() != "C":
                continue
            methyl_count = 0
            # Count substituents (neighbors not in the core) that are CH3 groups.
            for sub in nbr.GetNeighbors():
                if sub.GetIdx() in core_atoms:
                    continue
                # A methyl group: carbon attached only once (degree 1) and with three hydrogens.
                if sub.GetAtomicNum() == 6 and sub.GetDegree() == 1:
                    methyl_count += 1
            if methyl_count == 2:
                found_dimethyl = True
                break
        if found_dimethyl:
            break
    if not found_dimethyl:
        return False, "Missing two methyl substituents on the core carbon adjacent to the sulfur (expected at position 2)"

    # --- Step 3: Check for a carboxylate substituent attached directly to the core (position 3) ---
    # Look for carboxylate groups of the form C(=O)[O] or C(=O)[O-] that are attached to at least one core atom.
    carboxylate_smarts = Chem.MolFromSmarts("[CX3](=O)[O,OX1-]")
    if carboxylate_smarts is None:
        return False, "Error in carboxylate SMARTS"
    carboxylate_found = False
    for match in frag.GetSubstructMatches(carboxylate_smarts):
        # match[0] is the carboxylate carbon.
        c_idx = match[0]
        for nbr in frag.GetAtomWithIdx(c_idx).GetNeighbors():
            if nbr.GetIdx() in core_atoms:
                carboxylate_found = True
                break
        if carboxylate_found:
            break
    if not carboxylate_found:
        return False, "Missing a carboxylate substituent (C(=O)O or C(=O)[O-]) directly attached to the penicillin core"

    # --- Step 4: Check for a carboxamido fragment (N–C(=O)) attached directly to the core (position 6) ---
    # Use a SMARTS pattern for an amide (non-ring nitrogen) attached to a carbonyl.
    carboxamido_smarts = Chem.MolFromSmarts("[N]-C(=O)")
    if carboxamido_smarts is None:
        return False, "Error in carboxamido SMARTS"
    carboxamido_found = False
    for match in frag.GetSubstructMatches(carboxamido_smarts):
        # At least one atom in the amide fragment must be in the core
        if any(idx in core_atoms for idx in match):
            carboxamido_found = True
            break
    if not carboxamido_found:
        return False, "Missing a carboxamido group (N–C(=O)) directly attached to the penicillin core"

    return True, ("Molecule has a fused penicillin core with dimethyl at position 2, a carboxylate at position 3, "
                  "and a carboxamido at position 6")

# Example usage:
if __name__ == "__main__":
    # Test with Penicillin K SMILES (as one example)
    test_smiles = "CCCCCCCC(=O)N[C@H]1[C@H]2SC(C)(C)[C@@H](N2C1=O)C(O)=O"
    result, reason = is_penicillin(test_smiles)
    print("Is penicillin?", result)
    print("Reason:", reason)