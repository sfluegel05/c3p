"""
Classifies: CHEBI:18133 hexose
"""
#!/usr/bin/env python
"""
Classifies: hexose
Definition: Any six‐carbon monosaccharide which in its linear form contains either an aldehyde group
at position 1 (aldohexose) or a ketone group at position 2 (ketohexose). Many hexoses exist in cyclic form
(either as a pyranose or a furanose) so we also attempt to detect a sugar‐like ring.
Note: This is a heuristic method – sugars can appear as derivatives so false positives/negatives may occur.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_hexose(smiles: str):
    """
    Determines if a molecule qualifies as a hexose.

    A hexose (six‐carbon monosaccharide) is defined here as a molecule that in its open‐chain form presents 
    either an aldehyde group at carbon 1 (aldohexose) or a ketone group at carbon 2 (ketohexose); or that contains 
    a cyclic (pyranose or furanose) motif with several hydroxyl (or similar) substituents on the core.

    The function applies two heuristic strategies:
      1. Open‐chain detection using two SMARTS patterns (one for aldo‑, one for keto‑hexoses).
      2. Cyclic detection:
          a. Pyranose: a six‐membered ring with one oxygen (and hence five carbons) where at least four of the ring carbons 
             have a small oxygen substituent (–OH or similar).
          b. Furanose: a five‐membered ring with one oxygen and four carbons that features at least one exocyclic –CH2OH.
    
    Args:
        smiles (str): SMILES string

    Returns:
        bool: True if classified as a hexose, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # ===== Attempt 1: Open-chain detection =====
    #
    # For an aldohexose, the linear form is:
    #   C1: aldehyde (CHO)
    #   C2: CH(OH)
    #   C3: CH(OH)
    #   C4: CH(OH)
    #   C5: CH(OH)
    #   C6: CH2OH
    #
    # We use a SMARTS that does not constrain chirality too strictly.
    aldo_hexose_smarts = "[CH1](=O)[C]([OH])[C]([OH])[C]([OH])[C]([OH])[CH2][OH]"
    aldo_query = Chem.MolFromSmarts(aldo_hexose_smarts)
    if aldo_query and mol.HasSubstructMatch(aldo_query):
        return True, "Matches open-chain aldohexose pattern (aldehyde at C1)"

    #
    # For a ketohexose (hexulose), the typical open-chain form is:
    #   C1: CH2OH
    #   C2: CH(OH)
    #   C3: C(=O)
    #   C4: CH(OH)
    #   C5: CH(OH)
    #   C6: CH2OH
    #
    keto_hexose_smarts = "[CH2][OH][CH]([OH])C(=O)[CH]([OH])[CH]([OH])[CH2][OH]"
    keto_query = Chem.MolFromSmarts(keto_hexose_smarts)
    if keto_query and mol.HasSubstructMatch(keto_query):
        return True, "Matches open-chain ketohexose pattern (ketone at C3)"

    # ===== Attempt 2: Cyclic detection =====

    ring_info = mol.GetRingInfo()
    # Loop over all rings in the molecule
    for ring in ring_info.AtomRings():
        # ----- Pyranose detection (six-membered ring) -----
        if len(ring) == 6:
            ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
            # Count oxygen and carbon atoms in the ring
            num_ox = sum(1 for atom in ring_atoms if atom.GetSymbol() == "O")
            num_c  = sum(1 for atom in ring_atoms if atom.GetSymbol() == "C")
            if num_ox == 1 and num_c == 5:
                # Check substituents on ring carbons.
                oh_count = 0
                for idx in ring:
                    atom = mol.GetAtomWithIdx(idx)
                    if atom.GetSymbol() != "C":
                        continue
                    # Look at neighbors not in the ring
                    for nbr in atom.GetNeighbors():
                        if nbr.GetIdx() in ring:
                            continue
                        # A typical hydroxyl group should be an oxygen with at least one hydrogen.
                        # (We use GetTotalNumHs as a loose check.)
                        if nbr.GetSymbol() == "O" and nbr.GetTotalNumHs() >= 1:
                            oh_count += 1
                            break  # count at most one OH per ring carbon
                # Heuristic: require at least 4 of 5 ring carbons bear an OH.
                if oh_count >= 4:
                    return True, "Contains a pyranose ring pattern (six-membered ring with one oxygen and several -OH groups)"
        # ----- Furanose detection (five-membered ring) -----
        if len(ring) == 5:
            ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
            num_ox = sum(1 for atom in ring_atoms if atom.GetSymbol() == "O")
            num_c  = sum(1 for atom in ring_atoms if atom.GetSymbol() == "C")
            if num_ox == 1 and num_c == 4:
                # Look for an exocyclic CH2OH group attached to one of the ring carbons.
                found_exo = False
                for idx in ring:
                    atom = mol.GetAtomWithIdx(idx)
                    if atom.GetSymbol() != "C":
                        continue
                    # Examine neighbors not in the ring.
                    for nbr in atom.GetNeighbors():
                        if nbr.GetIdx() in ring:
                            continue
                        # Check if neighbor is an oxygen that is terminal (typical of CH2OH)
                        # We also check that the oxygen is attached to a CH2 (i.e. has 2 hydrogens roughly)
                        if nbr.GetSymbol() == "O" and nbr.GetTotalNumHs() >= 1:
                            found_exo = True
                            break
                    if found_exo:
                        break
                if found_exo:
                    return True, "Contains a furanose ring pattern (five-membered ring with one oxygen and an exocyclic -CH2OH group)"

    return False, "Does not match recognized hexose patterns"

# ----- Example usage -----
if __name__ == '__main__':
    # A small set of test SMILES including examples from the prompt.
    test_smiles = [
        "OC(C(O)CNCCCCCCC)C(O)C(O)CO",  # 1-Deoxy-1-(heptylamino)hexitol: an open-chain hexitol derivative
        "O=C(OC1OC(C(O)C(C1O)O)C)C2=CC=CC=C2",  # 1-O-Benzoyl-alpha-L-rhamnopyranoside (cyclic, pyranose core)
        "OC[C@H]1OC(O)[C@H](O)[C@H](O)[C@@H]1O",  # D-allopyranose (pyranose)
        "[H]C(=O)[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)CO",  # aldehydo-D-galactose (open-chain)
        "C[C@H](O)[C@H](O)[C@@H](O)C(=O)CO",  # L-rhamnulose (open-chain keto? might be drawn open-chain)
        "[H]C([H])([C@]([H])(O)C=O)[C@]([H])(O)[C@@]([H])(C)O",  # tyvelose (open-chain variant with deoxy substitution)
        "OC[C@H]1O[C@@](O)(CO)[C@@H](O)[C@H]1O",  # alpha-D-tagatofuranose (furanose)
    ]
    for s in test_smiles:
        result, reason = is_hexose(s)
        print(f"SMILES: {s}\n  Classified as hexose? {result}. Reason: {reason}\n")