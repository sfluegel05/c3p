"""
Classifies: CHEBI:18133 hexose
"""
#!/usr/bin/env python
"""
Classifies: hexose
Definition: Any six‐carbon monosaccharide which in its linear form contains either an aldehyde group
at position 1 (aldohexose) or a ketone group at position 2 (ketohexose). Many hexoses exist in cyclic form
(either as a pyranose or a furanose) so we also attempt to detect a sugar‐like ring.
Note: Heuristics are used. Many modifications (e.g. glycosidic linkages) may affect the match.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_hexose(smiles: str):
    """
    Determines if a molecule qualifies as a hexose.

    We use three heuristic strategies:
      1. Reject molecules with phosphorus atoms (commonly signaling nucleotides or related structures).
      2. Open–chain detection using loose SMARTS patterns for an aldohexose (terminal aldehyde)
         or ketohexose (internal ketone).
      3. Cyclic (ring) detection for sugar–like motifs:
            a. Pyranose candidate: a six–membered ring with exactly one oxygen. Of the five ring–carbons,
               at least 3 should carry an –OH (or similar) group.
            b. Furanose candidate: a five–membered ring with exactly one oxygen and at least one exocyclic –CH2OH group.
         In either case the candidate substructure should not be a very minor portion of a large molecule.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a hexose, otherwise False.
        str: Explanation for the classification decision.
    """
    # Parse SMILES and add explicit hydrogens.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # ----- Step 1: Quick rejection (e.g. nucleotides have phosphorus) -----
    if any(atom.GetSymbol() == "P" for atom in mol.GetAtoms()):
        return False, "Contains phosphorus which indicates a nucleotide or related structure"
    
    # ----- Step 2: Open-chain detection via SMARTS -----
    # We use a loose SMARTS for a six–carbon chain.
    # For an aldohexose the open-chain (linear) molecule should have a terminal aldehyde group.
    # The pattern below does not require all –OH groups to be unmodified.
    aldo_smarts = "[H][C;X3](=O)[C;X4][C;X4][C;X4][C;X4][CH2][OX2H]"
    aldo_query = Chem.MolFromSmarts(aldo_smarts)
    if aldo_query and mol.HasSubstructMatch(aldo_query):
        return True, "Matches open-chain aldohexose pattern (aldehyde at C1)"
    
    # For a ketohexose we expect the carbonyl to be internal.
    keto_smarts = "[CH2][OX2H][C;X4][C;X4]C(=O)[C;X4][CH2][OX2H]"
    keto_query = Chem.MolFromSmarts(keto_smarts)
    if keto_query and mol.HasSubstructMatch(keto_query):
        return True, "Matches open-chain ketohexose pattern (ketone in-chain)"
    
    # ----- Step 3: Cyclic (ring) detection for sugar rings -----
    ring_info = mol.GetRingInfo()
    total_atoms = mol.GetNumAtoms()
    candidate_found = False
    
    for ring in ring_info.AtomRings():
        ring_size = len(ring)
        # Get the atoms in the ring.
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        
        # Pyranose candidate: 6-membered ring with exactly one oxygen.
        if ring_size == 6:
            num_ox_in_ring = sum(1 for atom in ring_atoms if atom.GetSymbol() == "O")
            if num_ox_in_ring != 1:
                continue
            # Count ring-carbon substituents that are -OH.
            oh_count = 0
            ch2oh_found = False
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                # Consider only carbons in the ring.
                if atom.GetSymbol() != "C":
                    continue
                # Look at neighbors that are not in the ring.
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() in ring:
                        continue
                    # Look for hydroxyl: oxygen atom with at least one hydrogen.
                    if nbr.GetSymbol() == "O" and nbr.GetTotalNumHs() >= 1:
                        oh_count += 1
                        # Additionally, check if this oxygen is attached to a CH2 group (exocyclic –CH2OH)
                        for o_nbr in nbr.GetNeighbors():
                            if o_nbr.GetIdx() == atom.GetIdx():
                                continue
                            if o_nbr.GetSymbol() == "C" and o_nbr.GetTotalNumHs() >= 2:
                                ch2oh_found = True
                        break  # count only one substituent per ring carbon
            # Heuristic: require at least 3 of 5 ring carbons bear –OH.
            if oh_count >= 3:
                # Also, the sugar part should be a fair fraction of the whole molecule.
                if len(ring) < total_atoms * 0.8:  # if the ring is not just a minor decoration
                    candidate_found = True
                    return True, "Contains a pyranose ring pattern (6-membered ring with one oxygen and multiple -OH substituents)"
        
        # Furanose candidate: 5-membered ring with exactly one oxygen.
        if ring_size == 5:
            num_ox_in_ring = sum(1 for atom in ring_atoms if atom.GetSymbol() == "O")
            if num_ox_in_ring != 1:
                continue
            exo_ch2oh = False
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetSymbol() != "C":
                    continue
                # Check neighbors outside ring for a -CH2OH fragment.
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() in ring:
                        continue
                    if nbr.GetSymbol() == "C" and nbr.GetTotalNumHs() >= 2:
                        # Look for attached oxygen with at least one hydrogen.
                        for subnbr in nbr.GetNeighbors():
                            if subnbr.GetIdx() == atom.GetIdx():
                                continue
                            if subnbr.GetSymbol() == "O" and subnbr.GetTotalNumHs() >= 1:
                                exo_ch2oh = True
                                break
                    if exo_ch2oh:
                        break
                if exo_ch2oh:
                    break
            if exo_ch2oh:
                if len(ring) < total_atoms * 0.8:
                    candidate_found = True
                    return True, "Contains a furanose ring pattern (5-membered ring with one oxygen and an exocyclic -CH2OH group)"
    
    # If no strategy succeeded then return negative.
    if not candidate_found:
        return False, "Does not match recognized hexose patterns"

    # Fallback (should not reach here)
    return False, "Failed to classify molecule"

# ----- Example usage -----
if __name__ == '__main__':
    # Test a few examples (these include examples from true positives and known problematic cases)
    test_smiles = [
        "OC(C(O)CNCCCCCCC)C(O)C(O)CO",                         # 1-Deoxy-1-(heptylamino)hexitol (open-chain derivative)
        "O=C(OC1OC(C(O)C(C1O)O)C)C2=CC=CC=C2",                 # 1-O-Benzoyl-alpha-L-rhamnopyranoside (cyclic but with extra acyl group)
        "OC[C@H]1OC(O)[C@H](O)[C@H](O)[C@@H]1O",               # D-allopyranose
        "[H]C(=O)[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)CO",         # aldehydo-D-galactose (open-chain)
        "C[C@H](O)[C@H](O)[C@@H](O)C(=O)CO",                   # L-rhamnulose (open-chain keto variant)
        "[H]C([H])([C@]([H])(O)C=O)[C@]([H])(O)[C@@]([H])(C)O",  # tyvelose (open-chain with deoxy substitution)
        "OC[C@H]1O[C@@](O)(CO)[C@@H](O)[C@H]1O",               # alpha-D-tagatofuranose (cyclic furanose)
    ]
    for s in test_smiles:
        result, reason = is_hexose(s)
        print(f"SMILES: {s}\n  Classified as hexose? {result}. Reason: {reason}\n")