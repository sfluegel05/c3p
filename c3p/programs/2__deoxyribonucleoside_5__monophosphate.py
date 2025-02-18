"""
Classifies: CHEBI:18241 2'-deoxyribonucleoside 5'-monophosphate
"""
#!/usr/bin/env python3
"""
Classifies: CHEBI:2'-deoxyribonucleoside 5'-monophosphate

Definition: A 2'-deoxyribonucleoside monophosphate compound with the phosphate group in the 5'-position.
Heuristic (improved):
  1. The molecule must have exactly one phosphorus (P) and a molecular weight below ~500 Da.
  2. Look for a candidate furanose ring (a 5-membered ring with exactly 4 carbons and 1 oxygen).
  3. In that candidate ring:
       (a) There must be a nucleobase branch – a ring carbon having at least one exocyclic neighbor
           that is either an aromatic nitrogen or a carbon that itself is bonded to an aromatic nitrogen.
       (b) There must be a phosphate branch – a ring carbon attached to an exocyclic CH2 group (with ≥2 H)
           that in turn is bonded (via an oxygen) to a phosphorus atom.
       (c) Count free hydroxyl (-OH) groups on ring carbons (exocyclic oxygen atoms not involved in the phosphate branch)
           and expect exactly one free -OH (the 3'-OH in deoxyribose).
If one candidate ring meets all these conditions, the molecule is classified as a 2'-deoxyribonucleoside 5'-monophosphate.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2__deoxyribonucleoside_5__monophosphate(smiles: str):
    """
    Determines if a molecule is a 2'-deoxyribonucleoside 5'-monophosphate based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule fulfills the heuristic.
        str: Explanation of the classification.
    """
    # Parse SMILES and add explicit hydrogens to expose all –OH groups.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Global check: exactly one phosphorus and molecular weight below 500 Da.
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count != 1:
        return False, f"Molecule has {p_count} phosphorus atoms; expected exactly one for a monophosphate."
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 500:
        return False, f"Molecular weight too high ({mol_wt:.1f} Da) for a simple 2'-deoxyribonucleoside monophosphate."
        
    # Get rings available in the molecule.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # List of tuples (atom indices) for each ring.
    
    # We'll search for a candidate furanose ring: 5-membered with exactly 4 carbons and 1 oxygen.
    for ring in rings:
        if len(ring) != 5:
            continue
        count_C = 0
        count_O = 0
        for idx in ring:
            sym = mol.GetAtomWithIdx(idx).GetSymbol()
            if sym == "C":
                count_C += 1
            elif sym == "O":
                count_O += 1
        if count_C != 4 or count_O != 1:
            continue  # Not a proper furanose candidate.
        
        candidate_ring = set(ring)
        found_nucleobase = False
        found_phosphate_branch = False
        free_oh_count = 0
        
        # We will also record oxygen atoms that are part of a phosphate branch
        phosphate_oxy_ids = set()
        
        # For each atom in the ring, check its neighbors outside the ring.
        for idx in candidate_ring:
            atom = mol.GetAtomWithIdx(idx)
            # We focus on ring carbons for the branching attachments.
            if atom.GetSymbol() != "C":
                continue
            
            for nb in atom.GetNeighbors():
                if nb.GetIdx() in candidate_ring:
                    continue  # Skip atoms inside the ring.
                nb_sym = nb.GetSymbol()
                
                # (a) Identify nucleobase branch:
                # If the neighbor is an aromatic nitrogen, we mark nucleobase branch found.
                if nb_sym == "N" and nb.GetIsAromatic():
                    found_nucleobase = True
                # Alternatively, if neighbor is a carbon and (directly or via one bond)
                # connects to an aromatic nitrogen.
                elif nb_sym == "C":
                    # Check its neighbors (outside the ring) for an aromatic N.
                    for sub_nb in nb.GetNeighbors():
                        if sub_nb.GetIdx() in candidate_ring:
                            continue
                        if sub_nb.GetSymbol() == "N" and sub_nb.GetIsAromatic():
                            found_nucleobase = True
                            break
                
                # (b) Identify phosphate branch:
                # Look for a CH2-like (exocyclic) carbon: an atom that is not in the ring with at least 2 attached hydrogens.
                if nb_sym == "C" and (not nb.IsInRing()):
                    # Check if it has at least two hydrogens (explicit or implicit)
                    if nb.GetTotalNumHs() >= 2:
                        # Now look for an exocyclic oxygen attached to this carbon that in turn is bonded to a phosphorus.
                        for nb2 in nb.GetNeighbors():
                            if nb2.GetIdx() == atom.GetIdx():
                                continue  # Skip the bond back to the ring.
                            if nb2.GetSymbol() == "O" and (not nb2.IsInRing()):
                                # Check that this oxygen is not already counted as free –OH.
                                # Now check if nb2 connects to a phosphorus.
                                for nb3 in nb2.GetNeighbors():
                                    if nb3.GetIdx() in (nb.GetIdx(), nb2.GetIdx()):
                                        continue
                                    if nb3.GetSymbol() == "P":
                                        found_phosphate_branch = True
                                        phosphate_oxy_ids.add(nb2.GetIdx())
                                        break
                                if found_phosphate_branch:
                                    break
                    if found_phosphate_branch:
                        # We already found a valid phosphate branch from one of the ring carbons.
                        pass
                        
                # (c) Count free –OH groups:
                # A free –OH group comes from an oxygen not already used in the phosphate branch attached to a ring carbon.
                if nb_sym == "O":
                    # If this oxygen is part of the phosphate branch, skip counting.
                    if nb.GetIdx() in phosphate_oxy_ids:
                        continue
                    # Consider it as free –OH if it is directly attached to at least one hydrogen.
                    has_H = False
                    for sub in nb.GetNeighbors():
                        if sub.GetSymbol() == "H":
                            has_H = True
                            break
                    if has_H:
                        free_oh_count += 1
        
        # For a valid deoxyribonucleoside monophosphate, we expect:
        # - a found nucleobase branch,
        # - a found phosphate branch,
        # - exactly one free –OH on the sugar ring (i.e. the 3'-OH; the 2'-position should be deoxy).
        if found_nucleobase and found_phosphate_branch and free_oh_count == 1:
            reason = ("Candidate furanose ring (4 C, 1 O) detected with a nucleobase branch "
                      "(exocyclic aromatic N or nearby aromatic feature) and a phosphate branch "
                      "(exocyclic CH2 linked via O to P), with exactly one free –OH on the ring.")
            return True, reason

    return False, "Did not detect a valid 2'-deoxyribonucleoside 5'-monophosphate fragment based on the improved heuristic."

# For demonstration/testing:
if __name__ == "__main__":
    test_smiles = [
        "Nc1nc2n([C@H]3C[C@H](O)[C@@H](COP(O)(O)=O)O3)c(=O)[nH]c2c(=O)[nH]1",  # 8-oxo-dGMP
        "Nc1nc(=O)n(cc1CO)[C@H]1C[C@H](O)[C@@H](COP(O)(O)=O)O1",  # 5-hydroxymethyldeoxycytidylic acid
        "C12N(C(NC(C1(C)C3(C2N(C(NC3=O)=O)[C@@H]4O[C@@H]([C@H](C4)O)COP(=O)(O)O)C)=O)=O)[C@@H]5O[C@@H]([C@H](C5)O)COP(=O)(O)O",  # thymidine 5'-monophosphate dimer (should be rejected: 2 P atoms)
        "Cc1cn([C@H]2C[C@H](O)[C@@H](COP(O)([O-])=O)O2)c(=O)[nH]c1=O",  # dTMP(-)
        "O[C@H]1C[C@@H](O[C@@H]1COP(O)(O)=O)n1cc(F)c(=O)[nH]c1=O",  # 5-fluorodeoxyuridine monophosphate
        "Nc1ccn([C@H]2C[C@H](O)[C@@H](COP(O)(O)=O)O2)c(=O)n1",  # 2'-deoxycytosine 5'-monophosphate
        "Cc1cn([C@H]2C[C@H](O)[C@@H](COP(O)(=O)Oc3ccc(cc3)[N+]([O-])=O)O2)c(=O)[nH]c1=O",  # p-nitrophenyl thymidine 5'-monophosphate
        "N1(C2=C(C(N)=NC(=N2)N)N=C1)[C@@H]3O[C@H](COP([O-])([O-])=O)[C@H](C3)O",  # 2-amino-2'-deoxyadenosine 5'-monophosphate(2-)
        "O[C@H]1C[C@@H](O[C@@H]1COP(O)(O)=O)n1ccc(=O)[nH]c1=O",  # dUMP
        "C1(=C(O)C(NC(N1[C@@H]2O[C@H](COP(O)(=O)O)[C@H](C2)O)=O)=O)O",  # 5,6-dihydroxy-2'-deoxyuridine 5'-monophosphate
        "Nc1ncnc2n(cnc12)[C@H]1C[C@H](O)[C@@H](COP(O)(O)=O)O1",  # 2'-deoxyadenosine 5'-monophosphate
        "Cc1cn([C@@H]2C[C@@H](O)[C@H](COP(O)(O)=O)O2)c(=O)[nH]c1=O",  # dTMP (another representation)
    ]
    
    for s in test_smiles:
        valid, explanation = is_2__deoxyribonucleoside_5__monophosphate(s)
        print(f"SMILES: {s}\nResult: {valid}\nReason: {explanation}\n")