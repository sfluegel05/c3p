"""
Classifies: CHEBI:18241 2'-deoxyribonucleoside 5'-monophosphate
"""
#!/usr/bin/env python3
"""
Classifies chemical entities of the class 
2'-deoxyribonucleoside 5'-monophosphate.
Definition: A 2'-deoxyribonucleoside monophosphate compound with the phosphate group in the 5'-position.

Improved heuristic:
  1. The molecule must have exactly one phosphorus (P) and a molecular weight below ~500 Da.
  2. Search for a candidate furanose ring (a 5-membered ring with exactly 4 C and 1 O).
  3. In that candidate ring we require two distinct exocyclic branches:
       (a) A nucleobase attachment – one ring carbon must be directly bonded to an aromatic nitrogen.
       (b) A phosphate branch – a different ring carbon must be attached to a CH2-like (non‐ring sp3 carbon with at least 2 H)
           which then is bonded to an oxygen that in turn is attached to a phosphorus.
  4. Count free hydroxyl (-OH) groups on the sugar ring (on ring carbons not involved in the phosphate branch) 
     and expect exactly one free –OH (the 3'-OH in 2'-deoxyribose).
If one candidate ring meets all these conditions, the molecule is classified as a 2'-deoxyribonucleoside 5'-monophosphate.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2__deoxyribonucleoside_5__monophosphate(smiles: str):
    """
    Determines if a molecule is a 2'-deoxyribonucleoside 5'-monophosphate based on its SMILES string.
    
    The improved heuristic consists of:
      - Rejecting molecules that do not have exactly one phosphorus atom or exceed ~500 Da.
      - Scanning its ring systems for a candidate furanose ring (5-membered ring with exactly 4 C and 1 O).
      - In each candidate ring:
          * Identifying one ring carbon that is bonded directly to an aromatic nitrogen (indicative 
            of a nucleobase attachment).
          * Identifying a separate ring carbon that is attached to a non-ring (exocyclic) CH2–like carbon,
            which in turn bonds (directly via an oxygen) to a phosphorus atom (the 5'-phosphate branch).
          * Counting free hydroxyl groups (–OH) on ring carbons (excluding oxygen atoms that are part
            of the phosphate branch); for 2'-deoxyribose we expect exactly one free –OH.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule fulfills the heuristic.
        str: Explanation of the result.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to expose –OH groups.
    mol = Chem.AddHs(mol)
    
    # Global check: exactly one phosphorus and molecular weight below 500 Da.
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count != 1:
        return False, f"Molecule has {p_count} phosphorus atoms; expected exactly one for a monophosphate."
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 500:
        return False, f"Molecular weight too high ({mol_wt:.1f} Da) for a simple 2'-deoxyribonucleoside monophosphate."
    
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # List of atom index tuples for each ring.
    
    # Loop through each candidate 5-membered ring.
    for ring in rings:
        if len(ring) != 5:
            continue
        # Check ring composition: we require exactly 1 oxygen and 4 carbons.
        count_O = 0
        count_C = 0
        for idx in ring:
            sym = mol.GetAtomWithIdx(idx).GetSymbol()
            if sym == "O":
                count_O += 1
            elif sym == "C":
                count_C += 1
        if count_O != 1 or count_C != 4:
            continue  # Not a candidate furanose ring.
        
        candidate_ring = set(ring)
        
        found_nucleobase = False
        found_phosphate_branch = False
        free_oh_count = 0
        
        # We loop over ring atoms (preferably carbons) and inspect their exocyclic neighbors.
        for idx in candidate_ring:
            atom = mol.GetAtomWithIdx(idx)
            # For nucleobase and phosphate branch, we look only at ring carbons.
            if atom.GetSymbol() != "C":
                continue
            
            for nb in atom.GetNeighbors():
                # Only consider neighbors that are outside the ring.
                if nb.GetIdx() in candidate_ring:
                    continue
                nb_sym = nb.GetSymbol()
                
                # (a) Check for nucleobase attachment: immediate neighbor is aromatic nitrogen.
                # Also allow an aromatic carbon that is connected to an aromatic nitrogen.
                if nb_sym == "N" and nb.GetIsAromatic():
                    found_nucleobase = True
                elif nb_sym == "C" and nb.GetIsAromatic():
                    for sub_nb in nb.GetNeighbors():
                        if sub_nb.GetIdx() in candidate_ring:
                            continue
                        if sub_nb.GetSymbol() == "N" and sub_nb.GetIsAromatic():
                            found_nucleobase = True
                            break
                
                # (b) Check for phosphate branch.
                # We require that the branch comes off the ring via a CH2 group (an exocyclic carbon not in any ring,
                # with at least two explicit hydrogens) that in turn is attached (via an oxygen) to a phosphorus.
                if nb_sym == "C" and not nb.IsInRing():
                    # Check if this exocyclic carbon is "CH2‐like" (at least two explicit hydrogens)
                    # We use GetTotalNumHs(includeNeighbors=True) as a rough proxy.
                    if nb.GetTotalNumHs() >= 2:
                        # Look for a neighbor oxygen (not going back to the ring) that subsequently is linked to P.
                        for nb2 in nb.GetNeighbors():
                            if nb2.GetIdx() == atom.GetIdx():
                                continue
                            if nb2.GetSymbol() == "O":
                                # This oxygen should be exocyclic as well.
                                if nb2.IsInRing():
                                    continue
                                # Check if nb2 is bonded to a phosphorus.
                                for nb3 in nb2.GetNeighbors():
                                    if nb3.GetIdx() in (nb.GetIdx(), nb2.GetIdx()):
                                        continue
                                    if nb3.GetSymbol() == "P":
                                        found_phosphate_branch = True
                                        break
                                if found_phosphate_branch:
                                    break
                        if found_phosphate_branch:
                            # Once found on one ring carbon we do not reassign phosphate branch.
                            pass
                # (c) Count free hydroxyl groups.
                # A free -OH is an oxygen (directly attached to a ring carbon) that is not part of the phosphate branch.
                if nb_sym == "O":
                    # If any neighbor of this oxygen is phosphorus, then it is part of the phosphate branch.
                    if any(nbb.GetSymbol() == "P" for nbb in nb.GetNeighbors()):
                        continue
                    # Otherwise, if this oxygen has at least one hydrogen, assume it is a free -OH.
                    if any(nbb.GetSymbol() == "H" for nbb in nb.GetNeighbors()):
                        free_oh_count += 1
                        
        # To be a valid candidate, the ring must have a nucleobase branch, a phosphate branch (via a CH2 group),
        # and exactly one free -OH.
        if found_nucleobase and found_phosphate_branch and free_oh_count == 1:
            reason = ("Found candidate furanose ring (4 C, 1 O) with a nucleobase attachment (an aromatic N on a ring carbon), "
                      "a phosphate branch via an exocyclic CH2 group (leading to an O attached to P), "
                      "and exactly one free -OH on the ring (consistent with 2'-deoxyribose).")
            return True, reason
    
    return False, "Did not detect a valid 2'-deoxyribonucleoside 5'-monophosphate fragment based on the improved heuristic."

# For demonstration/testing:
if __name__ == "__main__":
    # Example SMILES strings (true positives and negatives based on the provided outcomes).
    test_smiles = [
        # True positives:
        "Nc1nc2n([C@H]3C[C@H](O)[C@@H](COP(O)(O)=O)O3)c(=O)[nH]c2c(=O)[nH]1",  # 8-oxo-dGMP
        "Nc1nc(=O)n(cc1CO)[C@H]1C[C@H](O)[C@@H](COP(O)(O)=O)O1",  # 5-hydroxymethyldeoxycytidylic acid
        "Cc1cn([C@H]2C[C@H](O)[C@@H](COP(O)([O-])=O)O2)c(=O)[nH]c1=O",  # dTMP(-)
        "O[C@H]1C[C@@H](O[C@@H]1COP(O)(O)=O)n1cc(F)c(=O)[nH]c1=O",  # 5-fluorodeoxyuridine monophosphate
        "Nc1ccn([C@H]2C[C@H](O)[C@@H](COP(O)(O)=O)O2)c(=O)n1",  # 2'-deoxycytosine 5'-monophosphate
        "Cc1cn([C@H]2C[C@H](O)[C@@H](COP(O)(=O)Oc3ccc(cc3)[N+]([O-])=O)O2)c(=O)[nH]c1=O",  # p-nitrophenyl thymidine 5'-monophosphate
        "N1(C2=C(C(N)=NC(=N2)N)N=C1)[C@@H]3O[C@H](COP([O-])([O-])=O)[C@H](C3)O",  # 2-amino-2'-deoxyadenosine 5'-monophosphate(2-)
        "O[C@H]1C[C@@H](O[C@@H]1COP(O)(O)=O)n1ccc(=O)[nH]c1=O",  # dUMP
        "C1(=C(O)C(NC(N1[C@@H]2O[C@H](COP(O)(=O)O)[C@H](C2)O)=O)=O)O",  # 5,6-dihydroxy-2'-deoxyuridine 5'-monophosphate
        "Nc1ncnc2n(cnc12)[C@H]1C[C@H](O)[C@@H](COP(O)(O)=O)O1",  # 2'-deoxyadenosine 5'-monophosphate
        "Cc1cn([C@@H]2C[C@@H](O)[C@H](COP(O)(O)=O)O2)c(=O)[nH]c1=O",  # dTMP (another representation)
        # Some false positives (should be rejected):
        "N(C1=NC=NC2=C1N=CN2[C@@H]3O[C@H](COP(=O)(O)O)[C@@H](O)C3)CC(=O)N",  # Wrongly classified earlier
    ]
    
    for s in test_smiles:
        res, msg = is_2__deoxyribonucleoside_5__monophosphate(s)
        print(f"SMILES: {s}\nResult: {res}\nReason: {msg}\n")