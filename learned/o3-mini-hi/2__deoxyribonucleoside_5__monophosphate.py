"""
Classifies: CHEBI:18241 2'-deoxyribonucleoside 5'-monophosphate
"""
#!/usr/bin/env python3
"""
Classifies chemical entities of the class 
2'-deoxyribonucleoside 5'-monophosphate.
Definition: A 2'-deoxyribonucleoside monophosphate compound with the phosphate group in the 5'-position.
The heuristic:
  (1) Find a furanose candidate ring (5 atoms, exactly 1 oxygen and 4 carbons).
  (2) Verify that one of the ring atoms is attached via an exocyclic bond (optionally via a methylene group) to an oxygen that in turn is attached to a phosphorus atom (phosphate branch).
  (3) Verify that one of the ring atoms is attached to an aromatic nitrogen (as part of a nucleobase).
  (4) Count free –OH groups attached to ring carbons (ignoring those that are part of the phosphate branch); exactly one free –OH is expected (the 3'-OH of deoxyribose).
If all conditions are met, we return True with a reason.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2__deoxyribonucleoside_5__monophosphate(smiles: str):
    """
    Determines if a molecule is a 2'-deoxyribonucleoside 5'-monophosphate based on its SMILES string.
    
    The improved heuristic consists of:
      - Identifying a candidate furanose sugar ring (5 atoms: 4 C, 1 O)
      - Checking for a branch from a ring atom that (directly or via a CH2 group) leads to an oxygen attached to a phosphorus atom.
      - Detecting a nucleobase attachment (heuristically, an aromatic nitrogen neighbor of a ring atom).
      - Counting free hydroxyl (-OH) groups on the ring carbons that are not part of the phosphate branch (expect exactly one).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule fulfills the heuristic, False otherwise.
        str: Explanation of the result.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to better detect -OH groups.
    mol = Chem.AddHs(mol)
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # each ring is a tuple of atom indices
    
    candidate_found = False
    reason = ""
    
    # Loop over all 5-membered rings to find a furanose candidate.
    for ring in rings:
        if len(ring) != 5:
            continue
            
        # Count atoms in the ring.
        count_O = 0
        count_C = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() == "O":
                count_O += 1
            elif atom.GetSymbol() == "C":
                count_C += 1
        if count_O != 1 or count_C != 4:
            continue  # not a typical furanose ring
        
        candidate_ring = set(ring)
        
        # --------- (1) Check for phosphate attachment ----------
        # We examine each atom in the candidate ring.
        phosphate_attached = False
        for idx in candidate_ring:
            atom = mol.GetAtomWithIdx(idx)
            # Look at neighbors not in the ring.
            for nb in atom.GetNeighbors():
                if nb.GetIdx() in candidate_ring:
                    continue
                # Case A: Direct O-attachment
                if nb.GetSymbol() == "O":
                    # Check if this oxygen is attached to a phosphorus atom.
                    for nb2 in nb.GetNeighbors():
                        if nb2.GetIdx() == atom.GetIdx():
                            continue
                        if nb2.GetSymbol() == "P":
                            phosphate_attached = True
                            break
                    if phosphate_attached:
                        break
                # Case B: Attachment via an exocyclic carbon (e.g. CH2 group)
                if nb.GetSymbol() == "C":
                    # Optionally check if this exocyclic carbon appears to be a CH2 (or CHx) by counting attached hydrogens.
                    # Then search its neighbors for an oxygen that is linked to a phosphorus.
                    for nb2 in nb.GetNeighbors():
                        if nb2.GetIdx() == atom.GetIdx():
                            continue
                        if nb2.GetSymbol() == "O":
                            for nb3 in nb2.GetNeighbors():
                                if nb3.GetIdx() in (nb.GetIdx(), nb2.GetIdx()):
                                    continue
                                if nb3.GetSymbol() == "P":
                                    phosphate_attached = True
                                    break
                        if phosphate_attached:
                            break
                    if phosphate_attached:
                        break
            if phosphate_attached:
                break
                
        if not phosphate_attached:
            # This ring candidate does not have a phosphate branch.
            continue
        
        # --------- (2) Check for nucleobase attachment ----------
        # Look for at least one ring atom that is bonded (outside the ring) to an aromatic nitrogen.
        nucleobase_found = False
        for idx in candidate_ring:
            atom = mol.GetAtomWithIdx(idx)
            for nb in atom.GetNeighbors():
                if nb.GetIdx() in candidate_ring:
                    continue
                # Heuristic: if the neighbor is nitrogen and aromatic, it may be part of a nucleobase.
                if nb.GetSymbol() == "N" and nb.GetIsAromatic():
                    nucleobase_found = True
                    break
                # Alternatively, if the neighbor is carbon but is part of an aromatic system containing at least one N,
                # check its neighbors.
                if nb.GetSymbol() == "C" and nb.GetIsAromatic():
                    for sub_nb in nb.GetNeighbors():
                        if sub_nb.GetSymbol() == "N" and sub_nb.GetIsAromatic():
                            nucleobase_found = True
                            break
                if nucleobase_found:
                    break
            if nucleobase_found:
                break
        if not nucleobase_found:
            # No good nucleobase attachment found.
            continue
        
        # --------- (3) Count free hydroxyl (-OH) groups on the ring carbons ----------
        free_oh_count = 0
        for idx in candidate_ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() != "C":
                continue  # we consider -OH only on carbons of the sugar ring
            for nb in atom.GetNeighbors():
                if nb.GetIdx() in candidate_ring:
                    continue
                if nb.GetSymbol() != "O":
                    continue
                # Exclude oxygens that are linked to phosphorus (i.e. part of a phosphate group) 
                # or that lead to the nucleobase (often a direct C-N bond is seen).
                attached_to_P = any(nbb.GetSymbol() == "P" for nbb in nb.GetNeighbors())
                if attached_to_P:
                    continue
                # Check that the oxygen bears at least one hydrogen.
                has_H = any(nbb.GetSymbol() == "H" for nbb in nb.GetNeighbors())
                if has_H:
                    free_oh_count += 1
        # For 2'-deoxyribose (with the 5'-OH phosphorylated and 1'-base attachment) we expect only one free -OH (3'-OH).
        if free_oh_count != 1:
            continue
        
        # All criteria for a deoxyribonucleoside 5'-monophosphate are met.
        candidate_found = True
        reason = ("Found candidate furanose ring (4 C, 1 O) that connects to a phosphate branch "
                  "(directly or via a CH2 group) and a nucleobase (via an aromatic N), "
                  "with exactly one free -OH on the ring (consistent with 2'-deoxyribose).")
        break  # no need to examine other rings
        
    if candidate_found:
        return True, reason
    else:
        return False, "Did not detect a 2'-deoxyribonucleoside 5'-monophosphate fragment based on the heuristic"


# For testing purposes – you can run this module to see outputs for some examples.
if __name__ == "__main__":
    test_smiles = [
        # True positive examples:
        "Nc1nc2n([C@H]3C[C@H](O)[C@@H](COP(O)(O)=O)O3)c(=O)[nH]c2c(=O)[nH]1",  # 8-oxo-dGMP (should be True)
        "Cc1cn([C@H]2C[C@H](O)[C@@H](COP(O)([O-])=O)O2)c(=O)[nH]c1=O",  # dTMP(-)
        # Examples that were previously missed:
        "N1(C2=C(C(N)=NC(=N2)N)N=C1)[C@@H]3O[C@H](COP([O-])([O-])=O)[C@H](C3)O",  # 2-amino-2'-deoxyadenosine 5'-monophosphate(2-)
        "Nc1nc(OCc2ccccc2)c2ncn([C@H]3C[C@H](O)[C@@H](COP(O)(O)=O)O3)c2n1",  # 6-O-benzyl-2'-deoxyguanosine 5'-monophosphate
        "Nc1ncnc2n(cnc12)[C@H]1C[C@H](O)[C@@H](COP(O)(O)=O)O1",  # 2'-deoxyadenosine 5'-monophosphate
        # A false positive example:
        "[C@H]12N3C=4N=C(NC(C4N=C3[C@H]([C@]([C@H](C1)O)(O2)[H])OP(O)(O)=O)=O)N",  # 8,5'-cyclo-2'-deoxyguanosine monophosphate (should be rejected)
    ]
    
    for s in test_smiles:
        res, msg = is_2__deoxyribonucleoside_5__monophosphate(s)
        print(f"SMILES: {s}\nResult: {res}\nReason: {msg}\n")