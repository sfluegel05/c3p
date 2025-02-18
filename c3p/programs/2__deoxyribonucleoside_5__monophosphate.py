"""
Classifies: CHEBI:18241 2'-deoxyribonucleoside 5'-monophosphate
"""
#!/usr/bin/env python3
"""
Classifies chemical entities of the class 
2'-deoxyribonucleoside 5'-monophosphate.
Definition: A 2'-deoxyribonucleoside monophosphate compound with the phosphate group in the 5'-position.
The improved heuristic:
  (1) The molecule must have exactly one phosphorus (P) and a molecular weight below ~500 Da.
  (2) Identify a candidate furanose sugar ring (a 5-membered ring with exactly 4 C and 1 O).
  (3) Verify that at least one atom of the ring has a branch (directly or via a CH2 group) that leads to an oxygen
      which is attached to a phosphorus atom.
  (4) Confirm that one of the ring atoms is attached (outside the ring) to an aromatic nitrogen (nucleobase).
  (5) Count free –OH groups on ring carbons (excluding any oxygen involved in the phosphate branch);
      for 2'-deoxyribose the count should be exactly one (the 3'-OH).
If all conditions are met for at least one candidate furanose ring, we return True.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2__deoxyribonucleoside_5__monophosphate(smiles: str):
    """
    Determines if a molecule is a 2'-deoxyribonucleoside 5'-monophosphate based on its SMILES string.
    
    The improved heuristic consists of:
      - Rejecting molecules that do not have exactly one phosphorus atom or have a molecular weight
        outside the expected range for simple nucleotides (< 500 Da).
      - Identifying a candidate furanose ring (5 atoms: 4 C, 1 O).
      - Checking that at least one ring atom has an exocyclic branch that (directly or via a CH2 group)
        leads to an oxygen attached to phosphorus (the phosphate branch).
      - Verifying that one of the ring atoms is attached to an aromatic nitrogen (nucleobase attachment).
      - Counting free hydroxyl (-OH) groups on ring carbons (ignoring oxygens that are part of the phosphate branch)
        and expecting exactly one free -OH (consistent with 2'-deoxyribose).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule fulfills the heuristic, False otherwise.
        str: Explanation of the result.
    """
    # Parse the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that -OH groups are explicit.
    mol = Chem.AddHs(mol)
    
    # Global check: exactly one phosphorus and a typical molecular weight.
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count != 1:
        return False, f"Molecule has {p_count} phosphorus atoms; expected exactly one for a monophosphate."
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 500:
        return False, f"Molecular weight too high ({mol_wt:.1f} Da) for a simple 2'-deoxyribonucleoside monophosphate."
    
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # list of tuples, each tuple is indices of atoms in one ring
    
    # Loop over candidate 5-membered rings.
    for ring in rings:
        if len(ring) != 5:
            continue
            
        # Check ring composition: should have exactly one oxygen and four carbons.
        count_O = 0
        count_C = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            sym = atom.GetSymbol()
            if sym == "O":
                count_O += 1
            elif sym == "C":
                count_C += 1
        if count_O != 1 or count_C != 4:
            continue  # Not a proper candidate furanose ring
        
        candidate_ring = set(ring)
        
        # (1) Check for phosphate branch: Look for a neighbor (or neighbor via a CH2 group)
        # from any ring atom that leads to an oxygen attached to a phosphorus.
        phosphate_attached = False
        for idx in candidate_ring:
            atom = mol.GetAtomWithIdx(idx)
            for nb in atom.GetNeighbors():
                if nb.GetIdx() in candidate_ring:
                    continue
                # Case A: Direct oxygen attachment.
                if nb.GetSymbol() == "O":
                    # Check if this oxygen has a neighbor that is P.
                    for nb2 in nb.GetNeighbors():
                        if nb2.GetIdx() == atom.GetIdx():
                            continue
                        if nb2.GetSymbol() == "P":
                            phosphate_attached = True
                            break
                    if phosphate_attached:
                        break
                # Case B: Attachment via an exocyclic carbon (e.g. CH2 group).
                if nb.GetSymbol() == "C":
                    # Optionally, we might want the exocyclic carbon to have only 1 heavy atom neighbor (a CH2-like pattern).
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
            continue  # This candidate ring does not carry a phosphate branch.
        
        # (2) Check for nucleobase attachment:
        nucleobase_found = False
        for idx in candidate_ring:
            atom = mol.GetAtomWithIdx(idx)
            for nb in atom.GetNeighbors():
                if nb.GetIdx() in candidate_ring:
                    continue
                # If neighbor is nitrogen and aromatic, consider it as part of a nucleobase.
                if nb.GetSymbol() == "N" and nb.GetIsAromatic():
                    nucleobase_found = True
                    break
                # Or if neighbor is a carbon that is part of an aromatic system containing nitrogen.
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
            continue  # No clear nucleobase attachment from this ring.
        
        # (3) Count free hydroxyl (-OH) groups on the ring carbons.
        free_oh_count = 0
        for idx in candidate_ring:
            atom = mol.GetAtomWithIdx(idx)
            # Only consider carbons of the sugar ring.
            if atom.GetSymbol() != "C":
                continue
            for nb in atom.GetNeighbors():
                if nb.GetIdx() in candidate_ring:
                    continue
                if nb.GetSymbol() != "O":
                    continue
                # Exclude oxygens that are attached to phosphorus (i.e. part of the phosphate branch).
                if any(nbb.GetSymbol() == "P" for nbb in nb.GetNeighbors()):
                    continue
                # Require that this oxygen has at least one hydrogen (i.e. it is an -OH group).
                if any(nbb.GetSymbol() == "H" for nbb in nb.GetNeighbors()):
                    free_oh_count += 1
        # For a 2'-deoxyribose, we expect exactly one free hydroxyl (located at 3'-position).
        if free_oh_count != 1:
            continue
        
        # If all criteria are met, we have found a valid candidate.
        reason = ("Found candidate furanose ring (4 C, 1 O) with a phosphate branch (directly or via a CH2 group) "
                  "and a nucleobase attachment (aromatic N), and exactly one free -OH on the ring "
                  "(consistent with 2'-deoxyribose).")
        return True, reason

    # If no candidate ring meets our heuristic:
    return False, "Did not detect a valid 2'-deoxyribonucleoside 5'-monophosphate fragment based on the heuristic."

# For testing purposes, run this module to examine some example SMILES.
if __name__ == "__main__":
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
        # False positives (should be rejected):
        "S(C(=O)C(CCCC(CCCC(C)C)C)C)CCNC(=O)CCNC(=O)C(O)C(COP(OP(OCC1OC(N2C3=NC=NC(N)=C3N=C2)C(O)C1OP(O)(O)=O)(O)=O)(O)=O)(C)C",  # a long-chain CoA derivative
    ]
    
    for s in test_smiles:
        res, msg = is_2__deoxyribonucleoside_5__monophosphate(s)
        print(f"SMILES: {s}\nResult: {res}\nReason: {msg}\n")