"""
Classifies: CHEBI:18241 2'-deoxyribonucleoside 5'-monophosphate
"""
"""
Classifies chemical entities of the class 
2'-deoxyribonucleoside 5'-monophosphate.
Definition: A 2'-deoxyribonucleoside monophosphate compound with the phosphate group in the 5'-position.
Improved heuristic:
  (1) Look for a 5-membered ring (furanose) containing 4 carbons and 1 oxygen.
  (2) Verify that an exocyclic carbon attached to the ring via a single bond leads (through an oxygen) to a phosphorus atom.
  (3) Check that a ring carbon is attached (outside the ring) to an aromatic heterocycle (nucleobase).
  (4) Count free hydroxyl (-OH) groups attached to ring carbons (ignoring those in a phosphate link).
      In 2'-deoxyribose (with the 5'-position phosphorylated) only one free –OH (at the 3'-position) should be present.
If all conditions are met, we classify the molecule as a 2'-deoxyribonucleoside 5'-monophosphate.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2__deoxyribonucleoside_5__monophosphate(smiles: str):
    """
    Determines if a molecule is a 2'-deoxyribonucleoside 5'-monophosphate based on its SMILES string.
    
    The improved heuristic consists of:
      - Identifying a candidate furanose sugar ring (5 atoms: 4 C and 1 O)
      - Verifying that a phosphate group (P, attached via an oxygen on an exocyclic carbon) is present
      - Verifying that the sugar is coupled to a nucleobase [heuristically, one of the ring carbons links to 
        an external aromatic heterocycle (ring size > 5 and containing at least one N)]
      - Checking that exactly one free hydroxyl (-OH) is attached to the ring carbons (free = not part of a phosphate link)
        (Note: In ribonucleotides one expects two free –OH; deoxyribose loses the 2'-OH.)
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets these features, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens for correct -OH detection.
    mol = Chem.AddHs(mol)
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # tuple of tuples of atom indices
    
    found_candidate = False
    reason = ""
    
    # We loop over all rings and keep track if one matches the deoxy sugar criteria.
    for ring in rings:
        if len(ring) != 5:
            continue
        
        # Count ring elements.
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        count_O = sum(1 for a in ring_atoms if a.GetSymbol() == "O")
        count_C = sum(1 for a in ring_atoms if a.GetSymbol() == "C")
        if count_O != 1 or count_C != 4:
            continue  # Not a typical furanose ring
        
        # Now check for phosphate attachment.
        phosphate_attached = False
        # We search over ring atoms: look for a neighbor not in the ring that is a carbon 
        # and via an oxygen is connected to a phosphorus.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nb in atom.GetNeighbors():
                if nb.GetIdx() in ring:
                    continue
                if nb.GetSymbol() != "C":
                    continue
                # nb is an exocyclic carbon; check its neighbors for an oxygen leading to P.
                for nb2 in nb.GetNeighbors():
                    if nb2.GetIdx() == atom.GetIdx():
                        continue
                    if nb2.GetSymbol() != "O":
                        continue
                    # Check if this oxygen is connected to a phosphorus atom.
                    for nb3 in nb2.GetNeighbors():
                        if nb3.GetSymbol() == "P":
                            # Found a phosphate attachment.
                            phosphate_attached = True
                            break
                    if phosphate_attached:
                        break
                if phosphate_attached:
                    break
            if phosphate_attached:
                break
        if not phosphate_attached:
            continue  # This ring does not lead to a phosphate group
        
        # Check for nucleobase attachment.
        has_nucleobase = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nb in atom.GetNeighbors():
                if nb.GetIdx() in ring:
                    continue
                # Heuristic: if a neighbor belongs to a ring of size > 5 and is aromatic and contains at least one N, we assume it is part of a nucleobase.
                for r in ring_info.AtomRings():
                    if nb.GetIdx() in r and len(r) > 5:
                        # Check that the neighbor has at least one nitrogen (or is itself N).
                        if nb.GetSymbol() == "N":
                            has_nucleobase = True
                            break
                        else:
                            # Alternatively, if any neighbor of nb is N (and nb is aromatic), use that.
                            if nb.GetIsAromatic():
                                for nb_sub in nb.GetNeighbors():
                                    if nb_sub.GetSymbol() == "N":
                                        has_nucleobase = True
                                        break
                        if has_nucleobase:
                            break
                if has_nucleobase:
                    break
            if has_nucleobase:
                break
        if not has_nucleobase:
            continue  # This sugar ring is not attached to a nucleobase
        
        # Count free hydroxyl groups (–OH) on ring carbons.
        # We want to count those oxygens that are directly attached to a ring carbon 
        # but which are not part of the phosphate link (i.e. oxygen that does NOT also neighbor a P).
        free_OH_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() != "C":
                continue
            for nb in atom.GetNeighbors():
                if nb.GetIdx() in ring:
                    continue
                if nb.GetSymbol() != "O":
                    continue
                # Check: if this oxygen is connected to a phosphorus, then it is in the phosphate linkage.
                attached_to_P = any(nbb.GetSymbol() == "P" for nbb in nb.GetNeighbors())
                if attached_to_P:
                    continue
                # Check that this oxygen has at least one hydrogen (i.e. is -OH).
                has_H = any(nbb.GetSymbol() == "H" for nbb in nb.GetNeighbors())
                if has_H:
                    free_OH_count += 1
        # In a 2'-deoxyribose with the 5'-OH replaced by phosphate and 1'-nucleobase,
        # normally only one free –OH (at the 3'-position) is present on the ring.
        if free_OH_count != 1:
            # We might be in the ribonucleotide case (which has 2 free OH groups) or mis-assigned.
            continue
        
        # All checks passed.
        found_candidate = True
        reason = ("Found furanose ring (4 C, 1 O) attached via an exocyclic carbon to a phosphate, "
                  "with a nucleobase substitution and exactly one free -OH on the ring (consistent with 2'-deoxyribose)")
        break

    if found_candidate:
        return True, reason
    else:
        return False, "Did not find the required 2'-deoxyribonucleoside 5'-monophosphate fragment"


# For testing purposes – you can run this module to see outputs for some examples.
if __name__ == "__main__":
    test_smiles = [
        # True examples
        "Nc1nc2n([C@H]3C[C@H](O)[C@@H](COP(O)(O)=O)O3)c(=O)[nH]c2c(=O)[nH]1",  # 8-oxo-dGMP
        "Cc1cn([C@H]2C[C@H](O)[C@@H](COP(O)([O-])=O)O2)c(=O)[nH]c1=O",  # dTMP(-)
        # A couple of false positive / negative test cases could be added:
        "NCCCCCCNc1ncnc2n(cnc12)[C@@H]1O[C@@H]2COP(O)(=O)O[C@H]2[C@H]1O",  # N(6)-(6-aminohexyl)-cAMP, should be rejected (ribose has 2 free OHs)
    ]
    
    for s in test_smiles:
        res, msg = is_2__deoxyribonucleoside_5__monophosphate(s)
        print(f"SMILES: {s}\nResult: {res}\nReason: {msg}\n")