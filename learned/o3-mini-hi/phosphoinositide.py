"""
Classifies: CHEBI:18179 phosphoinositide
"""
#!/usr/bin/env python
"""
Classifies: Phosphoinositide
Definition: Any phosphatidylinositol that is phosphorylated at one or more of the hydroxy groups of inositol.
A phosphoinositide must contain a glycerophosphatidylinositol “core” featuring two acyl chains attached via ester bonds,
and a phosphate group bridging the glycerol backbone to an inositol ring (a cyclohexane bearing several hydroxyl groups)
where at least one hydroxyl is phosphorylated.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    
    Heuristic criteria:
      1. The input SMILES must be valid.
      2. The molecule must contain at least two phosphorus atoms.
      3. The molecule must contain at least two acyl ester groups (C(=O)O motif) — suggesting the two acyl chains.
      4. The molecule must possess a phosphorylated inositol headgroup.
         This is verified by finding a neutral phosphate (P atom with formal charge 0) that is 
         connected via a single-bonded oxygen to a carbon. That carbon must belong to a six-membered ring
         (i.e. cyclohexane) made entirely of carbons. Finally, we count the number of hydroxyl (-OH) substituents
         on the ring atoms, requiring at least 3 to favor an inositol-like pattern.
         
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the criteria for a phosphoinositide, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # To better examine substituents (such as hydroxyl groups), add explicit hydrogens.
    mol = Chem.AddHs(mol)
    
    # Criterion 1: Check for at least two phosphorus atoms (atomic number 15).
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count < 2:
        return False, f"Only {p_count} phosphorus atom(s) found; phosphoinositides require at least 2."
    
    # Criterion 2: Look for at least two acyl ester groups ("C(=O)O").
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found only {len(ester_matches)} acyl ester group(s); need at least 2 acyl chains."
    
    # Criterion 3: Identify the phosphorylated inositol head.
    # We look for a neutral phosphate group (P with formal charge 0) that is connected via a single-bonded O 
    # to a carbon that is a member of a six-membered (cyclohexane) ring comprised solely of carbons.
    # Then we count the number of hydroxyl (-OH) substituents on that ring.
    found_inositol = False
    for p_atom in mol.GetAtoms():
        if p_atom.GetAtomicNum() == 15 and p_atom.GetFormalCharge() == 0:
            # Examine oxygen neighbors of the phosphate atom.
            for o_neigh in p_atom.GetNeighbors():
                if o_neigh.GetAtomicNum() != 8:
                    continue
                # Check that the bond between phosphorus and oxygen is a single bond.
                bond = mol.GetBondBetweenAtoms(p_atom.GetIdx(), o_neigh.GetIdx())
                if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                    continue
                # Now, from this oxygen, check its other neighbor(s) (other than the phosphorus).
                for carbon in o_neigh.GetNeighbors():
                    if carbon.GetIdx() == p_atom.GetIdx():
                        continue
                    if carbon.GetAtomicNum() != 6:
                        continue
                    # Check if this carbon is part of a six-membered ring.
                    ring_info = mol.GetRingInfo()
                    for ring in ring_info.AtomRings():
                        if len(ring) == 6 and carbon.GetIdx() in ring:
                            # Ensure the entire ring consists of carbons.
                            if not all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
                                continue
                            # Now count hydroxyl (-OH) groups attached to ring atoms.
                            oh_count = 0
                            for idx in ring:
                                ring_atom = mol.GetAtomWithIdx(idx)
                                for nbr in ring_atom.GetNeighbors():
                                    # Look for oxygen neighbors that (likely) represent hydroxyl groups.
                                    if nbr.GetAtomicNum() == 8:
                                        # We check if there is at least one hydrogen attached to this oxygen.
                                        # (We already added explicit hydrogens.)
                                        if any(nei.GetAtomicNum() == 1 for nei in nbr.GetNeighbors()):
                                            oh_count += 1
                                            break  # count one per ring atom
                            if oh_count >= 3:
                                found_inositol = True
                                break
                    if found_inositol:
                        break
                if found_inositol:
                    break
        if found_inositol:
            break

    if not found_inositol:
        return False, "No phosphorylated inositol head detected (missing inositol ring with sufficient -OH substituents connecting to a neutral phosphate)."
    
    return True, ("Molecule contains at least two phosphorus atoms, at least two acyl ester groups, and a phosphorylated inositol head "
                  "as indicated by a neutral phosphate bridging to a cyclohexane ring with multiple hydroxyl substituents.")

# Debug/test code (only executes when module is run directly)
if __name__ == "__main__":
    test_smiles = [
        # Example: PIP(18:0/16:0)
        "[C@@H]1(C(C(C([C@H](C1O)OP(=O)(O)O)O)O)O)OP(OC[C@](COC(CCCCCCCCCCCCCCCCC)=O)([H])OC(CCCCCCCCCCCCCCC)=O)(O)=O",
        # Example: 1-stearoyl-2-oleoyl-sn-glycero-3-phospho-1D-myo-inositol 4-phosphate
        "CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H](O)[C@H]1O)OC(=O)CCCCCCC\\C=C/CCCCCCCC"
    ]
    for smi in test_smiles:
        result, reason = is_phosphoinositide(smi)
        print("SMILES:", smi)
        print("Result:", result)
        print("Reason:", reason)
        print("----------")