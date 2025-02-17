"""
Classifies: CHEBI:18179 phosphoinositide
"""
#!/usr/bin/env python
"""
Classifies: Phosphoinositide
Definition: Any phosphatidylinositol that is phosphorylated at one or more of the hydroxy groups of inositol.
A phosphoinositide must contain a glycerophosphatidylinositol “core” with at least two acyl chains (detected via ester groups)
and a phosphate bridging the glycerol backbone to an inositol headgroup – a cyclohexane (six-membered carbon ring) decorated 
with oxygen substituents (–OH or –O–P(=O)(O) groups, where the phosphate is neutral).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string, using heuristic criteria.
    
    Heuristic criteria:
      1. The SMILES must be valid.
      2. The molecule must contain at least two phosphorus atoms.
      3. The molecule must contain at least two acyl ester groups (the C(=O)O motif) to indicate two acyl chains.
      4. The molecule must possess an inositol headgroup that is phosphorylated—that is,
         it contains a six-membered ring of carbons (typical of inositols) that is decorated
         with oxygen substituents. For our purposes the ring must have at least three such substituents,
         and at least one of them should be a phosphate substituent (i.e. an oxygen attached to a neutral P).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets our criteria for a phosphoinositide, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Add explicit hydrogens for better detection of -OH groups.
    mol = Chem.AddHs(mol)
    
    # Criterion 1: Check for at least two phosphorus atoms (atomic number 15).
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count < 2:
        return False, f"Only {p_count} phosphorus atom(s) found; need at least 2 for a phosphoinositide."
    
    # Criterion 2: Look for at least two acyl ester groups.
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found only {len(ester_matches)} acyl ester group(s); need at least 2 acyl chains."
    
    # Criterion 3: Identify a phosphorylated inositol head group.
    # We search for a six-membered ring whose atoms are all carbons (i.e. a cyclohexane),
    # and then check the ring atoms’ external substituents.
    # We count an external substituent if it is an oxygen that is either part of an -OH
    # (i.e. attached to at least one hydrogen) or it is attached to a neutral phosphorus atom.
    found_inositol = False
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        # Look only at 6-membered rings.
        if len(ring) != 6:
            continue
        # Check that every atom in the ring is carbon.
        if not all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
            continue
        
        # Now, for each ring atom, look at substituents that are not members of the ring.
        total_substituents = 0
        phosphate_substituents = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                # Skip neighbor if it is in the ring.
                if nbr.GetIdx() in ring:
                    continue
                # We are interested in oxygen substituents.
                if nbr.GetAtomicNum() == 8:
                    # Check if this oxygen is part of a hydroxyl group (i.e., has at least one H neighbor).
                    has_h = any(n.GetAtomicNum() == 1 for n in nbr.GetNeighbors())
                    # Alternatively, check if the oxygen is attached to a neutral phosphorus.
                    is_phosphate = False
                    for second_nbr in nbr.GetNeighbors():
                        if second_nbr.GetAtomicNum() == 15 and second_nbr.GetFormalCharge() == 0:
                            is_phosphate = True
                            break
                    if has_h or is_phosphate:
                        total_substituents += 1
                        if is_phosphate:
                            phosphate_substituents += 1
        # To count as an inositol head, require at least three such substituents,
        # with at least one being a phosphate substituent.
        if total_substituents >= 3 and phosphate_substituents >= 1:
            found_inositol = True
            break
    
    if not found_inositol:
        return False, ("No phosphorylated inositol head detected: "
                       "a six-membered carbon ring with at least three oxygen substituents (including at least one phosphate) was not found.")
    
    return True, ("Molecule contains at least two phosphorus atoms, at least two acyl ester groups, and a phosphorylated inositol head "
                  "as indicated by detection of a cyclohexane ring decorated with -OH and -O-P(=O)(O) substituents.")


# Debug/test code (executes only when run as a script)
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