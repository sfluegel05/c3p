"""
Classifies: CHEBI:18379 nitrile
"""
"""
Classifies: A compound having the structure RC#N (nitrile)
Definition: A nitrile is a compound in which a carbon atom forms a triple bond with nitrogen,
and that nitrile carbon is substituted with exactly one non‐hydrogen substituent.
Moreover, the nitrile nitrogen must be terminal (degree=1) and the bond from the nitrile carbon
to the R–group should be a single bond.
If that C–R bond is conjugated then we further require that the R–atom is not overly embedded in an extended π‐system.
"""

from rdkit import Chem

def is_nitrile(smiles: str):
    """
    Determines if a molecule is a nitrile based on its SMILES string.
    For a valid nitrile (RC#N), we require:
      (a) There is at least one triple bond connecting a carbon and a nitrogen.
      (b) The nitrile nitrogen is terminal (i.e. has a degree exactly equal to 1).
      (c) The nitrile carbon (the carbon of the C≡N bond) has exactly two neighbors:
          the nitrile nitrogen and one R–group (non-hydrogen).
      (d) The bond from the nitrile carbon to its R–group is a single bond.
      (e) If that C–R bond is conjugated then we require that the substituent atom (R)
          is not “overly embedded” in an extended π-system (we flag it if its heavy-atom degree is > 3 and it is not in a ring).
    Args:
        smiles (str): SMILES string for the input molecule.
    Returns:
        bool: True if at least one valid nitrile (RC#N) group is found that passes the checks.
        str: A reason message for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    valid_nitrile_found = False
    reasons = []

    # Iterate over all bonds and look explicitly for triple bonds.
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.rdchem.BondType.TRIPLE:
            continue
        
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        
        # Identify which atom is nitrogen and which is carbon, if possible.
        # (We expect one to be carbon (atomic #6) and the other to be nitrogen (#7).)
        if a1.GetAtomicNum() == 7 and a2.GetAtomicNum() == 6:
            nitrileN = a1
            nitrileC = a2
        elif a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 7:
            nitrileC = a1
            nitrileN = a2
        else:
            continue  # not a nitrile bond
        
        # (a) Check that the nitrile nitrogen is terminal - degree should be exactly 1.
        if nitrileN.GetDegree() != 1:
            continue

        # (b) Check that the nitrile carbon is substituted only with the nitrile N and one other heavy (non-H) atom.
        neighbors = nitrileC.GetNeighbors()
        if len(neighbors) != 2:
            continue  # It must only have two neighbors: the nitrile N and one R-group.
        
        # Identify the R-group atom: the neighbor that is not the nitrile nitrogen.
        r_neighbors = [nbr for nbr in neighbors if nbr.GetIdx() != nitrileN.GetIdx() and nbr.GetAtomicNum() > 1]
        if len(r_neighbors) != 1:
            continue
        r_atom = r_neighbors[0]
        
        # (c) The bond from the nitrile carbon to the R–group should be a single bond.
        r_bond = mol.GetBondBetweenAtoms(nitrileC.GetIdx(), r_atom.GetIdx())
        if r_bond is None or r_bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
            continue
        
        # (d) If the C–R bond is conjugated we require that the substituent is not 
        # excessively embedded: if conjugated, then flag if the R–atom has heavy-atom degree > 3 and is not in a ring.
        if r_bond.GetIsConjugated() and (r_atom.GetDegree() > 3 and not r_atom.IsInRing()):
            continue
        
        # Passed all tests for this nitrile group.
        valid_nitrile_found = True
        reasons.append("Found valid C≡N group with properly substituted carbon")
    
    if valid_nitrile_found:
        return True, "; ".join(reasons)
    return False, "No valid nitrile (RC#N) group found"

# Example usage and testing
if __name__ == "__main__":
    # A few test SMILES (including reported true positives and a known false-positive)
    test_smiles = [
        "COc1ccc(CC(C#N)C(\\C#N)=C\\c2ccc(OC)c(OC)c2)cc1",   # NK372135B (TP)
        "CN(C)CCC[C@]1(OCc2cc(ccc12)C#N)c1ccc(F)cc1",         # escitalopram (TP)
        "CC(C)(O)C#N",                                       # 2-hydroxy-2-methylpropanenitrile (TP)
        "C(#N)[C@H](O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H](O1)CO)",  # (R)-amygdalin (TP)
        "[C@]1(O[C@@H]([C@H]([C@H]1O)O)CO)(C#N)C=2N3N=CN=C(C3=CC2)N", # GS-441524 (TP)
        "CCOP(=O)(OCC)SCC(=O)NC(C)(C)C#N",                   # Cyanthoate (TP)
        "C[C@@H]1CCN(C[C@@H]1N(C)c1ncnc2[nH]ccc12)C(=O)CC#N", # tofacitinib (TP)
        "C1=CC=CC=C1CC[C@@H](C#N)O",                          # (2S)-2-hydroxy-4-phenylbutanenitrile (TP)
        "Clc1ccc(CC[C@](Cn2cncn2)(C#N)c2ccccc2)cc1",         # (S)-fenbuconazole (TP)
        "O=C(CC#N)c1ccccc1",                                 # benzoylacetonitrile (TP)
        "ClC(Cl)C#N",                                       # Dichloroacetonitrile (TP)
        # Example of a known false positive previously:
        "COCC(=O)N1[C@@H]([C@@H]([C@H]1C#N)C2=CC=C(C=C2)C#CCC3CCCC3)CO",
    ]
    
    for s in test_smiles:
        result, reason = is_nitrile(s)
        print(f"SMILES: {s}\nNitrile: {result}\nReason: {reason}\n")