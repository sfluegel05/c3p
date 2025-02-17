"""
Classifies: CHEBI:18379 nitrile
"""
"""
Classifies: A compound having the structure RC#N (nitrile)
Definition: A nitrile is a compound in which a carbon atom forms a triple bond with nitrogen,
and that nitrile carbon is substituted with exactly one non‐hydrogen substituent.
Moreover the nitrile nitrogen must be terminal and the bond from the nitrile carbon to the R-group should be a single bond.
If this bond is conjugated then we require that the substituent not be overly embedded in an extended pi-system.
"""

from rdkit import Chem

def is_nitrile(smiles: str):
    """
    Determines if a molecule is a nitrile based on its SMILES string.
    A nitrile is identified by a C≡N functionality.
    For a valid nitrile (RC#N), we require that:
      (a) the nitrile nitrogen is terminal (degree exactly 1),
      (b) the nitrile carbon has exactly one non-hydrogen neighbor aside from the nitrile nitrogen,
      (c) the bond from the nitrile carbon to the substituent is a single bond, and
      (d) if that bond is conjugated then the substituent should not be too highly connected 
          (we allow degree 3 or less, except if the substituent is in a ring).
    Args:
        smiles (str): SMILES string of the target molecule.
    Returns:
        bool: True if at least one valid nitrile (RC#N) group is found that passes these checks.
        str: A reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First, use a simple substructure match for a C triple-bond N.
    nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
    if nitrile_pattern is None:
        return False, "Failed to construct nitrile SMARTS"
    
    matches = mol.GetSubstructMatches(nitrile_pattern)
    if not matches:
        return False, "No C≡N pattern found"
    
    valid_nitrile_found = False
    reasons = []
    
    # Process each match separately.
    for (c_idx, n_idx) in matches:
        atom_c = mol.GetAtomWithIdx(c_idx)
        atom_n = mol.GetAtomWithIdx(n_idx)
        
        # Ensure we identify the nitrile C and nitrile N correctly.
        if atom_c.GetAtomicNum() != 6 or atom_n.GetAtomicNum() != 7:
            continue  # not a valid C≡N
        
        # (a) Check that the nitrile nitrogen is terminal: degree exactly 1.
        if atom_n.GetDegree() != 1:
            continue
        
        # (b) For the nitrile carbon, count its heavy (non-H) neighbors EXCLUDING the nitrile nitrogen.
        r_substituents = []
        for nbr in atom_c.GetNeighbors():
            if nbr.GetIdx() == atom_n.GetIdx():
                continue
            if nbr.GetAtomicNum() > 1:  # non-hydrogen
                r_substituents.append(nbr)
        if len(r_substituents) != 1:
            continue  # We require exactly one R-group attached aside from the nitrile N.
        r_atom = r_substituents[0]
        
        # (c) Confirm that the bond from the nitrile carbon to its substituent is a single bond.
        r_bond = mol.GetBondBetweenAtoms(atom_c.GetIdx(), r_atom.GetIdx())
        if r_bond is None or r_bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
            continue
        
        # (d) If the C–R bond is conjugated, then make sure the substituent is not overly embedded 
        # in an extended π-system. We check: if conjugated and the substituent's heavy-atom degree > 3,
        # and it is not part of a ring, then skip.
        if r_bond.GetIsConjugated() and (r_atom.GetDegree() > 3 and not r_atom.IsInRing()):
            continue
        
        valid_nitrile_found = True
        reasons.append("Found valid C≡N group with properly substituted carbon")
    
    if valid_nitrile_found:
        return True, "; ".join(reasons)
    return False, "No valid nitrile (RC#N) group found"


# Basic testing of the function.
if __name__ == "__main__":
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
        # Some examples that were false positives or negatives in previous version:
        "COCC(=O)N1[C@@H]([C@@H]([C@H]1C#N)C2=CC=C(C=C2)C#CCC3CCCC3)CO",  # false positive example
    ]
    
    for s in test_smiles:
        result, reason = is_nitrile(s)
        print(f"SMILES: {s}\nNitrile: {result}\nReason: {reason}\n")