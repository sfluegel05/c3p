"""
Classifies: CHEBI:18379 nitrile
"""
"""
Classifies: A compound having the structure RC#N (nitrile)
Definition: A nitrile is a compound in which a carbon atom forms a triple bond with nitrogen,
and that nitrile carbon is substituted with exactly one non‐hydrogen substituent.
Moreover, the nitrile nitrogen must be terminal (degree=1) and the bond from the nitrile carbon
to the R–group should be a single bond.
If that C–R bond is conjugated then we further require that the R–atom is not overly embedded 
in an extended π‐system (i.e. if its heavy-atom degree > 3 and it is not in a ring, we do not accept it).
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_nitrile(smiles: str):
    """
    Determines if a molecule is a nitrile based on its SMILES string.
    For a valid nitrile (RC#N), the following conditions must be met:
      (a) There exists at least one triple bond connecting a carbon and a nitrogen.
      (b) The nitrile nitrogen is terminal (its degree is exactly 1) and is sp hybridized and non‐aromatic.
      (c) The nitrile carbon (the carbon atom of the C≡N bond) has exactly two neighbors:
          the nitrile nitrogen and one R–group (a non-hydrogen substituent). Also the nitrile carbon must be sp hybridized and non‐aromatic.
      (d) The bond from the nitrile carbon to its R–group is a single bond.
      (e) If that C–R bond is conjugated then we require that the R–atom is not “overly embedded”
          in an extended π-system (we flag it if its heavy-atom degree > 3 and it is not in a ring).
    
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
    
    # Iterate over all bonds to find triple bonds.
    for bond in mol.GetBonds():
        if bond.GetBondType() != rdchem.BondType.TRIPLE:
            continue
        
        # Get the two atoms in the triple bond.
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        
        # We require one atom to be carbon (#6) and the other nitrogen (#7)
        if a1.GetAtomicNum() == 7 and a2.GetAtomicNum() == 6:
            nitrileN = a1
            nitrileC = a2
        elif a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 7:
            nitrileC = a1
            nitrileN = a2
        else:
            continue  # This triple bond is not between C and N.
        
        # Check that the nitrile nitrogen is terminal (degree==1), sp hybridized, and not aromatic.
        if nitrileN.GetDegree() != 1:
            continue
        if nitrileN.GetHybridization() != rdchem.HybridizationType.SP:
            continue
        if nitrileN.GetIsAromatic():
            continue

        # Check that the nitrile carbon has exactly 2 neighbors (nitrile N and one R–group),
        # is sp hybridized, and not aromatic.
        neighbors = nitrileC.GetNeighbors()
        if len(neighbors) != 2:
            continue
        if nitrileC.GetHybridization() != rdchem.HybridizationType.SP:
            continue
        if nitrileC.GetIsAromatic():
            continue

        # Identify the R–group atom (the neighbor that is not the nitrile nitrogen).
        r_atoms = [nbr for nbr in neighbors if nbr.GetIdx() != nitrileN.GetIdx() and nbr.GetAtomicNum() > 1]
        if len(r_atoms) != 1:
            continue
        r_atom = r_atoms[0]
        
        # Ensure the bond from the nitrile carbon to the R group is a SINGLE bond.
        r_bond = mol.GetBondBetweenAtoms(nitrileC.GetIdx(), r_atom.GetIdx())
        if r_bond is None or r_bond.GetBondType() != rdchem.BondType.SINGLE:
            continue

        # Check extra condition: if the C–R bond is conjugated, flag if the R-atom is overly embedded.
        # Here we require that if the bond is conjugated and the R-atom has heavy (non-H) degree > 3 and is not in a ring then we skip.
        if r_bond.GetIsConjugated() and (r_atom.GetDegree() > 3 and not r_atom.IsInRing()):
            continue

        # If we arrive here, we have a valid nitrile group.
        valid_nitrile_found = True
        reasons.append("Found valid C≡N group with properly substituted carbon (nitrile N is terminal and both atoms are sp-hybridized)")
    
    if valid_nitrile_found:
        return True, "; ".join(reasons)
    return False, "No valid nitrile (RC≡N) group found"


# Example usage and testing (the provided examples and one known false positive from previous attempt)
if __name__ == "__main__":
    test_smiles = [
        "COc1ccc(CC(C#N)C(\\C#N)=C\\c2ccc(OC)c(OC)c2)cc1",   # NK372135B (TP)
        "CN(C)CCC[C@]1(OCc2cc(ccc12)C#N)c1ccc(F)cc1",         # escitalopram (TP)
        "CC(C)(O)C#N",                                       # 2-hydroxy-2-methylpropanenitrile (TP)
        "C(#N)[C@H](O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H](O1)CO)",  # (R)-amygdalin (TP)
        "[C@]1(O[C@@H]([C@H]([C@H]1O)O)CO)(C#N)C=2N3N=CN=C(C3=CC2)N", # GS-441524 (TP)
        "C1[C@@H]2C([C@@H](N2)CN1CC3=COC=N3)C4=CC=C(C=C4)C5=CC=CC(=C5)C#N", # 3-[4-[(1S,5R)-3-(4-oxazolylmethyl)...]benzonitrile (TP)
        "CCC(=O)N1C[C@H]2[C@@H]([C@H](N2C(=O)C1)CO)C3=CC=C(C=C3)C4=CC=C(C=C4)C#N",  # 4-[4-[(6R,...]benzonitrile (TP)
        "C1CCN2[C@H](CN(C1)CC3=CC=CC=C3F)[C@@H]([C@@H]2CO)C4=CC=C(C=C4)C5=CC=CC(=C5)C#N", # 3-[4-[(8S,9S,10R)...]benzonitrile (TP)
        "CCOP(=O)(OCC)SCC(=O)NC(C)(C)C#N",                   # Cyanthoate (TP)
        "Nc1c(SC(F)(F)F)c(nn1-c1c(Cl)cc(cc1Cl)C(F)(F)F)C#N", # fipronil-sulfide (TP)
        "CC1=CN(CCN2CCN(CC2)C2=CC=CC=C2C2=CC(=CC(=C2)C#N)C(=O)NCCCN2CCCC2)C2=CC=CC=C12", # LLY-507 (TP)
        "C[C@@H]1CCN(C[C@@H]1N(C)c1ncnc2[nH]ccc12)C(=O)CC#N", # tofacitinib (TP)
        "C1=CC=CC=C1CC[C@@H](C#N)O",                          # (2S)-2-hydroxy-4-phenylbutanenitrile (TP)
        "Clc1ccc(CC[C@](Cn2cncn2)(C#N)c2ccccc2)cc1",         # (S)-fenbuconazole (TP)
        "O=C(CC#N)c1ccccc1",                                 # benzoylacetonitrile (TP)
        "ClC(Cl)C#N",                                       # Dichloroacetonitrile (TP)
        # One known false positive from previous attempt:
        "COCC(=O)N1[C@@H]([C@@H]([C@H]1C#N)C2=CC=C(C=C2)C#CCC3CCCC3)CO",
    ]
    
    for s in test_smiles:
        result, reason = is_nitrile(s)
        print(f"SMILES: {s}\nNitrile: {result}\nReason: {reason}\n")