"""
Classifies: CHEBI:18379 nitrile
"""
"""
Classifies: A compound having the structure RC#N (nitrile)
Definition: A nitrile is a compound in which a carbon atom forms a triple bond with nitrogen,
and that nitrile carbon is substituted with exactly one non‐hydrogen substituent.
We further require that the bond from that substituent to the nitrile C is a single bond.
Moreover, if that bond is conjugated (as in acrylonitrile derivatives) then we require that
the substituent is not highly connected (i.e. not part of an extended pi–system) since many
false positives came from such compounds.
"""

from rdkit import Chem

def is_nitrile(smiles: str):
    """
    Determines if a molecule is a nitrile based on its SMILES string.
    A nitrile is characterized by a C≡N functionality where the nitrile nitrogen
    is terminal (degree 1) and the nitrile carbon is bound to exactly one non‐hydrogen substituent.
    Additional checks on the nature of the bond from the nitrile carbon to the R-group (it must
    be a single bond and not extensively conjugated) are used to reduce false positives.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if at least one valid nitrile (RC#N) group is found, and it meets the extra criteria,
              False otherwise.
        str: A reason for the classification.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    valid_nitrile_found = False
    reasons = []
    
    # Iterate over all bonds in the molecule
    for bond in mol.GetBonds():
        # Only consider triple bonds—candidate for nitrile
        if bond.GetBondType() != Chem.rdchem.BondType.TRIPLE:
            continue
        
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        
        # Identify the nitrile C and N atoms
        if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 7:
            nitrile_carbon = atom1
            nitrile_nitrogen = atom2
        elif atom2.GetAtomicNum() == 6 and atom1.GetAtomicNum() == 7:
            nitrile_carbon = atom2
            nitrile_nitrogen = atom1
        else:
            continue  # Not a C≡N bond
        
        # (a) Ensure the nitrile nitrogen is terminal (degree exactly 1)
        if nitrile_nitrogen.GetDegree() != 1:
            continue
        
        # (b) Count the non-H substituents on the nitrile carbon (excluding the nitrile nitrogen)
        r_substituents = []
        for neighbor in nitrile_carbon.GetNeighbors():
            if neighbor.GetIdx() == nitrile_nitrogen.GetIdx():
                continue
            if neighbor.GetAtomicNum() != 1:  # non-hydrogen
                r_substituents.append(neighbor)
        
        if len(r_substituents) != 1:
            continue  # In a simple RC#N the C must have exactly one non-hydrogen substituent (besides the N)
        
        r_atom = r_substituents[0]
        
        # (c) Optionally, check that the bond from nitrile C to the substituent is a single bond.
        r_bond = mol.GetBondBetweenAtoms(nitrile_carbon.GetIdx(), r_atom.GetIdx())
        if r_bond is None or r_bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
            continue
        
        # (d) If the bond from nitrile carbon to the R atom is conjugated, then require that
        # the R atom is not highly connected (i.e. degree <= 2) so that the nitrile group is peripheral.
        if r_bond.GetIsConjugated() and r_atom.GetDegree() > 2:
            continue
        
        # If all conditions are met, mark as valid nitrile.
        valid_nitrile_found = True
        reasons.append("Found valid C≡N group with properly substituted carbon")
    
    if valid_nitrile_found:
        return True, "; ".join(reasons)
    
    return False, "No valid nitrile (RC#N) group found"


# You can do a quick test of this function when run as a script.
if __name__ == "__main__":
    test_smiles = [
        "COc1ccc(CC(C#N)C(\\C#N)=C\\c2ccc(OC)c(OC)c2)cc1",  # NK372135B (true positive)
        "CN(C)CCC[C@]1(OCc2cc(ccc12)C#N)c1ccc(F)cc1",         # escitalopram (true positive)
        "CC(C)(O)C#N",                                       # 2-hydroxy-2-methylpropanenitrile (true positive)
        "O=C(CC#N)c1ccccc1",                                 # benzoylacetonitrile (true positive)
        "ClC(Cl)C#N",                                       # Dichloroacetonitrile (true positive)
        "COc1ccc(cc1)\\C=C(/C#N)c1ccc(OC)cc1",              # 2,3-Bis(p-methoxyphenyl)acrylonitrile (false positive)
    ]
    for s in test_smiles:
        result, reason = is_nitrile(s)
        print(f"SMILES: {s}\nNitrile: {result}\nReason: {reason}\n")