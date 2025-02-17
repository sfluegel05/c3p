"""
Classifies: CHEBI:18379 nitrile
"""
"""
Classifies: A compound having the structure RC#N (nitrile)
Definition: A nitrile is a compound in which a carbon atom forms a triple bond with nitrogen,
where the carbon is substituted with exactly one non‐hydrogen substituent (i.e. a C‐substituted derivative of hydrocyanic acid).
"""

from rdkit import Chem

def is_nitrile(smiles: str):
    """
    Determines if a molecule is a nitrile based on its SMILES string.
    A nitrile is characterized by a carbon-nitrogen triple bond (C≡N) where the nitrile carbon 
    is bound to exactly one substituent that is not hydrogen.

    This improved version enforces that:
      - The triple bond is between a carbon (atomic number 6) and nitrogen (atomic number 7).
      - The nitrile nitrogen has degree exactly 1.
      - The nitrile carbon has exactly one neighbor (other than the nitrile nitrogen) that is not a hydrogen.
        This helps exclude cases where the C atom is over‐substituted.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if at least one valid nitrile (RC#N) group is found, False otherwise.
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
        # Check if the bond is a triple bond
        if bond.GetBondType() != Chem.rdchem.BondType.TRIPLE:
            continue

        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()

        # Determine which is carbon and which is nitrogen
        if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 7:
            nitrile_carbon = atom1
            nitrile_nitrogen = atom2
        elif atom2.GetAtomicNum() == 6 and atom1.GetAtomicNum() == 7:
            nitrile_carbon = atom2
            nitrile_nitrogen = atom1
        else:
            continue  # not a C≡N triple bond

        # Check that the nitrile nitrogen is terminal: degree should be exactly 1
        if nitrile_nitrogen.GetDegree() != 1:
            continue

        # Count non-hydrogen substituents on the nitrile carbon excluding the nitrile nitrogen.
        substituent_count = 0
        for neighbor in nitrile_carbon.GetNeighbors():
            if neighbor.GetIdx() == nitrile_nitrogen.GetIdx():
                continue
            # Only count non-hydrogen atoms
            if neighbor.GetAtomicNum() != 1:
                substituent_count += 1

        # For a true nitrile group RC#N, the nitrile carbon should have exactly one non-hydrogen substituent.
        if substituent_count == 1:
            valid_nitrile_found = True
            reasons.append("Found C≡N group with properly substituted carbon")
        # If there are more than one substituents the nitrile might be embedded in a more complex group
        # and should not be considered a simple nitrile.
    
    if valid_nitrile_found:
        return True, "; ".join(reasons)
    
    return False, "No valid nitrile (RC#N) group found"

# You can test the function with a few SMILES examples:
if __name__ == "__main__":
    test_smiles = [
        "COc1ccc(CC(C#N)C(\\C#N)=C\\c2ccc(OC)c(OC)c2)cc1",  # NK372135B
        "CC(C)(O)C#N",  # 2-hydroxy-2-methylpropanenitrile; expected True
        "O=C(CC#N)c1ccccc1",  # benzoylacetonitrile; expected True
        "ClC(Cl)C#N",  # Dichloroacetonitrile; expected True
    ]
    for s in test_smiles:
        result, reason = is_nitrile(s)
        print(f"SMILES: {s}\nNitrile: {result}\nReason: {reason}\n")