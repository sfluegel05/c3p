"""
Classifies: CHEBI:46761 dipeptide
"""
"""
Classifies: dipeptide
Definition: Any molecule that contains two amino‐acid residues connected by peptide linkages.
This function looks for a peptide (amide) bond and then checks (heuristically) for two alpha–carbons.
It also accepts cyclic dipeptides (diketopiperazines) where two amide bonds are present in a 6‐membered ring.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_dipeptide(smiles: str):
    """
    Determines if a given SMILES string corresponds to a dipeptide.
    For linear dipeptides, the molecule should contain exactly one peptide (amide) bond
    connecting two amino-acid residues; for cyclic dipeptides a 6-membered ring (diketopiperazine)
    with 2 peptide bonds is expected.
    
    Additionally, for linear dipeptides, we search (heuristically) for 2 alpha–carbons,
    defined as carbons that are connected to a nitrogen and to a carbonyl group (C=O).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if recognized as a dipeptide, False otherwise
        str: Reason for classification or failure message
    """
    # Parse SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Define SMARTS pattern for a peptide (amide) bond: C(=O)N
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    peptide_matches = mol.GetSubstructMatches(amide_pattern)
    
    # Check for cyclic dipeptide (diketopiperazine)
    # A typical diketopiperazine has a 6-membered ring that contains 2 amide bonds.
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) == 6:
            amide_in_ring = 0
            for match in peptide_matches:
                # match[0] is the carbonyl carbon and match[1] the peptide nitrogen.
                if match[0] in ring and match[1] in ring:
                    amide_in_ring += 1
            if amide_in_ring == 2:
                return True, "Cyclic dipeptide (diketopiperazine) detected with 2 peptide bonds in a 6-membered ring."
    
    # For a linear dipeptide, we expect exactly one peptide (amide) bond.
    if len(peptide_matches) != 1:
        return False, f"Expected 1 peptide bond for a linear dipeptide, found {len(peptide_matches)}."

    # Heuristically, count the number of alpha–carbons.
    # An alpha–carbon (of an amino acid residue) is typically connected to:
    #   (i) an amino group (N) and 
    #  (ii) a carbon that is part of a carbonyl group (C double-bonded to O).
    alpha_carbon_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:  # only consider carbon atoms
            continue
        neighbors = atom.GetNeighbors()
        has_nitrogen = any(nb.GetAtomicNum() == 7 for nb in neighbors)
        has_carbonyl = False
        # Look through neighbors that are carbons to see if any is a carbonyl carbon.
        for nb in neighbors:
            if nb.GetAtomicNum() == 6:
                for bond in nb.GetBonds():
                    # Check if the bond is a double bond to oxygen.
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        other = bond.GetOtherAtom(nb)
                        if other.GetAtomicNum() == 8:
                            has_carbonyl = True
                            break
                if has_carbonyl:
                    break
        if has_nitrogen and has_carbonyl:
            alpha_carbon_count += 1

    if alpha_carbon_count != 2:
        return False, f"Expected 2 alpha–carbons (one per amino acid residue), found {alpha_carbon_count}."

    return True, "Dipeptide detected: one peptide bond (linear) with 2 alpha–carbons (or cyclic diketopiperazine pattern)."


# Example usage (you can remove or comment these lines when using this function in production):
if __name__ == "__main__":
    # Some test examples (Ala-Val and a cyclic dipeptide are provided)
    test_smiles = [
        # Linear dipeptide example: Ala-Val
        "CC(C)[C@H](NC(=O)[C@H](C)N)C(O)=O",
        # Cyclic dipeptide (diketopiperazine) example:
        "C1C[C@@]2([H])C(N[C@H](C(N2C1)=O)CC(C)C)=O"
    ]
    for sm in test_smiles:
        result, reason = is_dipeptide(sm)
        print(f"SMILES: {sm}\nResult: {result}\nReason: {reason}\n")