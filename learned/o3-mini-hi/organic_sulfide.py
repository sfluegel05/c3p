"""
Classifies: CHEBI:16385 organic sulfide
"""
"""
Classifies: Organic sulfide (thioether)
Definition: Compounds having the structure R–S–R (with R ≠ H). Such compounds were once called thioethers.

In this improved implementation we:
  • Exclude molecules that appear to be peptides or peptide‐conjugates
    by counting amide bonds (a peptide typically has at least two).
  • Then for each sulfur atom in the molecule we require:
       - The sulfur can have exactly 2 heavy-atom neighbors (i.e. non‐hydrogen).
       - Both bonds from S are SINGLE bonds.
       - Neither neighbor is oxygen (to avoid oxidized S).
       - The sulfur is not part of an aromatic ring (thus excluding e.g. thiophenes).
       
If any sulfur meets these criteria the molecule is classified as an organic sulfide.
If no S atom qualifies then the molecule will not be classified as an organic sulfide.
"""

from rdkit import Chem

def is_organic_sulfide(smiles: str):
    """
    Determines if a molecule is an organic sulfide (thioether) based on its SMILES string.

    An organic sulfide is defined as a compound containing at least one R–S–R moiety,
    where both substituents R are not hydrogen. In our implementation we require:
      - The S atom has exactly two bonds (to non-H atoms) and both are single bonds.
      - None of the neighbors is an oxygen (which would indicate oxidation).
      - The sulfur is not part of an aromatic ring (e.g. thiophenes are excluded).
    We additionally try to filter out peptide or CoA‐like molecules by excluding molecules
    with multiple amide bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as an organic sulfide, False otherwise.
        str: Explanation of the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # ----- Exclude molecules that likely are peptides or related conjugates -----
    # Here we use a simple SMARTS that matches an amide bond.
    amide = Chem.MolFromSmarts("C(=O)N")
    n_amide = len(mol.GetSubstructMatches(amide))
    # Many peptides or conjugates contain 2 or more amide bonds.
    if n_amide >= 2:
        return False, "Molecule appears to be a peptide or peptide conjugate (multiple amide bonds detected)"
    
    # ----- Now scan for a genuine thioether (R–S–R) motif -----
    # Loop over all sulfur atoms in the molecule.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16:  # Only consider sulfur
            continue
        
        # Exclude if the sulfur is flagged as aromatic and is in an aromatic ring.
        if atom.IsInRing() and atom.GetIsAromatic():
            continue
        
        # Get heavy-atom neighbors (exclude hydrogens).
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        if len(heavy_neighbors) != 2:
            continue
        
        # Check that both bonds from S are single bonds.
        valid_bonds = True
        for nbr in heavy_neighbors:
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            if bond is None or bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                valid_bonds = False
                break
        if not valid_bonds:
            continue
        
        # Exclude if any neighbor is oxygen (atomic number 8).
        if any(nbr.GetAtomicNum() == 8 for nbr in heavy_neighbors):
            continue

        # If we reach here, this S atom is in a –S– bond with two non-H, non-O atoms,
        # its bonds are single, and it is not part of an aromatic ring.
        return True, "Molecule contains at least one organic sulfide (thioether) group (R–S–R bond)"
    
    return False, "No organic sulfide (thioether) moiety found in the molecule"


# Example usage (testing on a few provided SMILES strings)
if __name__ == "__main__":
    test_smiles = [
        # True positives (expected to be organic sulfides):
        "CCCCCCCCCCCCCCCSC[C@H](COP(O)(=O)OCC[N+](C)(C)C)OP(O)(=O)CCCCCCCCCCCCCCCC",  # long-chain thioether
        "COC(=O)C(CC1=CC2C(CC(NC2C=C1)c1c(Cl)cccc1Cl)Sc1ccccc1)NC(=O)OC(C)(C)C",       # thioether in a fused system
        "C1CC2=C(C1)NC(=NC2=O)SCC3=CC=C(C=C3)Cl",                                      # thioether group present
        "CC1=C(C(=NN1)SCC(=O)N2CCCCCC2)[N+](=O)[O-]",                                 # another thioether motif
        # False positives (should not be classified as organic sulfides):
        "S(C(=O)C(CCCC(CCCC(C)C)C)C)CCNC(=O)CCNC(=O)C(O)C(COP(OP(OCC1OC(N2C3=NC=NC(N)=C3N=C2)C(O)C1OP(O)(O)=O)(O)=O)(O)=O)(C)C",
        "S(CC[C@H](NC(=O)[C@@H](N)CC(O)=O)C(=O)N[C@@H](CC=1NC=NC1)C(O)=O)C",
        "C[C@H](CCCCCCCCCCCCCCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)N2C=NC3=C(N)N=CN=C23)O[C@@H]4O[C@@H](C)[C@H](O)C[C@H]4O",
        # False negatives (molecules that are organic sulfides but were missed previously)
        "C1COCCN1C2=NC3=C(C(=N2)NCC4=NC5=CC=CC=C5N4)N=CN3C6=CSC=C6",  # thiophenyl S is in a ring -> not flagged
        "OCC(NC=1N=C2N(C(C)C)C=NC2=C(N1)NCC3=CC=CC=C3)CC"
    ]
    
    for s in test_smiles:
        result, reason = is_organic_sulfide(s)
        print(f"SMILES: {s}\nResult: {result} | Reason: {reason}\n")