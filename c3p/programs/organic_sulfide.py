"""
Classifies: CHEBI:16385 organic sulfide
"""
"""
Classifies: Organic sulfide (thioether)
Definition: Compounds containing the structure R–S–R (with R ≠ H).
Such compounds were once called thioethers.

In this implementation we:
  • Remove exclusions that target peptide or amide‐rich molecules so as not to miss genuine thioethers.
  • Still filter out obvious small free amino acids (zwitterions) that might otherwise be mis‐classified.
  • Then, for each sulfur atom we require:
       - It is not aromatic.
       - It has exactly two neighbors.
       - Both sulfur–neighbor bonds are SINGLE bonds.
       - Neither neighbor is oxygen (to avoid oxidized S) and both are carbons.
       
If any sulfur meets these criteria the molecule is classified as an organic sulfide.
"""

from rdkit import Chem

def is_organic_sulfide(smiles: str):
    """
    Determines if a molecule is an organic sulfide (thioether) based on its SMILES string.
    
    An organic sulfide is defined as a compound containing at least one R–S–R moiety,
    where both substituents R are non‐hydrogen and the S is not oxidized (i.e. no S–O bonds)
    nor is part of an aromatic (e.g. fused) system.
    
    To avoid flagging small free amino acids, we filter molecules that are small (numAtoms <50)
    and match a typical zwitterionic pattern.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an organic sulfide, False otherwise.
        str: Explanation of the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # ----- Optional filter for very small zwitterionic amino acids -----
    zwitterion = Chem.MolFromSmarts("[C@H]([NH3+])C(=O)[O-]")
    if mol.GetNumAtoms() < 50 and mol.HasSubstructMatch(zwitterion):
        return False, "Molecule appears to be a free amino acid (not classified as organic sulfide)"
    
    # ----- Scan for a genuine thioether motif (R–S–R) -----
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16:  # only consider sulfur atoms
            continue
        
        # Check that the sulfur is not aromatic (avoids S in fused/heterocyclic systems)
        if atom.GetIsAromatic():
            continue
        
        # Require exactly 2 neighbors.
        neighbors = atom.GetNeighbors()
        if len(neighbors) != 2:
            continue
        
        # Both connecting bonds must be single bonds.
        bonds_valid = True
        for nbr in neighbors:
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            if bond is None or bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                bonds_valid = False
                break
        if not bonds_valid:
            continue
        
        # Both neighbors should be carbon atoms (atomic number 6).
        if not (neighbors[0].GetAtomicNum() == 6 and neighbors[1].GetAtomicNum() == 6):
            continue
        
        # Ensure the sulfur is not additionally bonded to any oxygen (avoids oxidized forms).
        if any(nbr.GetAtomicNum() == 8 for nbr in neighbors):
            continue
        
        # Passed all checks: we have found a genuine R–S–R (thioether) moiety.
        return True, "Molecule contains at least one organic sulfide (thioether) group (R–S–R bond)"
    
    return False, "No organic sulfide (R–S–R) moiety found in the molecule"


# Example usage (testing on a few provided SMILES strings):
if __name__ == "__main__":
    test_smiles = [
        # True positives (expected to be organic sulfides)
        "CCCCCCCCCCCCCCCSC[C@H](COP(O)(=O)OCC[N+](C)(C)C)OP(O)(=O)CCCCCCCCCCCCCCCC",  # long-chain thioether
        "COC(=O)C(CC1=CC2C(CC(NC2C=C1)c1c(Cl)cccc1Cl)Sc1ccccc1)NC(=O)OC(C)(C)C",       # thioether with aromatic groups
        "C1CC2=C(C1)NC(=NC2=O)SCC3=CC=C(C=C3)Cl",                                      # thioether group present
        "CC1=C(C(=NN1)SCC(=O)N2CCCCCC2)[N+](=O)[O-]",                                 # another thioether motif
        # A few cases that in previous attempts were flagged as false positives.
        "C1=CC=C(C=C1)OC2=CC=C(C=C2)NC(=O)COC(=O)CCCSC3=NC4=CC=CC=C4N3",               # example false positive from previous run
        "CSCCC(N)C(=O)NO",  # methioninehydroxamic acid (expected NOT to be classified as organic sulfide)
    ]
    for s in test_smiles:
        result, reason = is_organic_sulfide(s)
        print(f"SMILES: {s}\nResult: {result} | Reason: {reason}\n")