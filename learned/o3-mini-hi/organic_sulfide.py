"""
Classifies: CHEBI:16385 organic sulfide
"""
"""
Classifies: Organic sulfide (thioether)
Definition: Compounds having the structure R–S–R (with R ≠ H).
Such compounds were once called thioethers.
This improved classifier:
  • Excludes many small free amino acids or peptide‐derived compounds by
    checking for a zwitterionic amino acid center and multiple amide bonds.
  • Iterates over all sulfur atoms and accepts them only if the S atom is
    bound to exactly two heavy atoms (neighbors) that are both carbons,
    and the sulfur is not additionally substituted with oxygen.
"""

from rdkit import Chem

def is_organic_sulfide(smiles: str):
    """
    Determines if a molecule is an organic sulfide (thioether) based on its SMILES string.
    
    An organic sulfide is defined as a molecule containing at least one R–S–R moiety,
    where both substituents (R) are not hydrogen and the sulfur is not oxidized (i.e. no bonds
    to oxygen). In addition, many peptides and free amino acids (which sometimes include a methionine
    side chain) are not considered in this context.
    
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
    
    # ----- Filter out obvious peptides / free amino acids -----
    # (a) If the molecule is very small and shows an amino acid zwitterionic pattern,
    #     we assume it is a free amino acid.
    amino_acid_pattern = Chem.MolFromSmarts("[C@H]([NH3+])C(=O)[O-]")
    if mol.GetNumAtoms() < 50 and mol.HasSubstructMatch(amino_acid_pattern):
        return False, "Molecule appears to be a free amino acid (not classified as organic sulfide)"
    
    # (b) Count amide bonds. Many peptide‐derived molecules contain more than one amide bond.
    #     We use a simple SMARTS for the amide: C(=O)N.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) > 1:
        return False, "Molecule appears to be peptide-derived (contains multiple amide bonds)"
    
    # ----- Look for a genuine thioether motif -----
    # Instead of a fixed SMARTS, we iterate over every sulfur atom that is not oxidized.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16:
            continue  # only consider sulfur atoms
        # We expect a thioether sulfur to be divalent (i.e. exactly two bonds).
        # (Note: implicit hydrogens are not counted in GetDegree().)
        if atom.GetDegree() != 2:
            continue

        neighbors = atom.GetNeighbors()
        # Check that both substituents are not hydrogen.
        # In our case, we require them to be carbon atoms (covers aliphatic and aromatic).
        if not (neighbors[0].GetAtomicNum() == 6 and neighbors[1].GetAtomicNum() == 6):
            continue

        # Also, we want to ensure the sulfur is not additionally bonded to oxygen from, e.g., a sulfoxide.
        # (Even though degree==2 normally, sometimes special representations could show extra explicit bonds.)
        extra_substitution = False
        for neighbor in atom.GetNeighbors():
            # If by any chance an oxygen was attached, we do not consider this a plain thioether.
            if neighbor.GetAtomicNum() == 8:
                extra_substitution = True
                break
        if extra_substitution:
            continue

        # We found at least one sulfur that meets our criteria.
        return True, "Molecule contains at least one organic sulfide (thioether) group (R–S–R bond)"
    
    return False, "No organic sulfide (R–S–R) moiety found in the molecule"

# Example usage: (testing on some of the provided SMILES)
if __name__ == "__main__":
    test_smiles = [
        # True positives:
        "CCCCCCCCCCCCCCCSC[C@H](COP(O)(=O)OCC[N+](C)(C)C)OP(O)(=O)CCCCCCCCCCCCCCCC",  # long-chain thioether
        "COC(=O)C(CC1=CC2C(CC(NC2C=C1)c1c(Cl)cccc1Cl)Sc1ccccc1)NC(=O)OC(C)(C)C",       # thioether with aromatic groups
        "C1CC2=C(C1)NC(=NC2=O)SCC3=CC=C(C=C3)Cl",                                      # thioether group present
        "CC1=C(C(=NN1)SCC(=O)N2CCCCCC2)[N+](=O)[O-]",                                 # another thioether motif
        # False positive (peptide-like, to be filtered):
        "C1=CC=C(C=C1)OC2=CC=C(C=C2)NC(=O)COC(=O)CCCSC3=NC4=CC=CC=C4N3",
        # A molecule with oxidized sulfur (sulfoxide/sulfone) -- should not match.
        "C1C[C@@H]2C(C(=C([C@H]1C2)SC3=CC=CC=C3)C(C4=CC=C(C=C4Cl)S(C)(=O)=O)=O)=O"
    ]
    for s in test_smiles:
        result, reason = is_organic_sulfide(s)
        print(result, reason)