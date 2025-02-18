"""
Classifies: CHEBI:16385 organic sulfide
"""
"""
Classifies: Organic sulfide (thioether)
Definition: Compounds having the structure R–S–R (with R ≠ H).
Such compounds were once called thioethers.

Improvements compared with the previous classifier:
  • Excludes many peptides or free amino acids by checking for a typical amino acid zwitterion
    and by counting amide bonds.
  • For each sulfur atom, requires:
       - Exactly two explicit neighbors.
       - Each S–neighbor bond is a single bond.
       - Both neighbors are carbon (atomic number 6).
       - The sulfur is not flagged as aromatic (so that fused/sulfur‐in‐ring systems are avoided).
  • This more stringent approach helps to lower false positives.
"""

from rdkit import Chem

def is_organic_sulfide(smiles: str):
    """
    Determines if a molecule is an organic sulfide (thioether) based on its SMILES string.
    
    An organic sulfide is defined as a molecule containing at least one R–S–R moiety,
    where both substituents (R) are not hydrogen and the sulfur is not oxidized (e.g. no bonds
    to oxygen) and is not part of a fused aromatic moiety. In addition, we filter out
    many simple peptides or free amino acids.
    
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
    # (a) Check for a typical amino acid zwitterion (e.g. methionine side-chain in a free amino acid)
    amino_acid_pattern = Chem.MolFromSmarts("[C@H]([NH3+])C(=O)[O-]")
    if mol.GetNumAtoms() < 50 and mol.HasSubstructMatch(amino_acid_pattern):
        return False, "Molecule appears to be a free amino acid (not classified as organic sulfide)"
    
    # (b) Count amide bonds. Many peptide‐derived molecules contain more than one amide bond.
    # There is a risk that a lone amide is present in a non‐peptidic context
    # so we only disqualify if there is more than one.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) > 1:
        return False, "Molecule appears to be peptide-derived (contains multiple amide bonds)"
    
    # ----- Scan for a genuine thioether motif (R–S–R) -----
    # Iterate over every sulfur atom. For each sulfur, we demand:
    #   1. The S atom is not aromatic.
    #   2. It has exactly two neighbors (explicit heavy atoms).
    #   3. Both connecting bonds are SINGLE bonds.
    #   4. Both neighbors are carbon atoms.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16:  # Only consider sulfur atoms.
            continue
        
        # Check if the sulfur atom is aromatic.
        if atom.GetIsAromatic():
            continue  # likely part of a fused or heterocyclic system (not our desired thioether)
        
        # We require exactly 2 explicit neighbors.
        neighbors = atom.GetNeighbors()
        if len(neighbors) != 2:
            continue
        
        # Check each bond from S to neighbor is SINGLE.
        valid_bonds = True
        for nbr in neighbors:
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            if bond is None or bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                valid_bonds = False
                break
        if not valid_bonds:
            continue
        
        # Check that both neighbors are carbon atoms.
        if not (neighbors[0].GetAtomicNum() == 6 and neighbors[1].GetAtomicNum() == 6):
            continue

        # If we pass all the checks, then we have found a genuine thioether moiety.
        return True, "Molecule contains at least one organic sulfide (thioether) group (R–S–R bond)"
    
    return False, "No organic sulfide (R–S–R) moiety found in the molecule"

# Example usage: (testing on a subset of the provided SMILES)
if __name__ == "__main__":
    test_smiles = [
        # True positives (should classify as organic sulfide)
        "CCCCCCCCCCCCCCCSC[C@H](COP(O)(=O)OCC[N+](C)(C)C)OP(O)(=O)CCCCCCCCCCCCCCCC",  # long-chain thioether
        "COC(=O)C(CC1=CC2C(CC(NC2C=C1)c1c(Cl)cccc1Cl)Sc1ccccc1)NC(=O)OC(C)(C)C",       # thioether with aromatic groups
        "C1CC2=C(C1)NC(=NC2=O)SCC3=CC=C(C=C3)Cl",                                      # thioether group present
        "CC1=C(C(=NN1)SCC(=O)N2CCCCCC2)[N+](=O)[O-]",                                 # another thioether motif
        # False positive candidates (expected to be rejected due to peptide patterns or non-thioether S)
        "CN1C(C(=O)Nc2ccccn2)=C(O)c2sccc2S1(=O)=O",  # tenoxicam-like structure, S in heterocycle
        "CCOC(=O)C1C(C2=C(CCCC2=O)N=C1C)C3=CC=CS3"     # thiophen-containing, not a genuine R–S–R pattern
    ]
    for s in test_smiles:
        result, reason = is_organic_sulfide(s)
        print(f"SMILES: {s}\nResult: {result} | Reason: {reason}\n")