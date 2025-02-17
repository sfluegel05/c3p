"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
#!/usr/bin/env python
"""
Classifies: N-acetyl-amino acid
Definition: An N-acyl-amino acid that has acetyl (CH3C(=O)-) as the acyl group.
Improved Heuristic:
  1. The molecule must contain an N-acetyl amino acid backbone. This is recognized as a fragment where
     a nitrogen is bound to an acetyl group (CH3C(=O)-), and then to an alpha carbon that in turn is connected
     to a carboxyl group. For amino acids with a chiral alpha carbon, we expect a SMARTS match with [C@H] or [C@@H];
     for achiral glycine the pattern uses [CH2].
  2. We require that the entire molecule contains exactly one such backbone fragment.
"""

from rdkit import Chem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N‐acetyl‐amino acid based on its SMILES string.
    
    This improved heuristic checks for the presence of an N-acetyl amino acid backbone.
    The expected fragment is:
      - An acetyl group (CH3C(=O)-) attached to a nitrogen.
      - The nitrogen is bound to an alpha carbon that carries a carboxyl group (C(=O)[O;H1,O-]).
      - For most amino acids the alpha carbon is chiral ([C@H] or [C@@H]), but glycine is achiral and is represented as CH2.
    
    The function searches for any of three SMARTS patterns:
       1. "CC(=O)N[C@H]([*])C(=O)[O;H1,O-]"
       2. "CC(=O)N[C@@H]([*])C(=O)[O;H1,O-]"
       3. "CC(=O)N[CH2]C(=O)[O;H1,O-]"  (for glycine)
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an N-acetyl-amino acid, False otherwise.
        str: A reason explaining the decision.
    """
    # Parse the input molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for the N-acetyl amino acid backbone.
    # Pattern for chiral alpha-carbon (R-configured)
    smarts_chiral1 = "CC(=O)N[C@H]([*])C(=O)[O;H1,O-]"
    # Pattern for chiral alpha-carbon (S-configured)
    smarts_chiral2 = "CC(=O)N[C@@H]([*])C(=O)[O;H1,O-]"
    # Pattern for glycine (achiral alpha carbon CH2)
    smarts_glycine = "CC(=O)N[CH2]C(=O)[O;H1,O-]"
    
    try:
        frag_chiral1 = Chem.MolFromSmarts(smarts_chiral1)
        frag_chiral2 = Chem.MolFromSmarts(smarts_chiral2)
        frag_glycine = Chem.MolFromSmarts(smarts_glycine)
    except Exception as e:
        return False, "Error in constructing SMARTS patterns: " + str(e)
    
    # Get substructure matches for each fragment
    matches = set()
    for frag in (frag_chiral1, frag_chiral2, frag_glycine):
        if frag is None:
            return False, "Error in constructing one of the SMARTS fragments"
        # Using GetSubstructMatches returns a tuple of tuples (atom indices for each match)
        frag_matches = mol.GetSubstructMatches(frag, useChirality=True)
        for match in frag_matches:
            matches.add(match)  # add each unique match (order of atoms does not matter)
    
    n_matches = len(matches)
    
    # We expect exactly one occurrence of the N-acetyl amino acid backbone fragment.
    if n_matches == 0:
        return False, "No N-acetyl amino acid backbone fragment found"
    if n_matches > 1:
        return False, f"Multiple N-acetyl amino acid backbones detected ({n_matches}); likely a peptide chain or extra acyl groups present"
    
    return True, "Contains a single N-acetyl amino acid backbone (i.e., an acetyl group bound to the nitrogen of an amino acid with a carboxyl group)"

# Optional: Testing the function with several examples.
if __name__ == "__main__":
    test_examples = [
        # True positives:
        "CC(=O)N[C@@H](CCCNC(N)=N)C(O)=O",      # N(alpha)-acetyl-L-arginine
        "C[C@H](NC(C)=O)C(O)=O",                # N-acetyl-L-alanine (chiral) 
        "CC(=O)N1CCC[C@H]1C(O)=O",               # N-acetyl-L-proline
        "CC(=O)N[C@@H](Cc1ccccc1)C(O)=O",        # N-acetyl-L-phenylalanine
        "CC(=O)NCCC[C@H](N)C(O)=O",              # N(5)-acetyl-L-ornithine
        "CC(=O)NCC(O)=O",                       # N-acetylglycine (achiral alpha carbon)
        # Examples that should NOT be classified (false positives/negatives):
        "[H][C@]1(O[C@](O)(C[C@H](O)[C@H]1NC(C)=O)C(O)=O)[C@H](O)[C@H](O)CO",  # N-acetyl-alpha-neuraminic acid
        "OC(=O)[C@H](NC(=O)C)CCC(=O)N",         # N(2)-acetyl-D-glutamine (extra amide bond in sidechain)
    ]
    
    for sm in test_examples:
        result, reason = is_N_acetyl_amino_acid(sm)
        print(f"SMILES: {sm}\n  Classified? {result}\n  Reason: {reason}\n")