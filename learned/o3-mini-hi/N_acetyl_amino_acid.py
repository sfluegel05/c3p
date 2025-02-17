"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
#!/usr/bin/env python
"""
Classifies: CHEBI:??? N-acetyl-amino acid
Definition: An N-acyl-amino acid that has acetyl (CH3C(=O)-) as the acyl group.
Heuristic:
  1. Exactly one occurrence of an acetyl group (detected by either "[CH3]C(=O)N" or "N[C](=O)[CH3]").
  2. Exactly one free (terminal) carboxyl group (pattern "C(=O)[O;H1,O-]").
  3. Exactly one amide bond overall (pattern "[C](=O)N"); extra amide bonds indicate a peptide.
"""

from rdkit import Chem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl amino acid based on its SMILES string.

    The following conditions are enforced:
      1. The molecule must have exactly one acetyl fragment. This is checked by detecting the acetyl moiety either 
         as a substructure matching "[CH3]C(=O)N" or "N[C](=O)[CH3]". (Either order is allowed.)
      2. There must be exactly one free carboxyl group, represented by the SMARTS "C(=O)[O;H1,O-]". 
         (The carboxyl group of an amino acid is free; amide carbonyls in peptide bonds are excluded.)
      3. The molecule must have exactly one amide bond overall (SMARTS "[C](=O)N"). 
         (Additional amide bonds suggest the molecule is a di- or oligo-peptide.)
         
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if classified as a valid N-acetyl amino acid; False otherwise.
        str: A reason explaining the decision.
    """
    # Parse the SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for the acetyl fragment.
    # Two SMARTS are used so that the order of atoms in the SMILES does not lead to a mis-match.
    acetyl_smarts1 = "[CH3]C(=O)N"
    acetyl_smarts2 = "N[C](=O)[CH3]"
    try:
        acetyl_frag1 = Chem.MolFromSmarts(acetyl_smarts1)
        acetyl_frag2 = Chem.MolFromSmarts(acetyl_smarts2)
    except Exception as e:
        return False, "Error constructing acetyl SMARTS: " + str(e)
    if acetyl_frag1 is None or acetyl_frag2 is None:
        return False, "Error constructing one of the acetyl SMARTS patterns"
    
    # Find all matches for acetyl fragment (using union of the two SMARTS patterns)
    matches_acetyl = set(mol.GetSubstructMatches(acetyl_frag1, useChirality=True))
    for match in mol.GetSubstructMatches(acetyl_frag2, useChirality=True):
        matches_acetyl.add(match)
    
    acetyl_count = len(matches_acetyl)
    if acetyl_count == 0:
        return False, "No acetyl (CH3C(=O)-) group found"
    if acetyl_count > 1:
        return False, f"Multiple acetyl groups detected ({acetyl_count})"
    
    # Define SMARTS for a free carboxyl group.
    freeCO_smarts = "C(=O)[O;H1,O-]"
    try:
        freeCO_frag = Chem.MolFromSmarts(freeCO_smarts)
    except Exception as e:
        return False, "Error constructing free carboxyl SMARTS: " + str(e)
    if freeCO_frag is None:
        return False, "Error constructing free carboxyl SMARTS pattern"
    
    freeCO_matches = mol.GetSubstructMatches(freeCO_frag)
    freeCO_count = len(freeCO_matches)
    if freeCO_count == 0:
        return False, "No free carboxyl (COOH or COO-) group found"
    if freeCO_count > 1:
        return False, f"Multiple free carboxyl groups detected ({freeCO_count}); likely a peptide or diacid"
    
    # Count total amide bonds in the molecule. We use the SMARTS "[C](=O)N" to catch any amide.
    amide_smarts = "[C](=O)N"
    try:
        amide_frag = Chem.MolFromSmarts(amide_smarts)
    except Exception as e:
        return False, "Error constructing amide SMARTS: " + str(e)
    if amide_frag is None:
        return False, "Error constructing amide SMARTS pattern"
    
    amide_matches = mol.GetSubstructMatches(amide_frag)
    amide_count = len(amide_matches)
    # In a valid N-acetyl amino acid, the only amide bond is that from the acetyl group.
    if amide_count != 1:
        return False, f"Found {amide_count} amide bonds; expected exactly 1. This likely indicates a peptide chain or extra acyl groups"
    
    return True, ("Contains a single acetyl (CH3C(=O)-) group, one free carboxyl group, and one amide bond; "
                  "consistent with an N-acetyl amino acid structure")

# Optional: run some tests when executing this module directly.
if __name__ == "__main__":
    # A set of examples from the outcomes mentioned:
    examples = [
        ("CC(=O)N[C@@H](CCCNC(N)=N)C(O)=O", "N(alpha)-acetyl-L-arginine"),                   # True positive
        ("C[C@H](NC(C)=O)C(O)=O", "N-acetyl-L-alanine"),                                        # True positive
        ("CC(=O)N1CCC[C@H]1C(O)=O", "N-acetyl-L-proline"),                                       # True positive
        ("CC(=O)N[C@@H](Cc1ccccc1)C(O)=O", "N-acetyl-L-phenylalanine"),                          # True positive
        ("CC(=O)NCCC[C@H](N)C(O)=O", "N(5)-acetyl-L-ornithine"),                                  # Expected True now
        ("C(=O)(C(CC=1NC=NC1)NC(=O)C)O", "N-acetylhistidine"),                                   # Expected True now
        # False positives (peptides, etc.)
        ("O=C(NCC(=O)N[C@@H](CO)C(O)=O)[C@@H](N)C(C)C1NC=NC1", "His-Gly-Ser"),                  # Likely peptide
        ("O=C(N[C@@H](CCCN=C(N)N)C(O)=O)[C@@H](NC(=O)[C@@H](N)CC1NC=NC1)CC1NC=NC1", "His-His-Arg"),
    ]
    
    for sm, name in examples:
        result, reason = is_N_acetyl_amino_acid(sm)
        print("SMILES:", sm)
        print("NAME:", name)
        print("RESULT:", result)
        print("REASON:", reason)
        print("-" * 60)