"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
#!/usr/bin/env python
"""
Classifies: CHEBI:??? N-acetyl-amino acid
Definition: An N-acyl-amino acid that has acetyl (CH3C(=O)-) as the acyl group.
Heuristic:
  1. Exactly one occurrence of an acetyl group attached to the nitrogen.
  2. Exactly one free (terminal) carboxyl group (pattern "C(=O)[O;H1,O-]").
  3. Exactly one amide bond overall (pattern "[C](=O)N"); extra amide bonds indicate a peptide.
"""

from rdkit import Chem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl amino acid based on its SMILES string.
    
    Enforces:
      1. The molecule must have exactly one acetyl fragment attached to nitrogen. 
         Two SMARTS patterns are employed: "[CH3]C(=O)N" and "N[C](=O)[CH3]". We canonicalize
         the atom indices for each match (via sorting) to ensure the same fragment is not counted twice.
      2. Exactly one free carboxyl group, represented by the SMARTS "C(=O)[O;H1,O-]".
      3. Only one amide bond overall (SMARTS "[C](=O)N"). Additional amide bonds would indicate a peptide.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the structure is consistent with an N-acetyl amino acid.
        str: Explanation of the decision.
    """
    # Parse the SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for the acetyl fragment. Two patterns are used because the ordering 
    # of atoms in the SMILES might vary.
    acetyl_smarts1 = "[CH3]C(=O)N"
    acetyl_smarts2 = "N[C](=O)[CH3]"
    try:
        acetyl_frag1 = Chem.MolFromSmarts(acetyl_smarts1)
        acetyl_frag2 = Chem.MolFromSmarts(acetyl_smarts2)
    except Exception as e:
        return False, "Error constructing acetyl SMARTS: " + str(e)
    if acetyl_frag1 is None or acetyl_frag2 is None:
        return False, "Error constructing one of the acetyl SMARTS patterns"
    
    # Use a set to store unique acetyl fragment matches by canonicalizing the atom indices (sorted tuples).
    unique_acetyl_matches = set()
    for match in mol.GetSubstructMatches(acetyl_frag1, useChirality=True):
        unique_acetyl_matches.add(tuple(sorted(match)))
    for match in mol.GetSubstructMatches(acetyl_frag2, useChirality=True):
        unique_acetyl_matches.add(tuple(sorted(match)))
    
    acetyl_count = len(unique_acetyl_matches)
    if acetyl_count == 0:
        return False, "No acetyl (CH3C(=O)-) group found"
    if acetyl_count > 1:
        return False, f"Multiple acetyl groups detected ({acetyl_count})"
    
    # Define SMARTS for a free carboxyl group: terminal carboxylic acid (COOH) or carboxylate (COO-)
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
    
    # Define SMARTS for an amide bond which is expected to be exactly one (the one from acetylation).
    amide_smarts = "[C](=O)N"
    try:
        amide_frag = Chem.MolFromSmarts(amide_smarts)
    except Exception as e:
        return False, "Error constructing amide SMARTS: " + str(e)
    if amide_frag is None:
        return False, "Error constructing amide SMARTS pattern"
    
    amide_matches = mol.GetSubstructMatches(amide_frag)
    amide_count = len(amide_matches)
    if amide_count != 1:
        return False, f"Found {amide_count} amide bonds; expected exactly 1. This likely indicates a peptide chain or extra acyl groups"
    
    return True, ("Contains a single acetyl (CH3C(=O)-) group, one free carboxyl group, and one amide bond; "
                  "consistent with an N-acetyl amino acid structure")

# Optional testing block
if __name__ == "__main__":
    # A few examples to test:
    examples = [
        ("CC(=O)N[C@@H](CCCNC(N)=N)C(O)=O", "N(alpha)-acetyl-L-arginine"),
        ("C[C@H](NC(C)=O)C(O)=O", "N-acetyl-L-alanine"),
        ("CC(=O)N1CCC[C@H]1C(O)=O", "N-acetyl-L-proline"),
        ("CC(=O)N[C@@H](Cc1ccccc1)C(O)=O", "N-acetyl-L-phenylalanine"),
        ("CC(=O)NCCC[C@H](N)C(O)=O", "N(5)-acetyl-L-ornithine"),
    ]
    
    for sm, name in examples:
        result, reason = is_N_acetyl_amino_acid(sm)
        print("SMILES:", sm)
        print("NAME:", name)
        print("RESULT:", result)
        print("REASON:", reason)
        print("-" * 60)