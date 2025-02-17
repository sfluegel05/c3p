"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
#!/usr/bin/env python
"""
Classifies: N-acetyl-amino acid
Definition: An N-acyl-amino acid that has acetyl (CH3C(=O)-) as the acyl group.
Heuristic:
  1. The molecule must contain an N–acetyl group, identified as the fragment
     [CH3]C(=O)N (i.e. a methyl group directly attached to a carbonyl which in turn is attached to a nitrogen).
  2. The molecule must contain at least one carboxyl group (protonated or deprotonated),
     identified by the fragment "[C](=O)[O;H1,O-]".
  3. To avoid classifying peptides, we count all amide bonds (matching the substructure "C(=O)N").
     A genuine N–acetyl amino acid should contain exactly one amide bond (the N–acetyl linkage).
  4. If additional amide bonds exist, the molecule is likely a peptide or contains extra acyl groups.
"""

from rdkit import Chem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N‐acetyl‐amino acid based on its SMILES string.
    
    The method uses the following heuristic:
      1. It must have an N–acetyl group (the SMARTS "[CH3]C(=O)N").
      2. It must have at least one carboxyl group (SMARTS "[C](=O)[O;H1,O-]").
      3. It must have exactly one amide bond. That is, matching the SMARTS "C(=O)N".
         Peptides (with multiple amino acid units) will have additional amide bonds.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if classified as an N–acetyl–amino acid, False otherwise.
        str: A reason explaining the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for N-acetyl group: a methyl group bound to a carbonyl which is in turn bound to nitrogen.
    acetyl_smarts = "[CH3]C(=O)N"
    acetyl_frag = Chem.MolFromSmarts(acetyl_smarts)
    if acetyl_frag is None:
        return False, "Error in constructing acetyl SMARTS"
    acetyl_matches = mol.GetSubstructMatches(acetyl_frag)
    if not acetyl_matches:
        return False, "Missing N-acetyl group (expected pattern: CH3C(=O)N)"
    
    # Check for the presence of a carboxyl group: protonated or deprotonated.
    carboxyl_smarts = "[C](=O)[O;H1,O-]"
    carboxyl_frag = Chem.MolFromSmarts(carboxyl_smarts)
    if carboxyl_frag is None:
        return False, "Error in constructing carboxyl SMARTS"
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_frag)
    if not carboxyl_matches:
        return False, "Missing carboxylic acid group"
    
    # Count all amide bonds by looking for the pattern "C(=O)N".
    amide_smarts = "C(=O)N"
    amide_frag = Chem.MolFromSmarts(amide_smarts)
    if amide_frag is None:
        return False, "Error in constructing amide SMARTS"
    amide_matches = mol.GetSubstructMatches(amide_frag)
    n_amide = len(amide_matches)
    
    # For a single N-acetyl amino acid, we expect exactly one amide bond (the N-acetyl linkage).
    if n_amide == 0:
        return False, "No amide bond (N-acetyl linkage) found"
    if n_amide > 1:
        return False, f"Multiple amide bonds detected ({n_amide}); likely a peptide chain or additional acyl groups present"
    
    # If all conditions are met we classify the molecule as an N-acetyl-amino acid.
    return True, "Contains an N-acetyl group, a carboxyl group, and exactly one amide bond consistent with a single amino acid backbone"

# Optional: Testing the function with some examples.
if __name__ == "__main__":
    test_smiles = [
        # True positives:
        "CC(=O)N[C@@H](CCCNC(N)=N)C(O)=O",      # N(alpha)-acetyl-L-arginine
        "C[C@H](NC(C)=O)C(O)=O",                # N-acetyl-L-alanine
        "CC(=O)N1CCC[C@H]1C(O)=O",               # N-acetyl-L-proline
        "CC(=O)N[C@@H](Cc1ccccc1)C(O)=O",        # N-acetyl-L-phenylalanine
        "CC(=O)NCCC[C@H](N)C(O)=O",              # N(5)-acetyl-L-ornithine
        "CC(=O)NCCC(O)=O",                      # N-acetyl-beta-alanine
        # False positives should be rejected due to extra amide bonds (peptides):
        "O=C(NCC(=O)N[C@@H](CO)C(O)=O)[C@@H](N)CC=1NC=NC1",  # a di-/tripeptide example
        # Also test a false negative scenario if no appropriate backbone exists:
        "O=C(N[C@@H](CCCN=C(N)N)C(O)=O)[C@@H](NC(=O)[C@@H](N)CC=1NC=NC1)CC=2NC=NC2",  # a tripeptide
    ]
    
    for s in test_smiles:
        res, reason = is_N_acetyl_amino_acid(s)
        print(f"SMILES: {s}\n  Classified as N-acetyl-amino acid? {res}\n  Reason: {reason}\n")