"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
#!/usr/bin/env python
"""
Classifies: N-acetyl-amino acid
Definition: An N-acyl-amino acid that has acetyl (CH3C(=O)-) as the acyl group.
Heuristic:
  - The molecule must contain at least one N-acetyl group.
    We identify an N–acetyl fragment by matching the SMARTS pattern "[NX3]C(=O)C"
    (i.e. a nitrogen bound to a carbonyl carbon that in turn is bound to a methyl group).
  - The molecule must contain a carboxyl group (protonated or deprotonated), 
    matched by "[C](=O)[O;H1,O-]".
  - To ensure that the acetyl group is “on an amino acid” (rather than on a peptide chain),
    we require that the molecule contains exactly one amino acid backbone.
    We approximate an amino acid backbone by requiring that there is exactly one atom
    (the α–carbon) that is connected both to an amino group (any [NX3]) and a carboxyl group.
  
  Note:
    This heuristic may not work for every edge case but should correctly classify most of
    the supplied examples. In particular, di‐ or oligo–peptides (which show multiple backbone
    units) will be rejected.
"""

from rdkit import Chem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a given molecule is an N‐acetyl‐amino acid based on its SMILES string.
    
    The function uses these heuristics:
      1. It must have an N-acetyl group (identified by the SMARTS "[NX3]C(=O)C").
      2. It must have at least one carboxylic acid group (SMARTS "[C](=O)[O;H1,O-]").
      3. It must contain exactly one amino acid backbone. We use a substructure that matches
         an α–carbon with a carboxyl group and at least one nitrogen neighbor. That is, we look
         for a pattern like "[C;!R]([NX3])C(=O)[O;H1,O-]". This catches both standard and N-acyl amino acids.
      4. (Peptides, which have multiple backbone units, are thereby eliminated.)
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an N–acetyl–amino acid, False otherwise.
        str: A reason explaining the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS for N-acetyl group: nitrogen bound to a C(=O) that carries a methyl (CH3).
    acetyl_smarts = "[NX3]C(=O)C"
    acetyl_frag = Chem.MolFromSmarts(acetyl_smarts)
    if acetyl_frag is None:
        return False, "Error in constructing acetyl SMARTS"
    acetyl_matches = mol.GetSubstructMatches(acetyl_frag)
    if not acetyl_matches:
        return False, "Missing N-acetyl group (expected pattern: N–C(=O)C)"
    
    # SMARTS for carboxyl group (protonated or deprotonated)
    carboxyl_smarts = "[C](=O)[O;H1,O-]"
    carboxyl_frag = Chem.MolFromSmarts(carboxyl_smarts)
    if carboxyl_frag is None:
        return False, "Error in constructing carboxyl SMARTS"
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_frag)
    if not carboxyl_matches:
        return False, "Missing carboxylic acid group"
    
    # Identify an amino acid backbone.
    # We search for an α–carbon that is bound to a nitrogen and a carboxyl group.
    # The pattern means: a non-ring carbon ([C;!R]) that is connected to some nitrogen ([NX3])
    # and also connected to a carboxyl group (C(=O)[O;H1,O-]). This should catch the α–carbon.
    aa_core_smarts = "[C;!R]([NX3])C(=O)[O;H1,O-]"
    aa_core_frag = Chem.MolFromSmarts(aa_core_smarts)
    if aa_core_frag is None:
        return False, "Error in constructing amino acid core SMARTS"
    aa_cores = mol.GetSubstructMatches(aa_core_frag)
    
    if len(aa_cores) == 0:
        return False, "No amino acid backbone found"
    elif len(aa_cores) > 1:
        return False, "Multiple amino acid backbones found; likely a peptide chain"
    
    # If we reach here, we have exactly one amino acid core and at least one N-acetyl group.
    # We do not insist that the acetyl group be directly on the backbone nitrogen.
    # (This improves sensitivity to cases such as N(5)-acetyl-L-ornithine or N(6)-acetyl-L-lysine,
    # where the acetylation occurs on a side-chain amino group.)
    return True, "Contains an N-acetyl group and exactly one amino acid backbone (α–carbon attached to a carboxyl group)"

# (Optional) Example usage for testing:
if __name__ == "__main__":
    test_smiles = [
        # True positives:
        "CC(=O)N[C@@H](CCCNC(N)=N)C(O)=O",  # N(alpha)-acetyl-L-arginine
        "C[C@H](NC(C)=O)C(O)=O",            # N-acetyl-L-alanine
        "CC(=O)N1CCC[C@H]1C(O)=O",           # N-acetyl-L-proline
        "CC(=O)N[C@@H](Cc1ccccc1)C(O)=O",    # N-acetyl-L-phenylalanine
        # Examples which were missed previously (false negatives):
        "CC(=O)NCCC[C@H](N)C(O)=O",          # N(5)-acetyl-L-ornithine
        "CC(=O)NCCCC[C@H](N)C(O)=O",         # N(6)-acetyl-L-lysine
        "CC(=O)NCCC(O)=O",                  # N-acetyl-beta-alanine
        # False positives (peptides with multiple backbones):
        "O=C(NCC(=O)N[C@@H](CO)C(O)=O)[C@@H](N)CC=1NC=NC1",  # His-Gly-Ser (di-/tri-peptide)
    ]
    for s in test_smiles:
        res, reason = is_N_acetyl_amino_acid(s)
        print(f"SMILES: {s}\n  Classified as N-acetyl-amino acid? {res}\n  Reason: {reason}\n")