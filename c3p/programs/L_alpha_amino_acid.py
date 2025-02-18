"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
#!/usr/bin/env python
"""
Classifies: L-alpha-amino acid
Definition: Any alpha-amino acid having L-configuration at the alpha-carbon.
An L-alpha-amino acid should contain a free amino acid backbone, i.e.,
a free amino group (not amidated) attached to a chiral alpha-carbon that bears
a carboxylic acid group (with an –OH), and the computed CIP code for that carbon should be "S".
We also require that exactly one such backbone is present (to avoid peptides).
"""

from rdkit import Chem

def is_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an L-alpha-amino acid based on its SMILES string.
    Checks for:
      1. A free amino group (i.e. not amidated)
      2. A chiral alpha-carbon attached to a carboxylic acid (C(=O)O) that is explicitly protonated (has an -OH).
      3. Exactly one such amino acid backbone (to avoid peptides).
      4. The alpha-carbon has a computed CIP code of "S" (by convention corresponding to L configuration).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a free L-alpha-amino acid, False otherwise.
        str : Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Force RDKit to assign stereochemistry and CIP codes.
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)

    # We define two SMARTS patterns that capture a free amino acid backbone.
    # Explanation:
    #   - The amino nitrogen (N) must have at least one hydrogen (H1 or H2)
    #     and must not be directly bound to a carbonyl (avoiding amides).
    #   - The alpha carbon is chiral ([C@H] or [C@@H]) and bonded to some sidechain (any *)
    #   - The carboxyl group is defined as C(=O)O, ensuring that one oxygen is protonated.
    pattern1 = Chem.MolFromSmarts("[NX3;H2,H1;!$(N-C(=O))][C@H]([#6])C(=O)[O;H1]")
    pattern2 = Chem.MolFromSmarts("[NX3;H2,H1;!$(N-C(=O))][C@@H]([#6])C(=O)[O;H1]")
    
    # Find all backbone matches from both patterns.
    matches = mol.GetSubstructMatches(pattern1) + mol.GetSubstructMatches(pattern2)
    
    if not matches:
        return False, "Alpha-amino acid backbone (free amino group and free -COOH) not found"
    
    # If more than one backbone is found, the molecule is likely a peptide.
    if len(matches) > 1:
        return False, f"Found {len(matches)} amino acid backbones; likely a peptide, not a single amino acid"

    # We found exactly one match -- now check its stereochemistry.
    # Expect the pattern to assign the second atom in the SMARTS as the alpha-carbon.
    match = matches[0]
    # The SMARTS was defined so that the first atom is the amino nitrogen and the second is the alpha-carbon.
    alpha_idx = match[1]
    alpha_atom = mol.GetAtomWithIdx(alpha_idx)
    
    # Make sure the alpha-carbon has a CIP code assigned.
    if not alpha_atom.HasProp("_CIPCode"):
        return False, "Alpha-carbon has no CIP code; cannot determine configuration"
    cip = alpha_atom.GetProp("_CIPCode")
    
    if cip == "S":
        return True, "Found a free alpha-amino acid backbone with L (S CIP) configuration at the alpha-carbon"
    else:
        return False, f"Alpha-amino acid backbone found but alpha-carbon CIP code is '{cip}', not 'S'"

# Example usage:
# test_smiles = "N[C@@H](CC(=O)c1cc(O)ccc1N)C(O)=O"  # 5-hydroxy-L-kynurenine example
# result, reason = is_L_alpha_amino_acid(test_smiles)
# print(result, reason)