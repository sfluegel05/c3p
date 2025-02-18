"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
#!/usr/bin/env python
"""
Classifies: L-alpha-amino acid
Definition: Any alpha-amino acid having L-configuration at the alpha-carbon.
A free L-alpha-amino acid is defined here as a molecule with a free amino acid backbone,
i.e., an unmodified (non-amidated) amino group attached to a chiral alpha-carbon that bears
a protonated carboxylic acid group (-C(=O)O with an –OH), and with exactly one such backbone.
The alpha-carbon is assumed to have the CIP code "S" corresponding to the L configuration.
"""

from rdkit import Chem

def is_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is a free L-alpha-amino acid based on its SMILES string.
    Checks for:
      1. A free amino group that is not amidated.
      2. A chiral alpha-carbon attached to a protonated carboxylic acid (C(=O)[O;H1]).
      3. Exactly one such amino acid backbone in the molecule (to avoid peptides).
      4. The alpha-carbon has a CIP code of "S" (i.e. L-configuration).
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a free L-alpha-amino acid, False otherwise.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Assign stereochemistry and compute CIP codes.
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Define SMARTS patterns for the amino acid backbone:
    # The pattern looks for an amino group [NX3] (with one or two hydrogens) that is NOT bound
    # to a carbonyl (to ensure it is free and not in an amide). The amino is bonded to a chiral alpha-carbon,
    # which is in turn bonded to a carboxylic acid group (C(=O)[O;H1] ensures that one oxygen has a hydrogen).
    # We define two patterns to capture either possible chiral notations (@ or @@). In both,
    # the ordering is: atom0: amino nitrogen, atom1: alpha-carbon, atom2: carbonyl carbon.
    pattern1 = Chem.MolFromSmarts("[NX3;H2,H1;!$(N-C(=O))]-[C@H]([#6])-[C](=O)[O;H1]")
    pattern2 = Chem.MolFromSmarts("[NX3;H2,H1;!$(N-C(=O))]-[C@@H]([#6])-[C](=O)[O;H1]")
    
    # Get substructure matches (using useChirality so that pattern matching respects stereochemistry).
    matches1 = mol.GetSubstructMatches(pattern1, useChirality=True)
    matches2 = mol.GetSubstructMatches(pattern2, useChirality=True)
    all_matches = matches1 + matches2
    
    # Deduplicate matches by converting each match (a tuple of atom indices) into a sorted tuple.
    unique_matches = set()
    for m in all_matches:
        unique_matches.add(tuple(sorted(m)))
    
    if not unique_matches:
        return False, "Alpha-amino acid backbone (free amino group and free carboxylic acid) not found"
    
    if len(unique_matches) > 1:
        return False, f"Found {len(unique_matches)} amino acid backbones; likely a peptide, not a single amino acid"
    
    # We now have exactly one unique backbone. We need to check the chirality of the alpha-carbon.
    # We know that in our SMARTS, the second atom (index 1 in the match) is the chiral alpha-carbon.
    match = list(unique_matches)[0]
    # To preserve the SMARTS ordering, we retrieve the match using the SMARTS match from pattern1 if possible.
    # Otherwise, we use an arbitrary ordering from unique_matches.
    # NOTE: Since our patterns were defined with three atoms, we assume the center atom (alpha-carbon) is the second in our match.
    # Because we deduplicated by sorted order, we must find the corresponding match from the original list.
    target_alpha_idx = None
    for m in all_matches:
        if tuple(sorted(m)) == match:
            target_alpha_idx = m[1]  # second atom in the ordering from the SMARTS
            break
    if target_alpha_idx is None:
        return False, "Alpha-amino acid backbone detected but could not retrieve alpha-carbon information"
    
    alpha_atom = mol.GetAtomWithIdx(target_alpha_idx)
    if not alpha_atom.HasProp("_CIPCode"):
        return False, "Alpha-carbon has no CIP code; cannot determine configuration"
    cip = alpha_atom.GetProp("_CIPCode")
    if cip == "S":
        return True, "Found a free alpha-amino acid backbone with L (S CIP) configuration at the alpha-carbon"
    else:
        return False, f"Alpha-amino acid backbone found but alpha-carbon CIP code is '{cip}', not 'S'"

# Example usage:
# test_smiles = "N[C@@H](CC(=O)c1cc(O)ccc1N)C(O)=O"  # 5-hydroxy-L-kynurenine
# result, reason = is_L_alpha_amino_acid(test_smiles)
# print(result, reason)