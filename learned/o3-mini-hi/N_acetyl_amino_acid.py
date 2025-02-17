"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
#!/usr/bin/env python
"""
Classifies: N-acetyl-amino acid
Definition: An N-acyl-amino acid where the acyl group is specifically acetyl (CH3C(=O)-)
Heuristic:
  1. The molecule must contain an N-acetyl substructure. We identify this by matching the SMARTS pattern "[NX3]C(=O)C".
  2. The molecule must also contain a carboxylic acid group, matched by "[CX3](=O)[O;H1,O-]".
  3. To increase confidence that the N-acetyl group is on an amino acid moiety, we require that the nitrogen of the N-acetyl group is connected (aside from the acetyl carbon) to a candidate α‑carbon. That α‑carbon must be directly bonded to the carboxyl carbon (from the carboxylic acid group).
Note:
  This heuristic may not work for every edge case but should correctly classify most structures in the given list.
"""

from rdkit import Chem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.
    
    The function uses heuristics:
      - Checks for the presence of an N-acetyl group (N attached to a C(=O)CH3 fragment).
      - Checks for the presence of a carboxyl group (C(=O)O or C(=O)[O-]).
      - Verifies that the acetylated nitrogen is connected to an alpha carbon that is directly bonded to the carboxyl carbon.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an N-acetyl-amino acid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Define SMARTS for N-acetyl group: a nitrogen bound to a carbonyl carbon that has a methyl group.
    acetyl_smarts = Chem.MolFromSmarts("[NX3]C(=O)C")
    if acetyl_smarts is None:
        return False, "Error in acetyl SMARTS"
    
    if not mol.HasSubstructMatch(acetyl_smarts):
        return False, "Missing N-acetyl group (expected pattern: N–C(=O)C)"
    
    # Define SMARTS for a carboxyl group (protonated or deprotonated)
    carboxyl_smarts = Chem.MolFromSmarts("[CX3](=O)[O;H1,O-]")
    if carboxyl_smarts is None:
        return False, "Error in carboxyl SMARTS"
    
    if not mol.HasSubstructMatch(carboxyl_smarts):
        return False, "Missing carboxylic acid group"
    
    # Retrieve matches for the N-acetyl pattern and carboxyl pattern
    acetyl_matches = mol.GetSubstructMatches(acetyl_smarts)
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_smarts)
    
    # Build a set of carboxyl carbon atom indices from the carboxyl matches.
    # In the SMARTS, the first atom (index 0) is the carbonyl carbon.
    carboxyl_carbons = {match[0] for match in carboxyl_matches}
    
    # For each N-acetyl match, check connectivity: the N should be attached to an "alpha" carbon,
    # and that alpha carbon should be directly bonded to one of the carboxyl carbon(s).
    for match in acetyl_matches:
        acetyl_N = match[0]  # index of the nitrogen in the acetyl group
        acetyl_C = match[1]  # acetyl carbon (the carbonyl carbon)
        # Get the nitrogen atom from the molecule
        N_atom = mol.GetAtomWithIdx(acetyl_N)
        # Iterate over neighbors of N_atom:
        for neighbor in N_atom.GetNeighbors():
            # Exclude the acetyl carbon already identified
            if neighbor.GetIdx() == acetyl_C:
                continue
            # Consider this neighbor as a candidate for the alpha carbon.
            alpha_idx = neighbor.GetIdx()
            alpha_atom = neighbor
            # Check if the candidate alpha carbon has a neighbor that is one of the carboxyl carbons.
            for alpha_neighbor in alpha_atom.GetNeighbors():
                if alpha_neighbor.GetIdx() in carboxyl_carbons:
                    # Found an acetylated N attached to an alpha carbon that is bound to a carboxyl carbon.
                    return True, "Contains N-acetyl group on an amino acid backbone (alpha carbon connected to carboxyl group)"
    
    # If no connectivity fulfilling the heuristic is found, return false.
    return False, "N-acetyl group found but not appropriately connected to an amino acid (alpha carbon with carboxyl group)"

# (Optional) Example usage for testing:
if __name__ == "__main__":
    test_smiles = [
        "CC(=O)N[C@@H](CCCNC(N)=N)C(O)=O",  # N(alpha)-acetyl-L-arginine
        "C[C@H](NC(C)=O)C(O)=O",            # N-acetyl-L-alanine
        "CC(=O)N1CCC[C@H]1C(O)=O",           # N-acetyl-L-proline
        "CC(=O)NCCCC[C@H](N)C(O)=O",         # N(5)-acetyl-L-ornithine
        "CC(=O)NC(Cc1ccccc1)C(O)=O",         # N-acetyl-L-phenylalanine
        "CCCC",                            # A molecule that is not an amino acid
    ]
    for s in test_smiles:
        res, reason = is_N_acetyl_amino_acid(s)
        print(f"SMILES: {s}\n  Classified as N-acetyl-amino acid? {res}\n  Reason: {reason}\n")