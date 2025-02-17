"""
Classifies: CHEBI:16180 N-acylglycine
"""
"""
Classifies: N‐acylglycine
Definition: An N‐acyl‐amino acid in which the amino acid specified is glycine.
We require that the molecule contains at least one instance of the fragment:
    R‑C(=O)‐N‐CH2‐C(=O)[O]  
with the CH2 group enforced (via [CH2;H2]) to represent glycine.
Note: The previous version rejected molecules with extra amide bonds in the whole molecule,
but many valid N‑acylglycine derivatives occur in polyamide settings.
This version simply requires that a validated N‐acylglycine fragment is found.
"""

from rdkit import Chem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.
    
    An N-acylglycine is defined as an N-acyl amino acid in which the amino acid is glycine.
    The canonical fragment we require is:
          R-C(=O)-N-CH2-C(=O)[O]
    where the CH2 is enforced (by specifying two hydrogen atoms) so that the amino acid is glycine.
    
    This implementation:
      1. Parses the SMILES and adds explicit hydrogens (so that the [CH2;H2] constraint works better).
      2. Uses a SMARTS pattern which looks for [CX3](=O)[NX3][CH2;H2][C](=O)[O]
         The meaning is:
            [CX3](=O)   : a carbonyl carbon (the acyl carbon)
            [NX3]       : the amide nitrogen (trivalent)
            [CH2;H2]    : a methylene group with exactly two hydrogens (for glycine)
            [C](=O)[O]  : a carboxyl carbon (as acid or carboxylate)
      3. Does not (anymore) count overall amide bonds because some valid derivatives have multiple.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an N-acylglycine, False otherwise.
        str: Explanation of the decision.
    """
    # Parse SMILES and add explicit hydrogens so that the [CH2;H2] constraint is effective.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Define SMARTS for the N‐acylglycine fragment:
    # [CX3](=O)[NX3][CH2;H2][C](=O)[O]
    n_acylglycine_smarts = "[CX3](=O)[NX3][CH2;H2][C](=O)[O]"
    n_acylglycine_pattern = Chem.MolFromSmarts(n_acylglycine_smarts)
    if n_acylglycine_pattern is None:
        return False, "Failed to create SMARTS pattern for N-acylglycine"
    
    # Look for any substructure match corresponding to the N‐acylglycine fragment.
    matches = mol.GetSubstructMatches(n_acylglycine_pattern)
    if not matches:
        return False, "N-acylglycine substructure not found"
    
    # For additional confidence we can inspect each match.
    # Optionally, one might check that the glycine CH2 (match index 2) is only bonded to the amide N and the carboxyl carbon.
    valid_match_found = False
    for match in matches:
        # match is a tuple of atom indices corresponding to:
        # 0: acyl carbon (C in R-C(=O))
        # 1: amide nitrogen (N)
        # 2: glycine CH2 (CH2 with two hydrogens)
        # 3: carboxyl carbon (C in C(=O)[O])
        glycine_atom = mol.GetAtomWithIdx(match[2])
        # Count heavy neighbors that are not hydrogen.
        heavy_neighbors = [nbr for nbr in glycine_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        # For a true glycine CH2, the heavy neighbors should be exactly 2 (the amide nitrogen and the carboxyl carbon).
        if len(heavy_neighbors) == 2:
            valid_match_found = True
            break
            
    if not valid_match_found:
        return False, "N-acylglycine substructure found but the glycine CH2 environment is ambiguous."
    
    # If a valid match is observed, classify as N-acylglycine.
    return True, "Molecule contains the N-acylglycine substructure (R-C(=O)-N-CH2-C(=O)[O])"

# Example usage (for testing):
# test_smiles = "CC(=O)NCC(O)=O"  # N-acetylglycine; expected True.
# result, reason = is_N_acylglycine(test_smiles)
# print(result, reason)