"""
Classifies: CHEBI:16180 N-acylglycine
"""
"""
Classifies: N‐acylglycine
Definition: An N‐acyl‐amino acid in which the amino acid specified is glycine.
We require that the molecule contains at least one instance of the fragment:
    R‑C(=O)‐N‐CH2‐C(=O)[O]
with the CH2 enforced to have exactly two hydrogens.
We further check that the amide nitrogen is only connected to the expected two heavy atoms
(its acyl carbon and the glycine carbon) so that polyamide sequences are not wrongly classified.
"""

from rdkit import Chem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.
    
    An N-acylglycine is defined as an N-acyl amino acid in which the amino acid is glycine.
    The canonical fragment required is:
          R-C(=O)-N-CH2-C(=O)[O]
    where the methylene (CH2) is enforced by specifying that the carbon has exactly 2 hydrogens.
    
    This implementation:
      1. Parses the SMILES and adds explicit hydrogens.
      2. Uses a SMARTS pattern that looks for the fragment:
             [CX3](=O)[NX3][#6;X4;H2][CX3](=O)[O]
         which represents:
             - [CX3](=O)   : the acyl carbonyl group (R-C(=O))
             - [NX3]       : the amide nitrogen 
             - [#6;X4;H2]  : a tetrahedral carbon with exactly two hydrogen atoms 
                             (enforcing the glycine CH2 environment, regardless of chiral tag)
             - [CX3](=O)[O]: the carboxyl group.
      3. Examines each match to ensure that:
             a. The amide nitrogen is connected only to the acyl carbon and the glycine carbon.
             b. The glycine methylene carbon is connected only to the amide nitrogen and the carboxyl carbon.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an N-acylglycine, False otherwise.
        str: Explanation of the decision.
    """
    # Parse SMILES and add explicit hydrogens (for accurate hydrogen counting)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Define SMARTS for the N-acylglycine fragment:
    # [CX3](=O)[NX3][#6;X4;H2][CX3](=O)[O]
    n_acylglycine_smarts = "[CX3](=O)[NX3][#6;X4;H2][CX3](=O)[O]"
    n_acylglycine_pattern = Chem.MolFromSmarts(n_acylglycine_smarts)
    if n_acylglycine_pattern is None:
        return False, "Failed to create SMARTS pattern for N-acylglycine"
    
    matches = mol.GetSubstructMatches(n_acylglycine_pattern)
    if not matches:
        return False, "N-acylglycine substructure not found"
    
    valid_match_found = False
    for match in matches:
        # The match tuple atoms correspond to:
        # 0: acyl carbon (C of R-C(=O))
        # 1: amide nitrogen (N)
        # 2: glycine methylene carbon (should be CH2 with exactly 2 H's)
        # 3: carboxyl carbon (C of C(=O)[O])
        acyl_c = mol.GetAtomWithIdx(match[0])
        amide_n = mol.GetAtomWithIdx(match[1])
        glycine_c = mol.GetAtomWithIdx(match[2])
        carboxyl_c = mol.GetAtomWithIdx(match[3])
        
        # Check the heavy-neighbor connectivity:
        # For the amide nitrogen, expect exactly 2 heavy neighbors: 
        #   one should be the acyl carbon and one the glycine carbon.
        n_heavy_neighbors = [nbr for nbr in amide_n.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(n_heavy_neighbors) != 2:
            # Likely part of a larger peptide chain
            continue

        # For the glycine carbon, expect exactly 2 heavy neighbors:
        #   one from the amide nitrogen and one from the carboxyl carbon.
        glycine_heavy_neighbors = [nbr for nbr in glycine_c.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(glycine_heavy_neighbors) != 2:
            continue

        valid_match_found = True
        break
        
    if not valid_match_found:
        return False, "N-acylglycine substructure found but the key atom connectivity is ambiguous."

    return True, "Molecule contains the N-acylglycine substructure (R-C(=O)-N-CH2-C(=O)[O])"

# Example usage for testing (uncomment if needed):
# test_smiles = "CC(=O)NCC(O)=O"  # N-acetylglycine; expected True.
# result, reason = is_N_acylglycine(test_smiles)
# print(result, reason)