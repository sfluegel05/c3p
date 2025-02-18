"""
Classifies: CHEBI:16180 N-acylglycine
"""
"""
Classifies: N‐acylglycine
Definition: An N‐acyl‐amino acid in which the amino acid specified is glycine.
Canonical fragment: R-C(=O)-N-CH2-C(=O)[O]
We require that:
  • There is an acyl carbonyl group attached to an amide nitrogen.
  • The amide nitrogen is attached to a glycine –CH2– group,
    which in turn is attached to a carboxyl carbon.
  • The amide nitrogen must only be attached to the acyl carbon and the glycine carbon.
  • The glycine carbon is a CH2 (exactly 2 hydrogens) and has no extra heavy substitutes.
"""

from rdkit import Chem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.
    
    An N-acylglycine is an N-acyl amino acid where the amino acid is specifically glycine.
    The required fragment is:
          R-C(=O)-N-CH2-C(=O)[O]
    with the CH2 strictly enforced (exactly 2 hydrogen atoms) and proper connectivity.
    
    Process:
      1. Parse the SMILES and add explicit hydrogens.
      2. Search for the fragment using a SMARTS pattern:
             [CX3](=O)[NX3][CH2][CX3](=O)[O]
         which represents:
             • [CX3](=O)  : the acyl carbonyl group.
             • [NX3]      : the amide nitrogen.
             • [CH2]      : the glycine methylene (enforcing CH2, regardless of chiral tag).
             • [CX3](=O)[O] : the carboxyl group.
      3. For each match, verify:
             a. The amide nitrogen (match atom 1) is connected only to the acyl carbon and the glycine carbon.
             b. The glycine carbon (match atom 2) has exactly two heavy neighbors (the amide nitrogen and the carboxyl carbon).
             c. The glycine carbon indeed has exactly 2 hydrogens (using GetTotalNumHs).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an N-acylglycine, False otherwise
        str: Explanation of the decision.
    """
    # Parse SMILES and add explicit hydrogens for accurate hydrogen count.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Define SMARTS for the N-acylglycine fragment.
    # Here, we use [CH2] to ensure we get a methylene group.
    # Note: the final [O] in the carboxyl is expected to capture either –OH or –O^(-)
    n_acylglycine_smarts = "[CX3](=O)[NX3][CH2][CX3](=O)[O]"
    pattern = Chem.MolFromSmarts(n_acylglycine_smarts)
    if pattern is None:
        return False, "Failed to create SMARTS pattern for N-acylglycine"
    
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "N-acylglycine substructure not found"
    
    # Examine each match for proper connectivity.
    # Match tuple indices:
    # 0: acyl carbon; 1: amide nitrogen; 2: glycine CH2 carbon; 3: carboxyl carbon.
    for match in matches:
        acyl_c = mol.GetAtomWithIdx(match[0])
        amide_n = mol.GetAtomWithIdx(match[1])
        glycine_c = mol.GetAtomWithIdx(match[2])
        carboxyl_c = mol.GetAtomWithIdx(match[3])
        
        # For the amide nitrogen, expect exactly 2 heavy neighbors: 
        # one is the acyl carbon, one is the glycine carbon.
        n_heavy_neighbors = [nbr for nbr in amide_n.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(n_heavy_neighbors) != 2:
            continue
        if not (acyl_c in n_heavy_neighbors and glycine_c in n_heavy_neighbors):
            continue
        
        # For the glycine carbon: expect exactly 2 heavy neighbors:
        # one from the amide nitrogen and one from the carboxyl carbon.
        glycine_heavy_neighbors = [nbr for nbr in glycine_c.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(glycine_heavy_neighbors) != 2:
            continue
        if not (amide_n in glycine_heavy_neighbors and carboxyl_c in glycine_heavy_neighbors):
            continue
        
        # Verify that the glycine carbon actually has exactly 2 hydrogens.
        if glycine_c.GetTotalNumHs() != 2:
            continue
        
        return True, "Molecule contains the N-acylglycine substructure (R-C(=O)-N-CH2-C(=O)[O])"
    
    return False, "N-acylglycine substructure found but the key atom connectivity is ambiguous."
    
# Example usage (can be uncommented for testing):
# test_smiles = "CC(=O)NCC(O)=O"  # Expected True for N-acetylglycine.
# result, reason = is_N_acylglycine(test_smiles)
# print(result, reason)