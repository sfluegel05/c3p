"""
Classifies: CHEBI:76983 N-sulfonylurea
"""
"""
Classifies: N-sulfonylurea
Definition: An N-sulfonylurea is defined as a urea (–NH–C(=O)–NH–) in which one 
of the hydrogen atoms on a urea nitrogen is replaced by a sulfonyl group (–S(=O)(=O)–).
This motif is common in many herbicides and type 2 diabetes medications.
"""

from rdkit import Chem

def is_N_sulfonylurea(smiles: str):
    """
    Determines if a molecule is an N-sulfonylurea based on its SMILES string.
    
    Our strategy is to search for a urea group (N-C(=O)-N) in which one of the nitrogens
    is substituted with a sulfonyl group (S(=O)(=O)). We use two SMARTS patterns to cover 
    the two possibilities:
      Pattern 1 (left substitution): "N(-S(=O)(=O))C(=O)N"
      Pattern 2 (right substitution): "NC(=O)N(-S(=O)(=O))"
      
    To reduce false positives we also check that the substituted nitrogen shows the 
    appropriate hydrogen count: A regular (neutral) urea nitrogen would normally have two H’s;
    after substitution it should have one—and in some ionized cases it may lack an explicit hydrogen 
    if its formal charge is -1.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an N-sulfonylurea, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define two SMARTS patterns covering both possibilities:
    # Pattern 1: sulfonyl group on the left nitrogen: N(-S(=O)(=O))C(=O)N
    pattern1 = Chem.MolFromSmarts("N(-S(=O)(=O))C(=O)N")
    # Pattern 2: sulfonyl group on the right nitrogen: NC(=O)N(-S(=O)(=O))
    pattern2 = Chem.MolFromSmarts("NC(=O)N(-S(=O)(=O))")
    
    # Helper function: check if the substituted N (attached to S) has the proper H count.
    # For neutral nitrogen that lost one hydrogen it should have exactly one hydrogen.
    # For an ionized nitrogen (formal charge -1) it might be missing an explicit hydrogen.
    def valid_modified_n(atom):
        numHs = atom.GetTotalNumHs()  # implicit+explicit hydrogen count
        charge = atom.GetFormalCharge()
        # Accept if exactly one hydrogen OR (zero H and formal charge -1)
        return (numHs == 1) or (numHs == 0 and charge == -1)
    
    match_found = False
    reason_details = ""
    
    # Check pattern 1: in this match, the order is [N_modified, C, N] where the first N
    # is the one bonded to the sulfonyl group.
    matches1 = mol.GetSubstructMatches(pattern1)
    for match in matches1:
        # match[0] is the substituted nitrogen (the one with S(=O)(=O) attached)
        sub_n = mol.GetAtomWithIdx(match[0])
        if valid_modified_n(sub_n):
            match_found = True
            reason_details = "Molecule contains a valid N-sulfonylurea motif (pattern: N(-S(=O)(=O))C(=O)N)."
            break
    
    # If not found with pattern 1, try pattern 2.
    if not match_found:
        matches2 = mol.GetSubstructMatches(pattern2)
        for match in matches2:
            # For pattern 2, the substituted nitrogen is the second one in the pattern,
            # so we use match[2] (0-based indexing: atoms from NC(=O)N(-S(=O)(=O)) ).
            sub_n = mol.GetAtomWithIdx(match[2])
            if valid_modified_n(sub_n):
                match_found = True
                reason_details = "Molecule contains a valid N-sulfonylurea motif (pattern: NC(=O)N(-S(=O)(=O)))."
                break

    if match_found:
        return True, reason_details
    else:
        return False, "N-sulfonylurea moiety not found in the molecule."

# Example usage (remove when using as a module):
if __name__ == "__main__":
    # Test with glyburide SMILES as an example (should return True)
    test_smiles = "COc1ccc(Cl)cc1C(=O)NCCc1ccc(cc1)S(=O)(=O)NC(=O)NC1CCCCC1"
    result, reason = is_N_sulfonylurea(test_smiles)
    print("Result:", result)
    print("Reason:", reason)