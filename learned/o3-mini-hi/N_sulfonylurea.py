"""
Classifies: CHEBI:76983 N-sulfonylurea
"""
"""
Classifies: N-sulfonylurea
Definition: A urea in which one of the hydrogens attached to a nitrogen of the urea group is replaced by a sulfonyl group.
The N-sulfonylurea moiety is found in several herbicides and antidiabetic drugs.
"""

from rdkit import Chem

def is_N_sulfonylurea(smiles: str):
    """
    Determines if a molecule is an N-sulfonylurea based on its SMILES string.
    It searches for a urea motif (–NH–C(=O)–NH–) that is substituted on one nitrogen with a sulfonyl group (-S(=O)(=O)-).
    Two substructure SMARTS are used:
      Pattern 1: S(=O)(=O)[NX3]C(=O)[NX3]
      Pattern 2: [NX3]C(=O)[NX3]S(=O)(=O)
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as an N-sulfonylurea, False otherwise.
        str: Reason for the classification
    """
    # Parse the SMILES string to create a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for N-sulfonylurea motifs.
    # Pattern 1: sulfonyl group attached to the first nitrogen of a urea group
    smarts1 = "S(=O)(=O)[NX3]C(=O)[NX3]"
    # Pattern 2: sulfonyl group attached to the second nitrogen of a urea group
    smarts2 = "[NX3]C(=O)[NX3]S(=O)(=O)"
    
    pattern1 = Chem.MolFromSmarts(smarts1)
    pattern2 = Chem.MolFromSmarts(smarts2)
    
    # Check if either of the two patterns is found in the molecule.
    if mol.HasSubstructMatch(pattern1):
        return True, "Molecule contains an N-sulfonylurea motif (pattern: S(=O)(=O)[NX3]C(=O)[NX3])."
    elif mol.HasSubstructMatch(pattern2):
        return True, "Molecule contains an N-sulfonylurea motif (pattern: [NX3]C(=O)[NX3]S(=O)(=O))."
    
    # If neither pattern is found, the molecule is not an N-sulfonylurea.
    return False, "N-sulfonylurea moiety not found in the molecule."

# Example usage (you can remove these lines when using it as a module):
if __name__ == "__main__":
    # Test with glyburide SMILES as an example (it should return True)
    test_smiles = "COc1ccc(Cl)cc1C(=O)NCCc1ccc(cc1)S(=O)(=O)NC(=O)NC1CCCCC1"
    result, reason = is_N_sulfonylurea(test_smiles)
    print("Result:", result)
    print("Reason:", reason)