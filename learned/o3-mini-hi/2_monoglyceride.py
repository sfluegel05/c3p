"""
Classifies: CHEBI:17389 2-monoglyceride
"""
"""
Classifies: 2-monoglyceride
A 2-monoglyceride is defined as a monoglyceride in which the acyl substituent is located at position 2.
In other words, the glycerol backbone (normally HOCH2-CHOH-CH2OH) has its central hydroxyl replaced by an ester function.
We use a SMARTS pattern that requires a fragment with the connectivity:
    [CH2][OH]-[C;H1](O[C](=O)[*])-[CH2][OH]
which corresponds to a terminal CH2OH, a central CH that bears an ester (an O bonded to a C(=O) and an arbitrary R group), 
and another CH2OH.
"""

from rdkit import Chem

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    
    A 2-monoglyceride contains a glycerol backbone with an acyl (ester) group attached at the secondary (2-) position.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a 2-monoglyceride, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern to look for the key fragment.
    # The pattern looks for: CH2OH - CH(O-C(=O)[*]) - CH2OH, that is
    # a glycerol-like fragment in which the central carbon has an ester substituent.
    pattern = Chem.MolFromSmarts("[CH2][OH]-[C;H1](O[C](=O)[*])-[CH2][OH]")
    if pattern is None:
        return False, "Failed to generate SMARTS pattern"
    
    # Find all matches for the pattern.
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "No 2-monoglyceride (sn-2 ester on glycerol) fragment found"
    
    # In a true monoglyceride we expect one instance of the glycerol backbone with a single acyl group.
    if len(matches) > 1:
        return False, f"Found {len(matches)} glycerol ester fragments; expected exactly one for a monoglyceride"
    
    # Optionally, one could add further checks (e.g., ensuring no additional ester groups on the glycerol backbone).
    # For this example, matching the pattern is taken as sufficient.
    return True, "Contains a glycerol backbone with an acyl ester at the 2-position (2-monoglyceride)"

# Example usage:
if __name__ == "__main__":
    # A list of test SMILES (selected examples from the prompt):
    test_smiles = [
        "O(C(CO)CO)C(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC",  # MG(0:0/22:6(4Z,7Z,10Z,13Z,16Z,19Z)/0:0)
        "CCCCCCCCCCCCCCCC(=O)OC(CO)CO",  # 2-palmitoylglycerol
        "O(C(CO)CO)C(=O)CCCCCCCCCCCCCC",  # another example
        "CCO",  # obviously not a 2-monoglyceride
    ]
    
    for s in test_smiles:
        valid, reason = is_2_monoglyceride(s)
        print(f"SMILES: {s}\nResult: {valid}, Reason: {reason}\n")