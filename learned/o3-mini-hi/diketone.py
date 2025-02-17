"""
Classifies: CHEBI:46640 diketone
"""
"""
Classifies: A compound that contains two ketone functionalities (diketone)
Using a simple SMARTS pattern that detects ketone groups: [#6][C](=O)[#6]
This script uses RDKit to read the SMILES string and identify exactly 2 ketone groups.
"""

from rdkit import Chem

def is_diketone(smiles: str):
    """
    Determines if a molecule is a diketone based on its SMILES string.
    A diketone is defined as a compound that contains two ketone functionalities,
    where a ketone group is represented as a carbonyl (C=O) with the carbon connected
    to two other carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diketone, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a ketone functionality:
    # A carbonyl group (C=O) where the carbon is bound to two carbon (non-hydrogen) atoms.
    ketone_pattern = Chem.MolFromSmarts("[#6][C](=O)[#6]")
    
    # Find substructure matches for the ketone pattern.
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    num_ketones = len(ketone_matches)
    
    # Some molecules might have overlapping matches in fused systems,
    # however for typical diketones, we expect exactly 2 distinct ketone groups.
    if num_ketones != 2:
        return False, f"Found {num_ketones} ketone functionalities; exactly 2 are required for a diketone"
    
    return True, "Contains exactly two ketone functionalities"

# Example usage:
if __name__ == "__main__":
    # Test with 2,3-Undecanedione: O=C(CCCCCCCC)C(=O)C
    test_smiles = "O=C(CCCCCCCC)C(=O)C"
    result, reason = is_diketone(test_smiles)
    print(f"Result: {result}\nReason: {reason}")