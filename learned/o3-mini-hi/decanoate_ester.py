"""
Classifies: CHEBI:87658 decanoate ester
"""
"""
Classifies: decanoate ester
Definition: A fatty acid ester resulting from the condensation of decanoic acid with an alcohol or phenol.
"""
from rdkit import Chem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    A decanoate ester must have a decanoate (capric acid) acyl group—represented by a ten‐carbon fragment in the form
    of an ester linkage (either O–C(=O)CCCCCCCCC or its reverse orientation).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains a decanoate ester group, False otherwise.
        str: Reason for classification.
    """
    # Parse the input SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
  
    # Create two SMARTS patterns to account for both orientations of the decanoate ester group:
    # Pattern 1: Ester oxygen precedes the carbonyl (commonly seen in many esters)
    pattern1 = Chem.MolFromSmarts("OC(=O)CCCCCCCCC")
    # Pattern 2: Carbonyl precedes the oxygen (the acyl group appears on the left side)
    pattern2 = Chem.MolFromSmarts("CCCCCCCCC(=O)O")
    
    # Search for substructure matches using both patterns.
    match1 = mol.GetSubstructMatches(pattern1)
    match2 = mol.GetSubstructMatches(pattern2)
    
    if match1 or match2:
        return True, "Contains a decanoate ester functional group"
    else:
        return False, "No decanoate ester group found in the molecule"
    
# Example usage (for testing purposes):
if __name__ == '__main__':
    # Testing with ethyl decanoate: SMILES: O(C(CCCCCCCCC)=O)CC
    test_smiles = "O(C(CCCCCCCCC)=O)CC"
    result, reason = is_decanoate_ester(test_smiles)
    print(f"Test Molecule: {test_smiles}\nResult: {result}\nReason: {reason}")