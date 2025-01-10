"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrahydrofuranone derivative, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for tetrahydrofuranone
    # This pattern matches a five-membered ring with one oxygen (tetrahydrofuran)
    # containing a carbonyl group (=O), but allows for attached groups
    tetrahydrofuranone_pattern1 = Chem.MolFromSmarts("C1COC(=O)C1")
    tetrahydrofuranone_pattern2 = Chem.MolFromSmarts("C1C(=O)OCC1") # for variations

    # Check if the structure matches any of the tetrahydrofuranone patterns
    if mol.HasSubstructMatch(tetrahydrofuranone_pattern1) or mol.HasSubstructMatch(tetrahydrofuranone_pattern2):
        return True, "Contains a tetrahydrofuranone-like structure"
    else:
        return False, "Does not contain a tetrahydrofuranone-like structure"

# Example usage
# Test the function with a sample SMILES string of a known tetrahydrofuranone
test_smiles = "CCCCCCCC1OC(=O)C(=C)C1C(O)=O"
result, reason = is_tetrahydrofuranone(test_smiles)
print(f"Is Tetrahydrofuranone: {result}, Reason: {reason}")