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
        bool: True if molecule is a tetrahydrofuranone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for tetrahydrofuranone
    # This pattern matches a five-membered saturated ring with 4 carbons and 1 oxygen
    # and a carbonyl group in the ring.
    tetrahydrofuranone_pattern = Chem.MolFromSmarts("C1CCC(=O)O1")
    
    # Check if the structure matches the tetrahydrofuranone pattern
    if mol.HasSubstructMatch(tetrahydrofuranone_pattern):
        return True, "Contains a tetrahydrofuranone structure"
    else:
        return False, "Does not contain a tetrahydrofuranone structure"

# Example usage
# Test the function with a sample SMILES string of a known tetrahydrofuranone
test_smiles = "CCCCCCCC1OC(=O)C(=C)C1C(O)=O"
result, reason = is_tetrahydrofuranone(test_smiles)
print(f"Is Tetrahydrofuranone: {result}, Reason: {reason}")