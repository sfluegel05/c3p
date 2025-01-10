"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
from rdkit import Chem

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an iridoid monoterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the iridoid core (cyclopentane fused to a six-membered oxygen heterocycle)
    # This pattern is an approximation due to the complexity of variations
    iridoid_pattern = Chem.MolFromSmarts("C1CCC2OCCCC2C1")  # Simplified core structure pattern

    # Check for the iridoid pattern in the molecule
    if mol.HasSubstructMatch(iridoid_pattern):
        return True, "Contains iridoid monoterpenoid core structure"
    else:
        return False, "Does not match the typical iridoid monoterpenoid core structure"

# Example usage
smiles_example = 'OCC1=C(C=O)[C@H](C[C@@H]1OC(=O)\C=C/c1ccc(O)cc1)C(=C)C=O'  # Example SMILES of luzonidial B
result, reason = is_iridoid_monoterpenoid(smiles_example)
print(f"Result: {result}, Reason: {reason}")