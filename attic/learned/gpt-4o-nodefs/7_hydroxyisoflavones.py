"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
from rdkit import Chem

def is_7_hydroxyisoflavones(smiles: str):
    """
    Determines if a molecule is a 7-hydroxyisoflavone based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 7-hydroxyisoflavone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the 7-hydroxyisoflavone core structure
    # Here we attempt to generalize more based on provided examples
    general_pattern = Chem.MolFromSmarts("Oc1ccc2c(c1)occ(-c1ccc(O)c(O)c1)c2=O")

    # Check if the molecule matches the generic 7-hydroxyisoflavone pattern
    if not mol.HasSubstructMatch(general_pattern):
        # If it doesn't match, check a more relaxed core to improve sensitivity, given the diversity
        alternate_pattern = Chem.MolFromSmarts("Oc1ccc2c(c1)occ(-c1ccccc1)c2=O")  # Simplified allowances

        if not mol.HasSubstructMatch(alternate_pattern):
            return False, "Does not match 7-hydroxyisoflavone core structure"

    # Perform additional feature checks if needed
    # For example, checking for additional functional groups characteristic of known examples
    phenyl_ring_check = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(phenyl_ring_check):
        return False, "Missing characteristic phenyl ring structure found in 7-hydroxyisoflavones"
    
    # Further customization: Add specific substituents or configurations that strengthen the presence of a derivative

    return True, "Matches the 7-hydroxyisoflavone core structure"

# Example usage
smiles_examples = [
    "Oc1cc(O)c2c(c1)occ(-c1ccc3OCOc3c1)c2=O",  # Example 7-hydroxyisoflavone
    "COc1ccc(ccc1O)-c1coc2cc(O)ccc2c1=O",  # Example
    # Add more SMILES strings for testing
]

for smiles in smiles_examples:
    result, reason = is_7_hydroxyisoflavones(smiles)
    print(f"SMILES: {smiles} -> Is 7-hydroxyisoflavone? {result}. Reason: {reason}")