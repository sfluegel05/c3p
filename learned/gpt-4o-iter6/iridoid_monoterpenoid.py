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

    # Define a comprehensive SMARTS pattern for iridoid monoterpenoid structures
    # A cyclopentane with a fused oxygen-containing six-membered ring
    # This considers potential exocyclic double bonds and functionality adjustments
    iridoid_pattern = Chem.MolFromSmarts("C1CC[C@@]2(C1)OC=C2")  # Core cyclopentane and oxygen-heterocycle

    # Check for iridoid core along with functional groups common in iridoid monoterpenoids
    if mol.HasSubstructMatch(iridoid_pattern):
        # Check for significant functional motifs that might aid in identification
        known_functional_groups = [
            Chem.MolFromSmarts("[OH]"),  # Hydroxyl group
            Chem.MolFromSmarts("[O=C]"),  # Ketone group
            Chem.MolFromSmarts("[O=CO]"), # Ester linkage typically seen
        ]
        for fg in known_functional_groups:
            if mol.HasSubstructMatch(fg):
                return True, "Contains iridoid monoterpenoid core with characteristic functional groups"

    return False, "Does not match the typical iridoid monoterpenoid core structure"


# Example usage
smiles_example = 'OCC1=C(C=O)[C@H](C[C@@H]1OC(=O)\\C=C/c1ccc(O)cc1)C(=C)C=O'  # Example SMILES of luzonidial B
result, reason = is_iridoid_monoterpenoid(smiles_example)
print(f"Result: {result}, Reason: {reason}")