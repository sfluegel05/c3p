"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
from rdkit import Chem

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid based on its SMILES string.
    This includes core characteristics such as a cyclopentane ring fused to a six-membered oxygen heterocycle.

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
    
    # Defining core iridoid structure with a cyclopentane ring fused to a six-membered oxygen-containing ring
    cyclic_pattern = Chem.MolFromSmarts("C1C=CO[C@H]2[C@H]1OC=C2")
    
    # Additional functionalities seen in iridoids
    functional_patterns = [
        Chem.MolFromSmarts("C=O"),   # Aldehyde/Carbonyl
        Chem.MolFromSmarts("O=C(O)"), # Carboxylic acid/ester
        Chem.MolFromSmarts("O[C@H]"),  # Hydroxyl group
    ]

    # Verify the core heterocyclic structure
    if mol.HasSubstructMatch(cyclic_pattern):
        # Check for presence of at least one functional group
        for func_pattern in functional_patterns:
            if mol.HasSubstructMatch(func_pattern):
                return True, "Contains cyclopentane and oxygenated ring with functional groups typical of iridoid monoterpenoids"
        return True, "Contains the core cyclopentane and oxygenated ring structure typical of iridoids"
    
    return False, "Does not match the typical iridoid monoterpenoid structure"

# Example usage
smiles_example = 'OCC1=C(C=O)[C@H](C[C@@H]1OC(=O)\\C=C/c1ccc(O)cc1)'  # Example SMILES of luzonidial B
result, reason = is_iridoid_monoterpenoid(smiles_example)
print(f"Result: {result}, Reason: {reason}")