"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
from rdkit import Chem

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid based on its SMILES string.
    Characterized by a bicyclic system: a cyclopentane ring fused with a six-membered oxygen heterocycle.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an iridoid monoterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Core iridoid pattern with fused cyclopentane and six-membered ether ring
    core_pattern = Chem.MolFromSmarts("[C@H]1C[C@@H]2COC=C[C@@H]2C(=O)O1")
    
    # Seco-iridoid pattern where the cyclopentane may be cleaved, retaining key oxygens
    seco_pattern = Chem.MolFromSmarts("[C@@H]1COC=C[C@H]2O[C@H]2C(=O)[C@H]1")
    
    # Additional functional patterns often seen in iridoids
    functional_groups = [
        Chem.MolFromSmarts("C=O"),  # Carbonyl/aldehyde
        Chem.MolFromSmarts("O[C@H]"),  # Hydroxyl groups
        Chem.MolFromSmarts("OC(=O)"),  # Ester linkages
    ]

    # Check for core or seco-iridoid structures
    if mol.HasSubstructMatch(core_pattern) or mol.HasSubstructMatch(seco_pattern):
        for func_group in functional_groups:
            if mol.HasSubstructMatch(func_group):
                return True, "Contains iridoid core structure along with associated functional groups"
        return True, "Contains the core iridoid structure but lacks identifiable functional groups"
    
    return False, "Does not match the typical iridoid monoterpenoid structure"

# Example usage
smiles_example = 'OCC1=C(C=O)[C@H](C[C@@H]1OC(=O)\\C=C/c1ccc(O)cc1)'  # Example SMILES of luzonidial B
result, reason = is_iridoid_monoterpenoid(smiles_example)
print(f"Result: {result}, Reason: {reason}")