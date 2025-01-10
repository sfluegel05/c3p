"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
from rdkit import Chem

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid based on its SMILES string.
    Incorporates a diverse set of structural motifs typically found in iridoid monoterpenoids.

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
    
    # Define a more comprehensive set of structural patterns for iridoid monoterpenoids
    core_patterns = [
        Chem.MolFromSmarts("C1=CC2OC(O)C=C1OC2"),  # Simplified iridoid-like structure with heterocyclic nature
        Chem.MolFromSmarts("C1=CO[C@@H]2[C@@H]2C=C[C@H]1"),  # More complex bicyclic structures
        Chem.MolFromSmarts("O1C=CC2C(O)C1C2"),  # Considering secoridoid breakages
    ]
    
    # Include patterns for common functional groups
    functional_groups = [
        Chem.MolFromSmarts("[CH]=O"),  # Aldehyde group
        Chem.MolFromSmarts("[C@H](O)[C]=O"),  # Secondary alcohols with adjacent carbonyls
        Chem.MolFromSmarts("C(=O)OC"),  # Ester groups
        Chem.MolFromSmarts("C=C(O)C=O"),  # Enols with carbonyls common in breakdown products
    ]
    
    # Check for core patterns
    for pattern in core_patterns:
        if mol.HasSubstructMatch(pattern):
            # Check for additional identification via functional groups
            for fg in functional_groups:
                if mol.HasSubstructMatch(fg):
                    return True, "Contains iridoid monoterpenoid core with indicative functional groups"
            return True, "Contains iridoid monoterpenoid core"

    return False, "Does not match the typical iridoid monoterpenoid core structure"

# Example usage
smiles_example = 'OCC1=C(C=O)[C@H](C[C@@H]1OC(=O)\\C=C/c1ccc(O)cc1)'  # Example SMILES of luzonidial B
result, reason = is_iridoid_monoterpenoid(smiles_example)
print(f"Result: {result}, Reason: {reason}")