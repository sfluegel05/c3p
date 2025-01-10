"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
from rdkit import Chem

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid based on its SMILES string.
    Acknowledging structural diversity, especially of the core cyclopentane-fused to an oxygen-containing heterocycle.

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
    
    # Define potential structural patterns for iridoid monoterpenoids
    core_patterns = [
        Chem.MolFromSmarts("C1CCC2OCCC1O2"),  # Basic Iridoid core
        Chem.MolFromSmarts("C1C=CC2OC=CC1C2"), # Variations with double bond
        Chem.MolFromSmarts("C1CCC2(O)OCC1O2"), # Secoiridoid structures
    ]
    
    # Common functional groups in iridoids
    functional_groups = [
        Chem.MolFromSmarts("[CH]=O"),  # Aldehyde group
        Chem.MolFromSmarts("[C@H]([OH])[CH]=O"), # Hydroxyl paired with aldehydes
        Chem.MolFromSmarts("C[C@H](O)COC=O")  # Esters in secoiridoids
    ]
    
    # Check for core patterns
    for pattern in core_patterns:
        if mol.HasSubstructMatch(pattern):
            # Check for additional identification via functional groups
            for fg in functional_groups:
                if mol.HasSubstructMatch(fg):
                    return True, "Contains iridoid monoterpenoid core with characteristic functional groups"
            return True, "Contains iridoid monoterpenoid core"

    return False, "Does not match the typical iridoid monoterpenoid core structure"

# Example usage
smiles_example = 'OCC1=C(C=O)[C@H](C[C@@H]1OC(=O)\\C=C/c1ccc(O)cc1)'  # Example SMILES of luzonidial B
result, reason = is_iridoid_monoterpenoid(smiles_example)
print(f"Result: {result}, Reason: {reason}")