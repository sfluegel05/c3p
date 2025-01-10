"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
from rdkit import Chem

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine based on its SMILES string.
    
    A 1-O-acylglycerophosphoethanolamine is defined as a glycerophosphoethanolamine having 
    an unspecified O-acyl substituent at the 1-position of the glycerol fragment.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-O-acylglycerophosphoethanolamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Refined check for glycerol backbone with O-acyl substituent at 1-position
    # Allows for flexibility in the acyl chain length and configuration
    glycerol_acyl_pattern = Chem.MolFromSmarts("O[C@@H](COC(=O)C)CO")
    if not mol.HasSubstructMatch(glycerol_acyl_pattern):
        return False, "No glycerol backbone with O-acyl substituent at the 1-position found"
    
    # Refined check for phosphoethanolamine group
    # This pattern accounts for the presence of different ionic states and configurations
    phosphoethanolamine_pattern = Chem.MolFromSmarts("COP(=O)(O)OCCN")
    if not mol.HasSubstructMatch(phosphoethanolamine_pattern):
        return False, "No phosphoethanolamine group found"
    
    return True, "Contains glycerol backbone with O-acyl substituent at the 1-position and phosphoethanolamine group"

# Example SMILES strings of 1-O-acylglycerophosphoethanolamines can be tested using the above function.