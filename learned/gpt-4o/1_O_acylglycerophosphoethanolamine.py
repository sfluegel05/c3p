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

    # Generalized pattern for glycerol backbone with O-acyl substituent at the 1-position
    glycerol_acyl_pattern = Chem.MolFromSmarts("O[C@@H](CO)COC(=O)C")
    if not mol.HasSubstructMatch(glycerol_acyl_pattern):
        return False, "No glycerol backbone with O-acyl substituent at the 1-position found"

    # Pattern for phosphoethanolamine group, allowing for various ionic forms
    phosphoethanolamine_pattern = Chem.MolFromSmarts("OP(=O)(O)OCCN")
    if not mol.HasSubstructMatch(phosphoethanolamine_pattern):
        return False, "No phosphoethanolamine group found"
    
    # Ensure no conflicting groups like phosphocholine are present
    phosphocholine_pattern = Chem.MolFromSmarts("OP(=O)(O)OCC[N+](C)(C)C")
    if mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "Contains phosphocholine group, not a phosphoethanolamine"

    return True, "Contains glycerol backbone with O-acyl substituent at the 1-position and phosphoethanolamine group"