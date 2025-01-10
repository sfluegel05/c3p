"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
from rdkit import Chem

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 1-O-acylglycerophosphoethanolamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define pattern for 1-O-acylglycerol with a sn-glycero configuration and an acyl group at position 1
    glycerol_pattern = Chem.MolFromSmarts("[C@@H](O)COC(=O)C")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No sn-1-O-acylglycerol pattern found"

    # Define pattern for phosphoethanolamine group
    phosphoethanolamine_pattern = Chem.MolFromSmarts("COP(=O)(O)OCCN")
    if not mol.HasSubstructMatch(phosphoethanolamine_pattern):
        return False, "Phosphoethanolamine group not found"

    # Confirm the long chain of acyl group, at least 15 carbons excluding glycerol backbone
    carbon_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCCC=O")  # Pattern for a long carbon chain
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "No long acyl chain found"

    return True, "Contains 1-O-acylglycerophosphoethanolamine structure with appropriate groups"