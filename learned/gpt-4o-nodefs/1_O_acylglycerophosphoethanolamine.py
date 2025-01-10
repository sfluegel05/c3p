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
        bool: True if molecule is a 1-O-acylglycerophosphoethanolamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define pattern for 1-O-acylglycerol: a glycerol backbone with an acyl group at position 1
    acyl_glycerol_pattern = Chem.MolFromSmarts("OCC(O)COC(=O)C")
    if not mol.HasSubstructMatch(acyl_glycerol_pattern):
        return False, "No 1-O-acylglycerol pattern found"

    # Define pattern for phosphoethanolamine
    phosphoethanolamine_pattern = Chem.MolFromSmarts("COP(O)(=O)OCCN")
    if not mol.HasSubstructMatch(phosphoethanolamine_pattern):
        return False, "Phosphoethanolamine group not found"
    
    # Check carbon chain length of the linked acyl group (usually a long chain)
    # Count carbon atoms (excluding rings or erroneous interpretations)
    acyl_chain_pattern = Chem.MolFromSmarts("C=O")
    matches = mol.GetSubstructMatches(acyl_chain_pattern)
    if not matches:
        return False, "No acyl chain found"
    
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 15:
        return False, f"Not enough carbon atoms for a typical long acyl chain, found {carbon_count}"

    return True, "Contains 1-O-acylglycerophosphoethanolamine structure with appropriate groups"