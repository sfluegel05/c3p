"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
from rdkit import Chem

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine based on its SMILES string.
    A 1-O-acylglycerophosphoethanolamine has a glycerophosphoethanolamine structure with an O-acyl substituent 
    at the 1-position of the glycerol fragment.

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

    # Define enhanced pattern for glycerol backbone with O-acyl group at 1-position
    # It uses chiral [C@@H] to ensure stereochemistry, as needed.
    o_acyl_pattern = Chem.MolFromSmarts("[C@@H](O)COC(=O)[C]")
    if not mol.HasSubstructMatch(o_acyl_pattern):
        return False, "No O-acyl ester correctly positioned at the 1-position"

    # Enhanced pattern for the phosphoethanolamine group
    # Ensure connectivity is strictly matched.
    phosphoethanolamine_pattern = Chem.MolFromSmarts("P(=O)(O)OC[C@H](CN)")
    if not mol.HasSubstructMatch(phosphoethanolamine_pattern):
        return False, "No correctly oriented phosphoethanolamine group present"

    return True, "Contains glycerophosphoethanolamine structure with O-acyl substituent at the 1-position"