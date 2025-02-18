"""
Classifies: CHEBI:59644 oxo fatty acid
"""
"""
Classifies: oxo fatty acid (CHEBI:76224)
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmiles, MolFromSmarts

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.
    An oxo fatty acid contains at least one aldehydic or ketonic group in addition to the carboxylic acid group.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an oxo fatty acid, False otherwise
        str: Reason for classification
    """
    mol = MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES"
    
    # Check for at least one carboxylic acid group (-COOH)
    carboxylic_acid = MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group"
    
    # Find all carbonyl groups (keto or aldehyde), excluding the carboxylic acid's
    # Aldehyde: CH=O (but not part of COOH)
    aldehyde = MolFromSmarts("[CX3H1](=O)[!#8]")
    # Ketone: C(=O)C where not part of COOH
    ketone = MolFromSmarts("[CX3](=O)[CX4]")
    
    # Count matches for aldehyde/ketone not in carboxylic acid
    oxo_matches = len(mol.GetSubstructMatches(aldehyde)) + len(mol.GetSubstructMatches(ketone))
    if oxo_matches < 1:
        return False, "No additional oxo groups"
    
    # Basic check for sufficient chain length (at least 8 carbons in main chain)
    # This is a heuristic and might need adjustment
    carbon_chain = MolFromSmarts("[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]")
    if not mol.HasSubstructMatch(carbon_chain):
        return False, "Insufficient chain length for fatty acid"
    
    return True, "Contains carboxylic acid with additional oxo group(s)"