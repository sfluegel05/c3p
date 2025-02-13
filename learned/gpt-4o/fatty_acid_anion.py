"""
Classifies: CHEBI:28868 fatty acid anion
"""
from rdkit import Chem

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    A fatty acid anion is characterized by a deprotonated carboxylic acid group 
    and a typically long hydrocarbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid anion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylate group pattern [CX3](=O)[O-]
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=O)[O-]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate group found (required for fatty acid anion)"
    
    # Check for long hydrocarbon chain (at least 8 continuous carbon atoms)
    carbon_chain_pattern = Chem.MolFromSmarts("[C][C][C][C][C][C][C][C]")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "No long hydrocarbon chain found, not indicative of a fatty acid anion"

    return True, "Contains a carboxylate group and a long hydrocarbon chain characteristic of a fatty acid anion"