"""
Classifies: CHEBI:28868 fatty acid anion
"""
from rdkit import Chem

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    A fatty acid anion is characterized by a deprotonated carboxylic acid group (-C([O-])=O)
    and a hydrocarbon chain, which may include double bonds, rings, or functional groups,
    and should have sufficient chain length.

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

    # Define carboxylate SMARTS pattern
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    if carboxylate_pattern is None:
        return (None, "Error creating carboxylate pattern")

    # Detect the presence of a carboxylate group
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate group (-C([O-])=O) found"

    # Count carbon atoms sufficient to indicate a fatty acid type structure
    carbon_pattern = Chem.MolFromSmarts("[CH,CH2,CH3]")
    carbon_matches = mol.GetSubstructMatches(carbon_pattern)
    
    if len(carbon_matches) < 6:
        return False, "Not enough carbons for a typical fatty acid chain"

    return True, "Contains a deprotonated carboxylate group with a sufficiently long hydrocarbon chain"