"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
from rdkit import Chem

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid based on its SMILES string.
    A 2-hydroxydicarboxylic acid has one hydroxyl group at the alpha-carbon relative to two
    carboxyl groups. 

    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is a 2-hydroxydicarboxylic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Create SMARTS pattern for 2-hydroxydicarboxylic acid
    # Matches a carbon with hydroxyl and two neighboring carboxylic acid groups.
    pattern = Chem.MolFromSmarts("[#6]([CH](O))[#6](C(=O)O)[#6](C(=O)O)")
    
    # Check if the molecule matches the pattern
    if mol.HasSubstructMatch(pattern):
        return True, "Contains 2-hydroxy group and two carboxylic acid groups in a suitable position"
    
    return False, "Does not have the structural characteristics of a 2-hydroxydicarboxylic acid"