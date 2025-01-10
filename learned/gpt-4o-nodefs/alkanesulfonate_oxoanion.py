"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
from rdkit import Chem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is an alkanesulfonate oxoanion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define sulfonate group SMARTS
    sulfonate_pattern = Chem.MolFromSmarts("S([O-])(=O)=O")
    if not mol.HasSubstructMatch(sulfonate_pattern):
        return False, "Missing sulfonate group S([O-])(=O)=O"

    # Check for aliphatic carbon connection
    # Typically, an alkanesulfonate should have a connection to an aliphatic carbon chain
    aliphatic_carbon = Chem.MolFromSmiles("C")
    for match in mol.GetSubstructMatches(sulfonate_pattern):
        sulfonate_atom = mol.GetAtomWithIdx(match[0])  # Sulfur atom index
        for neighbor in sulfonate_atom.GetNeighbors():
            if neighbor.GetSymbol() == 'C' and mol.HasSubstructMatch(aliphatic_carbon, atoms=[neighbor.GetIdx()]):
                return True, "Sulfonate group attached to an aliphatic structure"
    
    return False, "Sulfonate group not connected to an aliphatic structure"