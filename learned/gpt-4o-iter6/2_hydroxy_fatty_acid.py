"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
from rdkit import Chem

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    A 2-hydroxy fatty acid has a hydroxy group in the alpha- or 2-position relative to the carboxylic acid group.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a 2-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for the 2-hydroxy fatty acid: [-][OH]C-[C](=O)O (saturated and unsaturated considered)
    # Match carboxylic acid group (C(=O)O) with hydroxyl on adjacent carbon
    # The use of wildcard neighbors allows for the capture of long chains and potential stereocenters
    pattern = Chem.MolFromSmarts("[CX4,CX3][C@@H1|C@H1;!R](O)[CX3](=O)O")  # Specifying chiral centers should capture R/S forms
    if mol.HasSubstructMatch(pattern):
        return True, "Contains 2-hydroxy group with carboxylic acid at the 2-position"
    
    return False, "No 2-hydroxy group found at the alpha position"

# Example usage
print(is_2_hydroxy_fatty_acid("CCCCCCC[C@@H](O)C(O)=O"))  # (R)-2-hydroxynonanoic acid