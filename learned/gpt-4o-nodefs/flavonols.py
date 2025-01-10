"""
Classifies: CHEBI:28802 flavonols
"""
from rdkit import Chem

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol based on its SMILES string.
    Flavonols have a 3-hydroxyflavone backbone: a 15-carbon skeleton with a ketone group on the 
    heterocyclic C-ring and a hydroxyl group at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for a flavonol core structure
    flavonol_pattern = Chem.MolFromSmarts("c1cc(O)cc2c(c1)oc(c(=O)c2O)-c1ccccc1")
    
    # Check if the molecule has a flavonol core structure
    if not mol.HasSubstructMatch(flavonol_pattern):
        return False, "Flavonol core structure not found"
    
    return True, "Contains a flavonol core structure with the required functional groups"

# Example of usage:
# result, reason = is_flavonols("COc1cc(O)c2c(c1)oc(-c1ccc(O)cc1)c(O)c2=O") 
# print(result, reason)  # Expects to print: True, "Contains a flavonol core structure with the required functional groups"