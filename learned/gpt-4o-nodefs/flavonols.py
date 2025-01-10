"""
Classifies: CHEBI:28802 flavonols
"""
from rdkit import Chem

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol based on its SMILES string.
    Flavonols are characterized by the presence of a 3-hydroxyflavone backbone
    with varying substituents.

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
    
    # Define a more flexible SMARTS pattern for a flavonol core structure
    # Flavonol core: 3-hydroxyflavone
    flavonol_pattern = Chem.MolFromSmarts("Oc1cc2c(cc1)-c1oc(=O)cc(c1O)c2")  # core flavonol structure derived
    
    # Check if the molecule has a flavonol core structure
    if not mol.HasSubstructMatch(flavonol_pattern):
        return False, "Flavonol core structure not found"

    # Check for typical substituents, often hydroxyl or methoxy groups attached to the flavonoid rings
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    methoxy_pattern = Chem.MolFromSmarts("CO")
    # Reason for possible substituents count
    substituents_count = len(mol.GetSubstructMatches(hydroxyl_pattern)) + len(mol.GetSubstructMatches(methoxy_pattern))
    # Many flavonols have multiple such substituents
    if substituents_count < 3:  # Assume at least 3 for typical flavonol complexity
        return False, f"Insufficient typical substituents for flavonol classification, found {substituents_count}"
    
    return True, "Contains a flavonol core structure with expected substituents"

# Example of usage:
# result, reason = is_flavonols("COc1cc(O)c2c(c1)oc(-c1ccc(O)cc1)c(O)c2=O")
# print(result, reason)  # Expecting: True, "Contains a flavonol core structure with expected substituents"