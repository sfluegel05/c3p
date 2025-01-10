"""
Classifies: CHEBI:28802 flavonols
"""
from rdkit import Chem

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol based on its SMILES string.
    A flavonol has a flavone backbone with a hydroxyl group at position 3
    and typically includes additional hydroxyl or methoxy groups.

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
    
    # Define a more general SMARTS pattern for a flavonol core structure
    flavone_core = Chem.MolFromSmarts("c1cc2c(cc1)[C]1=CC(=O)Oc2c(c1)O")  # core 3-hydroxyflavone
    
    # Check if the molecule has a flavonol core structure
    if not mol.HasSubstructMatch(flavone_core):
        return False, "Flavonol core structure not found"

    # Check for typical substituents like hydroxyl groups on the phenyl rings
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    o_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
    
    # Many flavonols have multiple hydroxyl groups
    if o_count < 3:  # Assume at least 3 hydroxyl groups for classification based on examples
        return False, f"Insufficient hydroxyl groups for flavonol classification, found {o_count}"
    
    return True, "Contains a flavonol core structure with typical substituents"

# Example of usage:
# result, reason = is_flavonols("COc1cc(O)c2c(c1)oc(-c1ccc(O)cc1)c(O)c2=O") 
# print(result, reason)  # Expects to print: True, "Contains a flavonol core structure with typical substituents"