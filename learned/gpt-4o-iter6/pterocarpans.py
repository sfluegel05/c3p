"""
Classifies: CHEBI:26377 pterocarpans
"""
from rdkit import Chem

def is_pterocarpans(smiles: str):
    """
    Determines if a molecule is a pterocarpan based on its SMILES string.
    Pterocarpans are characterized by a benzofurochromene skeleton with
    additional potential substructures.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a pterocarpan, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the core structure of pterocarpans
    # Below SMARTS pattern is a generic representation of the benzofurochromene skeleton
    # Considering the 6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene core and variants
    # Example SMARTS pattern for this structure
    pterocarpan_pattern = Chem.MolFromSmarts("c1cc2oc3cccc4occc(c1)c2c34")
    if pterocarpan_pattern is None:
        return None, "Error in defining SMARTS pattern"

    # Check if the core structure is present
    if not mol.HasSubstructMatch(pterocarpan_pattern):
        return False, "No pterocarpan core found"

    # Additional functional group checks, such as hydroxyl, methoxy, etc.
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    methoxy_pattern = Chem.MolFromSmarts("CO")
    
    if mol.HasSubstructMatch(hydroxyl_pattern) and mol.HasSubstructMatch(methoxy_pattern):
        return True, "Contains pterocarpan core with common functional groups"
    
    return True, "Contains pterocarpan core structure"