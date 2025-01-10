"""
Classifies: CHEBI:26377 pterocarpans
"""
from rdkit import Chem

def is_pterocarpans(smiles: str):
    """
    Determines if a molecule is a pterocarpan based on its SMILES string.
    Pterocarpans are characterized by a 6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene skeleton
    with additional potential substitutions.

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
    # Adjust the SMARTS pattern to accommodate the 6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene structure
    pterocarpan_pattern = Chem.MolFromSmarts("c1cc2oc3c(ccc4c3Cc3ccc(O)cc3O4)c(c1)CC2")
    if pterocarpan_pattern is None:
        return None, "Error in defining SMARTS pattern"

    # Check if the core structure is present
    if not mol.HasSubstructMatch(pterocarpan_pattern):
        return False, "No pterocarpan core found"

    # Additional checks for common substituents like hydroxyl and methoxy groups
    # These checks aren't mandatory for the core but can enhance identification
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    methoxy_pattern = Chem.MolFromSmarts("CO")
    
    has_hydroxyl = mol.HasSubstructMatch(hydroxyl_pattern)
    has_methoxy = mol.HasSubstructMatch(methoxy_pattern)

    if has_hydroxyl or has_methoxy:
        return True, f"Contains pterocarpan core with {'hydroxyl and methoxy groups' if (has_hydroxyl and has_methoxy) else 'hydroxyl group' if has_hydroxyl else 'methoxy group'}"

    return True, "Contains pterocarpan core structure"