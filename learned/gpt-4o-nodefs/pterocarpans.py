"""
Classifies: CHEBI:26377 pterocarpans
"""
from rdkit import Chem

def is_pterocarpans(smiles: str):
    """
    Determines if a molecule is a pterocarpan based on its SMILES string.
    Pterocarpans have a core structure often described as 6a,11a-dihydrobenzofuro[3,2-c]chromene.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool, str: True and reason if the molecule is a pterocarpan, False and reason otherwise
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a core structure resembling pterocarpans
    # This needs to capture the polycyclic nature including furan and chromene types
    pterocarpan_pattern = Chem.MolFromSmarts("c1c2c(O3)ccc(O3)c2c4c1OCC4")
    
    if pterocarpan_pattern is None:
        # If the pattern fails to compile, return a diagnostic message
        return (None, "Error in constructing SMARTS pattern")

    # Check if molecule contains the pterocarpan core structure
    if mol.HasSubstructMatch(pterocarpan_pattern):
        return True, "Contains the pterocarpan core structure"
    
    return False, "Does not contain the pterocarpan core structure"