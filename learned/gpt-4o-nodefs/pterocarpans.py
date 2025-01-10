"""
Classifies: CHEBI:26377 pterocarpans
"""
from rdkit import Chem

def is_pterocarpans(smiles: str):
    """
    Determines if a molecule is a pterocarpan based on its SMILES string.
    Pterocarpans have a polycyclic structure often described with a
    benzofurochromene-like core.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool, str: True and reason if the molecule is a pterocarpan, False and reason otherwise
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Attempt a more inclusive pattern for pterocarpans
    # This pattern needs to be validated against diverse examples
    # of pterocarpans noted for their unique aromatic fused rings.
    pterocarpan_pattern = Chem.MolFromSmarts("C1Oc2ccccc2C3=C(O1)C=CC=C3")  # Hypothetical improved pattern
    
    if pterocarpan_pattern is None:
        return (None, "Error in constructing SMARTS pattern")

    # Check for presence of the core pterocarpan structure
    if mol.HasSubstructMatch(pterocarpan_pattern):
        return True, "Contains the pterocarpan core structure"
    
    return False, "Does not contain the pterocarpan core structure"