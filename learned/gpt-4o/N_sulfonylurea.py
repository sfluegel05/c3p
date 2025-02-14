"""
Classifies: CHEBI:76983 N-sulfonylurea
"""
from rdkit import Chem

def is_N_sulfonylurea(smiles: str):
    """
    Determines if a molecule is an N-sulfonylurea based on its SMILES string.
    N-sulfonylureas are characterized by a sulfonyl group attached to a nitrogen of a urea moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-sulfonylurea, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the N-sulfonylurea pattern
    # - A urea group 'NC(=O)N' where one of the nitrogens is bonded to a sulfonyl group 'S(=O)(=O)'
    nsulfonylurea_pattern = Chem.MolFromSmarts("NC(=O)NS(=O)(=O)")
    
    # Check if the molecule has this pattern
    match = mol.HasSubstructMatch(nsulfonylurea_pattern)
    if match:
        return True, "Contains N-sulfonylurea moiety"
    else:
        return False, "Does not contain N-sulfonylurea moiety"