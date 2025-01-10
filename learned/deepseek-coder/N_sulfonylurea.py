"""
Classifies: CHEBI:76983 N-sulfonylurea
"""
"""
Classifies: N-sulfonylurea
Definition: A urea in which one of the hydrogens attached to a nitrogen of the urea group is replaced by a sulfonyl group.
"""
from rdkit import Chem

def is_N_sulfonylurea(smiles: str):
    """
    Determines if a molecule is an N-sulfonylurea based on its SMILES string.
    An N-sulfonylurea has a urea group (-NC(=O)N-) with one of the nitrogen hydrogens replaced by a sulfonyl group (-S(=O)(=O)-).

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

    # Define the N-sulfonylurea pattern: urea group with one nitrogen attached to a sulfonyl group
    n_sulfonylurea_pattern = Chem.MolFromSmarts("[NX3][C](=[OX1])[NX3][S](=[OX1])(=[OX1])")
    
    # Check if the molecule matches the N-sulfonylurea pattern
    if mol.HasSubstructMatch(n_sulfonylurea_pattern):
        return True, "Contains a urea group with one nitrogen attached to a sulfonyl group"
    else:
        return False, "No N-sulfonylurea pattern found"