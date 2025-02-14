"""
Classifies: CHEBI:76983 N-sulfonylurea
"""
from rdkit import Chem

def is_N_sulfonylurea(smiles: str):
    """
    Determines if a molecule is an N-sulfonylurea based on its SMILES string.
    An N-sulfonylurea is a urea group where one nitrogen is bonded to a sulfonyl group.

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

    # Define the SMARTS pattern for N-sulfonylurea
    n_sulfonylurea_pattern = Chem.MolFromSmarts("[NX2][CX3](=[OX1])[NX2][S](=[OX2])(=[OX2])")


    # Check for the substructure
    if mol.HasSubstructMatch(n_sulfonylurea_pattern):
      return True, "Contains N-sulfonylurea group"
    else:
      return False, "Does not contain N-sulfonylurea group"