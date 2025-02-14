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

    # Define the SMARTS pattern for the N-sulfonylurea group
    # This pattern represents a urea (N-C(=O)-N) with one N connected to a sulfonyl group (S(=O)(=O)).
    # The ~ means "any bond" and the nitrogen connected to the sulfonyl can have any connections.
    n_sulfonylurea_pattern = Chem.MolFromSmarts("[NX3]([CX3](=[OX1])[NX3])~[S](=[OX2])(=[OX2])")

    # Check if the molecule contains the N-sulfonylurea group
    if mol.HasSubstructMatch(n_sulfonylurea_pattern):
        return True, "Contains N-sulfonylurea group"
    else:
        return False, "Does not contain N-sulfonylurea group"