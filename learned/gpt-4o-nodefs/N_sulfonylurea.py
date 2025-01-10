"""
Classifies: CHEBI:76983 N-sulfonylurea
"""
from rdkit import Chem

def is_N_sulfonylurea(smiles: str):
    """
    Determines if a molecule is an N-sulfonylurea based on its SMILES string.
    An N-sulfonylurea contains a sulfonyl group (-SO2-) attached to a nitrogen atom,
    and a urea moiety (-NH-C(=O)-NH-) in its structure, specifically structured.

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

    # Adjusted Pattern for N-sulfonylurea: To account for variations in linkage
    # Contains the wider pattern: may include substitutions on the nitrogen or alternative linkages
    n_sulfonylurea_pattern = Chem.MolFromSmarts("[NX3]S(=O)(=O)N")
    urea_pattern = Chem.MolFromSmarts("NC(=O)N")

    # Match both the sulfonyl binding to nitrogen and the urea group presence
    if mol.HasSubstructMatch(n_sulfonylurea_pattern) and mol.HasSubstructMatch(urea_pattern):
        return True, "Contains N-sulfonyl group and urea moiety"

    return False, "Missing specific N-sulfonylurea linkage in the structure"