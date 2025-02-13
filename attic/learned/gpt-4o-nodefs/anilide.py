"""
Classifies: CHEBI:13248 anilide
"""
from rdkit import Chem

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.
    An anilide is defined by a phenyl group bonded to an amide (-C(=O)N-) linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anilide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a phenyl group connected to an amide linkage
    anilide_pattern = Chem.MolFromSmarts("c1ccccc1NC(=O)")
    if mol.HasSubstructMatch(anilide_pattern):
        return True, "Contains phenyl group connected to an amide linkage"

    return False, "Does not contain the characteristic anilide structure"