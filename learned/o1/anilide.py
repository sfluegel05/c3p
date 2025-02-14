"""
Classifies: CHEBI:13248 anilide
"""
"""
Classifies: anilide
"""
from rdkit import Chem

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.
    An anilide is any aromatic amide obtained by acylation of aniline,
    meaning it has an amide group where the nitrogen is directly attached
    to an aromatic ring.

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
    
    # Define anilide pattern: aromatic ring bonded to nitrogen of an amide group
    anilide_smarts = "[c][N;X3][C](=O)"
    anilide_pattern = Chem.MolFromSmarts(anilide_smarts)
    
    if mol.HasSubstructMatch(anilide_pattern):
        return True, "Contains an anilide moiety (aromatic amide from acylation of aniline)"
    else:
        return False, "No anilide moiety found (amide nitrogen not attached to aromatic ring)"