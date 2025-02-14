"""
Classifies: CHEBI:13248 anilide
"""
from rdkit import Chem

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.
    An anilide is an aromatic amide obtained by acylation of aniline.

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

    # Define the SMARTS pattern for the anilide core.
    #   - [c:1] is a carbon in an aromatic ring
    #   - [N:2H] is the nitrogen connected to the ring, with at least one H (aniline nitrogen)
    #   - [C:3](=[O:4]) is the carbonyl carbon of the amide and the O
    anilide_pattern = Chem.MolFromSmarts("[c:1]1[c:1][c:1][c:1][c:1][c:1]1[N:2H][C:3](=[O:4])")

    if mol.HasSubstructMatch(anilide_pattern):
         return True, "Correct Anilide core detected"
    
    return False, "Anilide core not found"