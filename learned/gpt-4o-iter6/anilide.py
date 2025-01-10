"""
Classifies: CHEBI:13248 anilide
"""
from rdkit import Chem

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.
    An anilide is an aromatic amide obtained by acylation of aniline,
    which typically includes a phenyl group directly attached to a nitrogen of an amide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anilide, False otherwise.
        str: Reason for classification.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Enhanced anilide pattern: phenyl group directly bonded to amide nitrogen with potential additional substituents
    # In order to mitigate false positives, we add some complexity to the SMARTS pattern that accounts for common cases
    # But preserves the primary structure of aniline-based acylation.
    anilide_pattern = Chem.MolFromSmarts("c1ccc(cc1)NC(=O)C")

    # Check for anilide pattern
    match = mol.HasSubstructMatch(anilide_pattern)
    
    if not match:
        return False, "No anilide pattern found (phenyl group directly bonded to amide nitrogen)"

    return True, "Contains phenyl group directly bonded to amide nitrogen, fitting anilide definition"