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

    # Improved SMARTS pattern for anilides: Ensures that a phenyl is directly connected to an amide nitrogen
    # The pattern checks for a phenyl group connected to the N of an amide group which is a common representation of anilide
    anilide_pattern = Chem.MolFromSmarts("[cH]1[cH][cH][cH][cH][cH]1-N-C(=O)")
    
    match = mol.HasSubstructMatch(anilide_pattern)

    if not match:
        return False, "No anilide pattern found (phenyl group directly bonded to amide nitrogen)"

    return True, "Contains phenyl group directly bonded to amide nitrogen, fitting anilide definition"