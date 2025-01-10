"""
Classifies: CHEBI:13248 anilide
"""
from rdkit import Chem

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.
    An anilide is defined as an aromatic amide obtained by acylation of aniline.

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

    # Define aniline substructure - phenyl group bound to a nitrogen
    aniline_pattern = Chem.MolFromSmarts("c1ccccc1N")

    # Define acyl group substructure (C=O attached to nitrogen)
    acyl_pattern = Chem.MolFromSmarts("NC(=O)")

    # Check for aniline substructure
    if not mol.HasSubstructMatch(aniline_pattern):
        return False, "No aniline moiety found (phenyl ring with nitrogen)"

    # Check for acyl group substructure
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No acyl (C=O) group attached to nitrogen"

    return True, "Contains aniline moiety with acyl group, fitting anilide definition"