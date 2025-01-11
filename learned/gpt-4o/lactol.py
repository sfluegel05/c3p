"""
Classifies: CHEBI:38131 lactol
"""
from rdkit import Chem

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol based on its SMILES string.
    A lactol is a cyclic hemiacetal formed through intramolecular addition
    of a hydroxy group to a carbonyl group, forming a cyclic ether.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lactol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define improved SMARTS pattern for lactol
    lactol_patterns = [
        Chem.MolFromSmarts("[C!H0]1[O][C!H0][C]([OH])[C][O]1"),  # Generic cyclic hemiacetal
        Chem.MolFromSmarts("[C;R](O)([O;R])[C;R]")               # Cyclic structure with a hydroxyl and an ether linkage
    ]

    for pattern in lactol_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains cyclic hemiacetal (lactol) structure"

    return False, "No cyclic hemiacetal (lactol) pattern found"

# This function takes a SMILES string as input and returns True with a reason if it matches
# the lactol structure, or False if it does not.