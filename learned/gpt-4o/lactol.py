"""
Classifies: CHEBI:38131 lactol
"""
from rdkit import Chem

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol based on its SMILES string.
    A lactol is a cyclic hemiacetal formed through intramolecular addition
    of a hydroxy group to a carbonyl group.

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

    # Define SMARTS pattern for lactol: cyclic structure with hemiacetal in a ring
    lactol_pattern = Chem.MolFromSmarts("C1OC([OH])C(C)O1")
    if not mol.HasSubstructMatch(lactol_pattern):
        return False, "No cyclic hemiacetal (lactol) pattern found"

    return True, "Contains cyclic hemiacetal (lactol) structure"

# This function takes a SMILES string as input and returns True with a reason if it matches 
# the lactol structure, or False if it does not.