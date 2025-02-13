"""
Classifies: CHEBI:32863 secondary amine
"""
from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule is a secondary amine based on its SMILES string.
    A secondary amine has a nitrogen atom connected to two carbon atoms and not
    engaged in amide, imide, or other specific exclusionary groups or rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a pattern for secondary amine
    sec_amine_pattern = Chem.MolFromSmarts("[NX3;R0;!$([NX3][CX3](=[OX1])C)]([C&R0])[C&R0]")  # Secondary amine, not in ring, not amide

    # Check for the pattern in the molecule
    if mol.HasSubstructMatch(sec_amine_pattern):
        return True, "Nitrogen is bonded to two carbon atoms and not part of amide or rings, indicative of a secondary amine"

    return False, "No nitrogen atom bonded to exactly two non-cyclic carbons in a secondary amine structure found"