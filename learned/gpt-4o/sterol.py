"""
Classifies: CHEBI:15889 sterol
"""
from rdkit import Chem

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.
    Sterols typically consist of a 3-hydroxy steroid skeleton closely related to cholestan-3-ol with side chains.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) indicating whether the molecule is a sterol and the reason for classification.
    """

    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for sterane skeleton: ABCD ring system
    sterane_pattern = Chem.MolFromSmarts('C1CCC2(C)C(CC3C2C4CC(CC3C4)C)C1')
    if not mol.HasSubstructMatch(sterane_pattern):
        return False, "No sterane skeleton found"

    # SMARTS pattern to identify a hydroxyl group on a tertiary carbon
    # This should match a carbon atom that is part of a ring system and has exactly one bonded oxygen (hydroxyl)
    hydroxy_smarts = Chem.MolFromSmarts('[C&r3][C&r]([O])[C&r]')
    if not mol.HasSubstructMatch(hydroxy_smarts):
        return False, "No 3-hydroxy group found"

    # Check for additional structural features typical of sterols
    # Allow for branching side chains indicative of sterol complexity
    branch_smarts = Chem.MolFromSmarts('[C&r3]C([C&!r])([C&!r])([C&!r])')
    if not mol.HasSubstructMatch(branch_smarts):
        return False, "No variable side chains or functional group combinations found"

    return True, "Contains sterane skeleton with 3-hydroxy group and potential side chains"

# Testing the function with examples provided can help verify the functionality.