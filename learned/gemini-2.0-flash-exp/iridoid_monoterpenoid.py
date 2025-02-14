"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
from rdkit import Chem

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple(bool, str): (True, reason) if molecule is an iridoid monoterpenoid,
                         (False, reason) otherwise.
                         (None, None) if error.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # General Iridoid core structure (cyclopentane fused with 6-membered ring with oxygen)
    # Allowing for various substitutions (X) and double bonds (=) and any bond ~.
    # This pattern is for iridoid.
    iridoid_core_pattern = Chem.MolFromSmarts("[CX4]1~[CX4]~[CX4]~[OX2]~[CX4]2~[CX4]1~[CX4]~[CX4]2")

    # Secoiridoid core structure (cleaved bond in the cyclopentane ring)
    # Defining a pattern with explicit broken bond using disconnected parts.
    secoiridoid_core_pattern = Chem.MolFromSmarts("[CX4]1~[CX4]~[CX4]~[OX2]~[CX4]2~[CX4]~[CX4]2.[CX4]1")
    

    if not (mol.HasSubstructMatch(iridoid_core_pattern) or
            mol.HasSubstructMatch(secoiridoid_core_pattern)):
        return False, "No core iridoid or secoiridoid structure found"
    
    return True, "Iridoid or secoiridoid core structure found"