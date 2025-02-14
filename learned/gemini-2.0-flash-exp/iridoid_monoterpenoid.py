"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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
    
    # More specific SMARTS for iridoid core
    iridoid_core_pattern = Chem.MolFromSmarts("[CH2]1[CH]=[CH][CH2][C]2([O][CH2]1)[CH]=[CH][CH2]2") # Core ring with typical double bonds
    iridoid_core_pattern2 = Chem.MolFromSmarts("[CH2]1[CH][CH][CH2][C]2([O][CH2]1)[CH]=[CH][CH2]2") # Core ring with one double bonds
    iridoid_core_pattern3 = Chem.MolFromSmarts("[CH2]1[CH]=[CH][CH2][C]2([O][CH2]1)[CH][CH][CH2]2")  # Core ring with one double bonds
    iridoid_core_pattern4 = Chem.MolFromSmarts("[CH2]1[CH][CH][CH2][C]2([O][CH2]1)[CH][CH][CH2]2") # Core ring without double bonds
    # Secoiridoid core structure (cleaved bond in the cyclopentane ring)
    secoiridoid_core_pattern = Chem.MolFromSmarts("[CX4]1[CX4]~[CX4]~[OX2]~[CX4]2[CX4]([CX4]1)~[CX4]~[CX4]2")

    if not (mol.HasSubstructMatch(iridoid_core_pattern) or
            mol.HasSubstructMatch(iridoid_core_pattern2) or
             mol.HasSubstructMatch(iridoid_core_pattern3) or
             mol.HasSubstructMatch(iridoid_core_pattern4) or
            mol.HasSubstructMatch(secoiridoid_core_pattern)):
        return False, "No core iridoid or secoiridoid structure found"

    # Check for 10 carbon atoms, as a monoterpenoid.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
        return False, f"Fewer than 10 carbon atoms found ({c_count})"
    
    return True, "Iridoid or secoiridoid core structure found with at least 10 carbons"