"""
Classifies: CHEBI:38131 lactol
"""
"""
Classifies: CHEBI:51475 lactol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol based on its SMILES string.
    A lactol is a cyclic hemiacetal formed by intramolecular addition of a hydroxy group
    to an aldehydic or ketonic carbonyl group.

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
    
    # Look for intramolecular hemiacetal pattern formed by addition of OH to carbonyl
    hemiacetal_pattern = Chem.MolFromSmarts("[OX2r3][CX4r3][CX3](=[OX1])")
    hemiacetal_match = mol.GetSubstructMatches(hemiacetal_pattern)
    
    if not hemiacetal_match:
        return False, "No intramolecular hemiacetal group found"
    
    # Check if the carbonyl is an aldehyde or ketone
    carbonyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[C,O]")
    carbonyl_match = mol.GetSubstructMatches(carbonyl_pattern)
    
    if not carbonyl_match:
        return False, "Carbonyl group is not an aldehyde or ketone"
    
    # Count rotatable bonds - lactols typically have few rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable > 6:
        return False, "Too many rotatable bonds for a cyclic lactol"
    
    # Check for ring size - lactols are typically small rings
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(x) for x in ring_info.AtomRings()]
    
    if not any(size >= 4 and size <= 8 for size in ring_sizes):
        return False, "Ring size not typical for lactols"
    
    return True, "Contains an intramolecular hemiacetal group formed by addition of a hydroxyl group to an aldehydic or ketonic carbonyl group"