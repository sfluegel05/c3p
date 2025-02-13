"""
Classifies: CHEBI:26195 polyphenol
"""
"""
Classifies: CHEBI:16641 polyphenol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polyphenol(smiles: str):
    """
    Determines if a molecule is a polyphenol based on its SMILES string.
    A polyphenol is defined as a member of the class of phenols that contains 2 or more benzene rings, each of which is substituted by at least one hydroxyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyphenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns
    aromatic_ring_pattern = Chem.MolFromSmarts("c1ccccc1")
    hydroxyl_group_pattern = Chem.MolFromSmarts("[OX2H]")
    
    # Find all aromatic rings
    aromatic_rings = mol.GetSubstructMatches(aromatic_ring_pattern)
    
    # Check if there are at least 2 aromatic rings
    if len(aromatic_rings) < 2:
        return False, "Fewer than 2 aromatic rings"
    
    # Check if each aromatic ring has at least one hydroxyl group
    for ring in aromatic_rings:
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        ring_has_hydroxy = any(mol.HasSubstructMatch(hydroxyl_group_pattern, atoms=ring_atoms))
        if not ring_has_hydroxy:
            return False, "At least one aromatic ring lacks a hydroxyl group"
    
    # Check for at least one hydroxyl group outside of the aromatic rings
    non_ring_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIdx() not in [idx for ring in aromatic_rings for idx in ring]]
    has_non_ring_hydroxy = mol.HasSubstructMatch(hydroxyl_group_pattern, atoms=non_ring_atoms)
    
    if has_non_ring_hydroxy:
        return True, "Contains 2 or more aromatic rings, each with at least one hydroxyl group, and at least one additional hydroxyl group"
    else:
        return False, "No hydroxyl groups outside of the aromatic rings"