"""
Classifies: CHEBI:27325 xanthophyll
"""
"""
Classifies: CHEBI:37208 xanthophyll

A subclass of carotenoids consisting of the oxygenated carotenes.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    A xanthophyll is an oxygenated carotenoid, containing at least one oxygen atom
    and having a specific structural pattern of alternating single and double bonds
    in a long conjugated carbon chain, often with ring structures at one or both ends.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a xanthophyll, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for oxygen atoms
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count == 0:
        return False, "No oxygen atoms found, not a xanthophyll"
    
    # Check for conjugated carbon chain pattern
    carotenoid_pattern = Chem.MolFromSmarts("[cH2]=[c;!r][c;!r]=[c;!r][c;!r]=[c;!r][c;!r]=[c;!r][c;!r]=[c;!r]")
    carotenoid_matches = mol.GetSubstructMatches(carotenoid_pattern)
    if not carotenoid_matches:
        return False, "No carotenoid-like conjugated carbon chain found, not a xanthophyll"
    
    # Check for ring structures at the ends of the carbon chain
    ring_pattern = Chem.MolFromSmarts("[R1]~[R1]~[R1]~[R1]~[R1]~[R1]")
    ring_matches = mol.GetSubstructMatches(ring_pattern)
    
    # Check molecular weight and atom counts
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400 or mol_wt > 1000:
        return False, "Molecular weight outside typical range for xanthophylls"
    
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 20 or carbon_count > 60:
        return False, "Carbon count outside typical range for xanthophylls"
    
    # If all conditions are met, classify as a xanthophyll
    if ring_matches:
        reason = "Contains oxygenated conjugated carbon chain with ring structures, characteristic of xanthophylls"
    else:
        reason = "Contains oxygenated conjugated carbon chain, characteristic of xanthophylls"
    
    return True, reason