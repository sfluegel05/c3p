"""
Classifies: CHEBI:27325 xanthophyll
"""
"""
Classifies: xanthophylls (oxygenated carotenoids)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    Xanthophylls are oxygenated carotenoids.

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

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Xanthophylls are typically C40 terpenoids
    if c_count < 35 or c_count > 45:
        return False, f"Carbon count {c_count} outside typical range for xanthophylls (35-45)"
    
    # Must have oxygen atoms
    if o_count == 0:
        return False, "No oxygen atoms found - xanthophylls must be oxygenated"
        
    # Look for conjugated polyene pattern (alternating single-double bonds)
    conjugated_pattern = Chem.MolFromSmarts("C=CC=CC=CC=C")
    if not mol.HasSubstructMatch(conjugated_pattern):
        return False, "Missing characteristic conjugated polyene system"
    
    # Count conjugated double bonds
    conjugated_matches = len(mol.GetSubstructMatches(conjugated_pattern))
    if conjugated_matches < 2:
        return False, "Insufficient conjugation for xanthophyll"
    
    # Look for methyl branches (common in carotenoids)
    methyl_pattern = Chem.MolFromSmarts("CC(C)C")
    if not mol.HasSubstructMatch(methyl_pattern):
        return False, "Missing characteristic methyl branching"

    # Check for oxygen-containing functional groups
    hydroxy_pattern = Chem.MolFromSmarts("OH")
    keto_pattern = Chem.MolFromSmarts("C(=O)C")
    epoxy_pattern = Chem.MolFromSmarts("OC(C)C")
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    
    oxygen_groups = (
        len(mol.GetSubstructMatches(hydroxy_pattern)) +
        len(mol.GetSubstructMatches(keto_pattern)) +
        len(mol.GetSubstructMatches(epoxy_pattern)) +
        len(mol.GetSubstructMatches(ester_pattern))
    )
    
    if oxygen_groups == 0:
        return False, "No oxygen-containing functional groups found"

    # Count rings - most xanthophylls have at least one
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count == 0:
        return False, "No rings found - most xanthophylls have cyclic end groups"
        
    # Calculate degree of unsaturation
    double_bonds = rdMolDescriptors.CalcNumDoubleBonds(mol)
    if double_bonds < 8:
        return False, "Insufficient double bonds for xanthophyll conjugated system"

    return True, f"Oxygenated carotenoid with {oxygen_groups} oxygen-containing groups, {double_bonds} double bonds, and {ring_count} rings"