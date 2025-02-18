"""
Classifies: CHEBI:23849 diterpenoid
"""
"""
Classifies: CHEBI:23849 diterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES string.
    A diterpenoid is derived from a diterpene, typically with a C20 skeleton, which may be modified.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Diterpenoids typically have around 20 carbons, but can vary due to modifications
    if c_count < 15 or c_count > 25:
        return False, f"Carbon count ({c_count}) is outside the typical range for diterpenoids"

    # Check for the presence of multiple rings or long carbon chains
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 1:
        return False, "No rings found, which is uncommon for diterpenoids"

    # Check for functional groups commonly found in diterpenoids
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    
    if not (len(ester_matches) > 0 or len(hydroxyl_matches) > 0 or len(carbonyl_matches) > 0):
        return False, "No common functional groups (ester, hydroxyl, carbonyl) found"

    # Check for long carbon chains or complex ring systems
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 3:
        return False, "Too few rotatable bonds for a typical diterpenoid"

    # Check molecular weight - diterpenoids typically have a molecular weight between 250 and 500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 500:
        return False, f"Molecular weight ({mol_wt:.2f} Da) is outside the typical range for diterpenoids"

    return True, "Molecule has characteristics consistent with a diterpenoid (C20 skeleton, functional groups, rings, and rotatable bonds)"