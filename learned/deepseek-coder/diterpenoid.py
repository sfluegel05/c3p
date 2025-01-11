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
    A diterpenoid is derived from a diterpene (C20 skeleton) and may have undergone modifications.

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
    
    # Diterpenoids are derived from a C20 skeleton, but modifications may reduce or increase the count
    if c_count < 15 or c_count > 35:
        return False, f"Carbon count ({c_count}) is outside the typical range for diterpenoids (15-35)"

    # Check for terpenoid-like structure (rings and branching)
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings < 1:
        return False, "No rings found, unlikely to be a diterpenoid"

    # Check for typical functional groups in diterpenoids (e.g., alcohols, ketones, esters, ethers)
    functional_groups = ["[OH]", "[C]=O", "[O]", "[C](=O)[O]", "[O][C]", "[C](=O)[C]"]
    has_functional_group = False
    for group in functional_groups:
        pattern = Chem.MolFromSmarts(group)
        if pattern is not None and mol.HasSubstructMatch(pattern):
            has_functional_group = True
            break
    if not has_functional_group:
        return False, "No typical diterpenoid functional groups found"

    # Check molecular weight - diterpenoids typically have MW between 250 and 500
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 500:
        return False, f"Molecular weight ({mol_wt}) is outside the typical range for diterpenoids (250-500)"

    # Check for some degree of branching
    n_branches = sum(1 for atom in mol.GetAtoms() if len(atom.GetNeighbors()) > 2)
    if n_branches < 2:
        return False, "Insufficient branching for a diterpenoid structure"

    # Check for a terpenoid backbone (e.g., isoprene units)
    isoprene_pattern = Chem.MolFromSmarts("CC(=C)C")
    if isoprene_pattern is not None and mol.HasSubstructMatch(isoprene_pattern):
        return True, "Contains a modified C20 skeleton with typical diterpenoid features (rings, functional groups, branching, isoprene units)"

    return True, "Contains a modified C20 skeleton with typical diterpenoid features (rings, functional groups, branching)"