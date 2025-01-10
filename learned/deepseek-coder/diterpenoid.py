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
    
    # Diterpenoids are derived from a C20 skeleton, but modifications may reduce the count
    if c_count < 15 or c_count > 25:
        return False, f"Carbon count ({c_count}) is outside the typical range for diterpenoids (15-25)"

    # Check for terpenoid-like structure (isoprene units or similar)
    # Terpenoids often have a high degree of branching and rings
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings < 1:
        return False, "No rings found, unlikely to be a diterpenoid"

    # Check for typical functional groups in diterpenoids (e.g., alcohols, ketones, esters)
    functional_groups = ["[OH]", "[C=O]", "[O]", "[C(=O)O]"]
    has_functional_group = any(mol.HasSubstructMatch(Chem.MolFromSmarts(group)) for group in functional_groups)
    if not has_functional_group:
        return False, "No typical diterpenoid functional groups found"

    # Check for isoprene-like patterns (C5H8 units)
    isoprene_pattern = Chem.MolFromSmarts("CC(=C)C")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 2:
        return False, "Fewer than 2 isoprene-like patterns found"

    # Check molecular weight - diterpenoids typically have MW between 250 and 400
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 400:
        return False, f"Molecular weight ({mol_wt}) is outside the typical range for diterpenoids (250-400)"

    return True, "Contains a C20 skeleton with typical diterpenoid features (rings, functional groups, isoprene-like patterns)"