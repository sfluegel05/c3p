"""
Classifies: CHEBI:26660 sesterterpenoid
"""
"""
Classifies: CHEBI:36541 sesterterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    A sesterterpenoid is derived from a sesterterpene (C25 skeleton), which may be rearranged or modified.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesterterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Sesterterpenoids typically have around 25 carbons, allow some modifications
    if not (20 <= c_count <= 40):
        return False, f"Carbon count ({c_count}) is not consistent with sesterterpenoid (C25 skeleton)"

    # Check molecular weight (sesterterpenoids are typically large molecules)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for a sesterterpenoid"

    # Check for terpenoid-like structure (multiple rings and functional groups)
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 1:
        return False, "Too few rings for a sesterterpenoid"

    # Check for complex ring systems
    if max(len(r) for r in ring_info.AtomRings()) < 4:
        return False, "Ring systems too simple for a sesterterpenoid"

    # Check for typical terpenoid functional groups
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 1:
        return False, "Too few oxygen atoms for a sesterterpenoid"

    # Check for terpenoid-like branching pattern (multiple methyl groups)
    methyl_pattern = Chem.MolFromSmarts("[CH3]")
    methyl_matches = mol.GetSubstructMatches(methyl_pattern)
    if len(methyl_matches) < 2:
        return False, "Missing terpenoid-like branching pattern (multiple methyl groups)"

    # Check for long carbon chains or complex structure
    if rdMolDescriptors.CalcNumRotatableBonds(mol) < 2:
        return False, "Molecule is too rigid for a sesterterpenoid"

    return True, "Contains a C25 skeleton with terpenoid-like features (rings, functional groups, and complexity)"