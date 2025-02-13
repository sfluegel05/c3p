"""
Classifies: CHEBI:26935 tetraterpenoid
"""
"""
Classifies: CHEBI:25853 tetraterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    Tetraterpenoids are derived from tetraterpenes, which have a C40 skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetraterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Tetraterpenoids typically have around 40 carbon atoms
    if c_count < 35 or c_count > 45:
        return False, f"Carbon count ({c_count}) is not consistent with a tetraterpenoid (expected ~40)"

    # Check for long conjugated systems (characteristic of tetraterpenoids)
    conjugated_system_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    conjugated_matches = mol.GetSubstructMatches(conjugated_system_pattern)
    if len(conjugated_matches) < 8:
        return False, "Insufficient conjugated double bonds for a tetraterpenoid"

    # Check for oxygen-containing functional groups (common in tetraterpenoids)
    oxygen_pattern = Chem.MolFromSmarts("[OX2]")
    oxygen_matches = mol.GetSubstructMatches(oxygen_pattern)
    if len(oxygen_matches) < 1:
        return False, "No oxygen-containing functional groups found"

    # Check for cyclic structures (common in tetraterpenoids)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 1:
        return False, "No cyclic structures found, which are common in tetraterpenoids"

    # Check molecular weight (tetraterpenoids are typically large molecules)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for a tetraterpenoid"

    return True, "Molecule has a C40-like skeleton with conjugated double bonds and oxygen-containing functional groups, consistent with a tetraterpenoid"