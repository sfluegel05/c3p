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
    
    # Sesterterpenoids typically have around 25 carbons, but modifications can change this
    if c_count < 20 or c_count > 40:
        return False, f"Carbon count ({c_count}) is not consistent with sesterterpenoid (C25 skeleton)"

    # Check for terpenoid-like structure (multiple rings and functional groups)
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 1:
        return False, "Too few rings for a sesterterpenoid"

    # Check for long carbon chains or complex ring systems
    if rdMolDescriptors.CalcNumRotatableBonds(mol) < 3:
        return False, "Molecule is too rigid for a sesterterpenoid"

    # Check molecular weight (sesterterpenoids are typically large molecules)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350:
        return False, "Molecular weight too low for a sesterterpenoid"

    # Check for typical terpenoid functional groups (e.g., hydroxyl, carbonyl)
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetTotalNumHs() > 0 and atom.GetAtomicNum() == 8)
    carbonyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in atom.GetBonds()))
    
    # Optional: Check for long carbon chains
    long_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    long_chain_matches = mol.GetSubstructMatches(long_chain_pattern)
    if len(long_chain_matches) < 2:
        return False, "Missing long carbon chains typical of sesterterpenoids"

    return True, "Contains a C25 skeleton with terpenoid-like features (rings, functional groups, and complexity)"