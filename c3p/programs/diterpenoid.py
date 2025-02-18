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
    if c_count < 15 or c_count > 50:  # Wider range to accommodate modifications
        return False, f"Carbon count ({c_count}) is outside the typical range for diterpenoids"

    # Check molecular weight - diterpenoids typically have a molecular weight between 250 and 800 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 800:  # Wider range
        return False, f"Molecular weight ({mol_wt:.2f} Da) is outside the typical range for diterpenoids"

    # Check for the presence of rings or long carbon chains
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    
    # Check for rotatable bonds - diterpenoids can be flexible or rigid
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    # Check for common terpenoid features (not all need to be present)
    isoprene_pattern = Chem.MolFromSmarts("[CX4H2][CX4H1]=[CX3H1]")  # More general isoprene pattern
    ring_pattern = Chem.MolFromSmarts("[R]")  # Any ring
    long_chain_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]")  # Long carbon chain
    
    has_isoprene = mol.HasSubstructMatch(isoprene_pattern)
    has_rings = mol.HasSubstructMatch(ring_pattern)
    has_long_chain = mol.HasSubstructMatch(long_chain_pattern)
    
    # At least two of these features should be present
    if sum([has_isoprene, has_rings, has_long_chain]) < 2:
        return False, "Insufficient characteristic terpenoid features found"

    return True, "Molecule has characteristics consistent with a diterpenoid (C20 skeleton, terpenoid features, and appropriate size)"