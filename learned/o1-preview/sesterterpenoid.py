"""
Classifies: CHEBI:26660 sesterterpenoid
"""
"""
Classifies: CHEBI:26661 sesterterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    
    A sesterterpenoid is a terpenoid derived from a sesterterpene, which typically has 25 carbons
    built from five isoprene units (C5 units). The molecule may be rearranged or modified by 
    the removal of one or more skeletal atoms (generally methyl groups).
    
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
    
    # Count number of carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # For sesterterpenoids, expect around 25 carbons, but accept a wider range due to modifications
    if c_count < 20 or c_count > 55:
        return False, f"Carbon count is {c_count}, which is not typical for a sesterterpenoid (expected ~25 carbons)"
    
    # Check for presence of multiple double bonds (indicative of isoprene units)
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE)
    if double_bond_count < 2:
        return False, f"Found {double_bond_count} double bonds; sesterterpenoids typically have multiple double bonds due to isoprene units"
    
    # Check for presence of ring structures
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 1:
        return False, "No ring structures found; sesterterpenoids often contain rings"
    
    # Check for methyl groups attached to carbons (common in terpenoids)
    methyl_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetDegree() == 1:
            # Methyl group carbon
            neighbors = atom.GetNeighbors()
            if neighbors and neighbors[0].GetAtomicNum() == 6:
                methyl_count += 1
    if methyl_count < 1:
        return False, "No methyl groups found attached to carbons; sesterterpenoids often contain methyl groups"
    
    return True, "Molecule meets criteria for a sesterterpenoid"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:26661',
        'name': 'sesterterpenoid',
        'definition': 'Any terpenoid derived from a sesterterpene. The term includes compounds in which the C25 skeleton of the parent sesterterpene has been rearranged or modified by the removal of one or more skeletal atoms (generally methyl groups). Sometimes sesterterpenoids are erroneously referred to as sesterpenoids.',
        'parents': ['CHEBI:26873', 'CHEBI:26872']
    },
    'message': None,
    'success': True
}