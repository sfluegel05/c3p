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
    if c_count < 20 or c_count > 30:
        return False, f"Carbon count is {c_count}, which is not typical for a sesterterpenoid (expected ~25 carbons)"

    # Look for isoprene units (C=C-C-C=C)
    isoprene_pattern = Chem.MolFromSmarts("C=C-C-C=C")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 2:
        return False, f"Found {len(isoprene_matches)} isoprene units, less than expected for sesterterpenoid (at least 2)"

    # Check for oxygen atoms (common in terpenoids)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count == 0:
        return False, "No oxygen atoms found, not a terpenoid"
    
    # Check for ring structures
    ring_info = mol.GetRingInfo()
    if not ring_info.IsInitialized() or ring_info.NumRings() == 0:
        return False, "No ring structures found, unlikely to be a terpenoid"
    
    # Additional terpenoid features (optional)
    # Check for hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) == 0:
        return False, "No hydroxyl groups found, which are common in terpenoids"
    
    # Passed all checks
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