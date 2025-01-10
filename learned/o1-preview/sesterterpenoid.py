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
    if c_count < 22 or c_count > 28:
        return False, f"Carbon count is {c_count}, which is not typical for a sesterterpenoid (expected ~25 carbons)"
    
    # Check for the presence of rings (common in terpenoids)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 1:
        return False, "No ring structures found, unlikely to be a terpenoid"
    
    # Check for oxygen atoms (common in terpenoids)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1:
        return False, "No oxygen atoms found, terpenoids often contain oxygen"
    
    # Check for common functional groups in terpenoids
    # Hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl groups found, which are common in terpenoids"
    
    # Carbonyl group
    carbonyl_pattern = Chem.MolFromSmarts("C=O")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "No carbonyl groups found, which are common in terpenoids"
    
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 600:
        return False, f"Molecular weight is {mol_wt:.2f}, which is not typical for a sesterterpenoid"
    
    # Check for terpenoid-like skeleton using fingerprint similarity (optional)
    # This can be implemented if more advanced methods are allowed
    
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