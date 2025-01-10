"""
Classifies: CHEBI:36615 triterpenoid
"""
"""
Classifies: triterpenoids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a triterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Triterpenoids typically have around 30 carbons (allowing some variation)
    if c_count < 25 or c_count > 35:
        return False, f"Carbon count {c_count} outside typical range for triterpenoids (25-35)"
    
    # Count rings
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    
    if ring_count < 4:
        return False, f"Too few rings ({ring_count}) for a triterpenoid"
        
    # Look for fused ring systems
    sssr = Chem.GetSymmSSSR(mol)
    if len(sssr) < 4:
        return False, "Insufficient ring systems for a triterpenoid"
        
    # Count methyl groups (-CH3)
    methyl_pattern = Chem.MolFromSmarts('[CH3]')
    methyl_count = len(mol.GetSubstructMatches(methyl_pattern))
    
    if methyl_count < 4:
        return False, f"Too few methyl groups ({methyl_count}) for a triterpenoid"
        
    # Check molecular weight - should be substantial
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for triterpenoid"
        
    # Look for characteristic carbon framework
    # Pattern for 6-membered ring
    six_memb_ring = Chem.MolFromSmarts('C1CCCCC1')
    if not mol.HasSubstructMatch(six_memb_ring):
        return False, "Missing characteristic ring structure"
        
    # Count sp3 carbons (should be many in a triterpenoid)
    sp3_c_pattern = Chem.MolFromSmarts('[CX4]')
    sp3_c_count = len(mol.GetSubstructMatches(sp3_c_pattern))
    if sp3_c_count < 15:
        return False, f"Too few sp3 carbons ({sp3_c_count}) for triterpenoid skeleton"
        
    # Calculate number of rings and rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    # Typical triterpenoid characteristics
    if ring_count >= 4 and methyl_count >= 4 and sp3_c_count >= 15:
        reason = (f"Has characteristic triterpenoid features: {ring_count} rings, "
                 f"{methyl_count} methyl groups, {c_count} carbons, "
                 f"{sp3_c_count} sp3 carbons")
        return True, reason
        
    return False, "Does not match typical triterpenoid characteristics"