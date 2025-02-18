"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_sesquiterpenoid(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Calculate number of carbons (sesquiterpene base is C15, modifications may alter)
    num_c = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (12 <= num_c <= 20):  # Wider range to accommodate modifications
        return False, f"Carbon count {num_c} not in 12-20"
    
    # Check for methyl groups or branching (common in terpenoids)
    methyl_pattern = Chem.MolFromSmarts("[CH3]")
    methyl_matches = len(mol.GetSubstructMatches(methyl_pattern))
    
    # Require at least one methyl group (common in terpenoids)
    if methyl_matches < 1:
        return False, f"Only {methyl_matches} methyl groups"
    
    # Check for terpenoid features: rings or multiple double bonds
    ring_info = mol.GetRingInfo()
    rings = ring_info.NumRings()
    double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    
    # Allow for acyclic structures with multiple double bonds
    if rings == 0 and double_bonds < 2:
        return False, "Insufficient double bonds for acyclic sesquiterpenoid"
    
    # Check for isoprene-like patterns (approximate)
    isoprene_pattern = Chem.MolFromSmarts("CC(=C)C")  # Basic isoprene unit
    if not mol.HasSubstructMatch(isoprene_pattern):
        # May not work for rearranged structures, but a heuristic
        return False, "No isoprene-like subunits detected"
    
    return True, f"{num_c} carbons, {methyl_matches} methyl groups, terpenoid features"