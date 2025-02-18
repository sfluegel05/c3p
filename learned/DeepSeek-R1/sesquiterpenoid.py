"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_sesquiterpenoid(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Calculate number of carbons
    num_c = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (13 <= num_c <= 17):
        return False, f"Carbon count {num_c} not in 13-17"
    
    # Check for methyl groups or branching indicative of terpenoids
    methyl_pattern = Chem.MolFromSmarts("[CH3]")
    methyl_matches = len(mol.GetSubstructMatches(methyl_pattern))
    
    # Check for at least two methyl groups or branching (approximate)
    if methyl_matches < 2:
        return False, f"Only {methyl_matches} methyl groups"
    
    # Additional check for terpenoid-like structure (e.g., rings or double bonds)
    rings = Chem.GetSSSR(mol)
    double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    
    if rings < 1 and double_bonds < 2:
        return False, "Insufficient rings or double bonds for terpenoid"
    
    return True, f"{num_c} carbons, {methyl_matches} methyl groups, terpenoid features"