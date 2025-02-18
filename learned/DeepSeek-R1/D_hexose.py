"""
Classifies: CHEBI:4194 D-hexose
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    A D-hexose has six carbons and D-configuration at position 5.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for exactly 6 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 6:
        return False, f"Expected 6 carbons, found {c_count}"
    
    # Check for multiple hydroxyl groups (at least 4)
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() >= 1)
    if hydroxyl_count < 4:
        return False, f"Expected at least 4 hydroxyls, found {hydroxyl_count}"
    
    # Find chiral centers
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if not chiral_centers:
        return False, "No chiral centers found"
    
    # Assuming the last chiral center is position 5 (heuristic)
    # For D-hexose, the configuration at position 5 should be R (example-based)
    # This is a simplification and may not be accurate for all cases
    last_center = chiral_centers[-1]
    atom_idx, cip = last_center
    if cip != 'R':
        return False, f"Chiral center at atom {atom_idx} is {cip}, expected R"
    
    return True, "D-configuration at position 5 (heuristic check)"