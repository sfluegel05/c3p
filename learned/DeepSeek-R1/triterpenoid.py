"""
Classifies: CHEBI:36615 triterpenoid
"""
"""
Classifies: CHEBI:XXXXX triterpenoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid based on its SMILES string.
    A triterpenoid is derived from a triterpene (C30 skeleton), possibly modified.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Count carbon atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 25:
        return False, f"Only {carbon_count} carbons, insufficient for triterpenoid"
    
    # Calculate molecular weight (approximate)
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecular weight ({mol_wt:.1f}) too low for triterpenoid"
    
    # Check for complex ring system (at least 3 rings)
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 3:
        return False, f"Only {num_rings} rings, insufficient for triterpenoid structure"
    
    # Optional: Check for isoprene-related patterns (simplified)
    # This is challenging; instead, rely on the above checks
    
    return True, "Meets triterpenoid criteria: sufficient carbons, molecular weight, and ring complexity"