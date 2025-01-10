"""
Classifies: CHEBI:36615 triterpenoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdchem import HybridizationType

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid based on its SMILES string.
    A triterpenoid generally contains around 30 carbon atoms and various modifications
    such as hydroxyl, ketone, and carboxylic acid groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Check if carbon count is close to 30 (allow some flexibility)
    if not (25 <= c_count <= 35):
        return False, f"Carbon count {c_count} not in triterpenoid range (25-35)"
    
    # Check for presence of oxygen and other functional groups
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    functional_groups_present = any(atom.GetHybridization() in {HybridizationType.SP2, HybridizationType.SP3} for atom in mol.GetAtoms())
    
    if o_count < 1:
        return False, "No oxygen atoms found; typical triterpenoids have oxygenated functional groups"
    
    if not functional_groups_present:
        return False, "No functional groups indicative of triterpenoid modifications (e.g., -OH, =O)"
    
    # Check for ring structures - triterpenoids usually have cyclic backbones
    ring_info = mol.GetRingInfo()
    if not ring_info.IsInitialized() or ring_info.NumRings() < 1:
        return False, "No ring structures found; triterpenoids usually have rings"
    
    # If all checks pass, classify as a triterpenoid
    return True, "Classification as triterpenoid based on carbon count and presence of functional groups and rings"