"""
Classifies: CHEBI:16219 cucurbitacin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin based on its SMILES string.
    Cucurbitacins are tetracyclic triterpenoids derived from cucurbitane.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a cucurbitacin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the tetracyclic triterpenoid structure, generically as 4 rings
    ring_info = mol.GetRingInfo()
    if not ring_info.IsInitialized():
        return False, "Ring information could not be determined"
    
    num_rings = len(ring_info.AtomRings())
    if num_rings < 4:
        return False, "Less than 4 rings, does not meet tetracyclic criteria"
    
    # Simplistic checks for key functional groups
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetAtomMapNum() is None)
    if hydroxyl_count < 1:
        return False, "Does not have expected hydroxyl groups"

    carbonyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetAtomMapNum() == 1)
    if carbonyl_count < 1:
        return False, "Does not have expected carbonyl groups"
    
    # Count total carbon atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 20:
        return False, "Too few carbons to represent common cucurbitacins"
    
    # Optionally, match specific substructure patterns if known
    
    return True, "Contains structural features typical of cucurbitacins"