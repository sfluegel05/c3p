"""
Classifies: CHEBI:15889 sterol
"""
"""
Classifies: sterol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.
    A sterol is any 3-hydroxy steroid whose skeleton is closely related to cholestan-3-ol
    (additional carbon atoms may be present in the side chain).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Extract the Murcko scaffold (core steroid backbone)
    core = MurckoScaffold.GetScaffoldForMol(mol)
    if core is None:
        return False, "Cannot extract molecular scaffold"
    
    # Get ring information
    ring_info = core.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings != 4:
        return False, f"Scaffold does not have 4 rings, found {num_rings}"
    
    # Get sizes of the rings
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    if sorted(ring_sizes) != [5, 6, 6, 6]:
        return False, f"Ring sizes are not [5, 6, 6, 6], found {sorted(ring_sizes)}"
    
    # Check if rings are fused correctly (i.e., form the steroid nucleus)
    if not Chem.FusedRingUtils.IsFused(core):
        return False, "Rings are not fused correctly to form steroid nucleus"
    
    # Define 3-hydroxy group attached to ring A
    hydroxy_pattern = Chem.MolFromSmarts('[#6]-1(-[#8H])-[#6]=[#6]-[#6]-[#6]-1')  # 3-hydroxy group on ring A
    if hydroxy_pattern is None:
        return False, "Invalid hydroxy group SMARTS pattern"

    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "Does not have 3-hydroxy group at position 3"
    
    return True, "Contains steroid backbone with 3-hydroxy group characteristic of sterols"

__metadata__ = {
    'chemical_class': {
        'name': 'sterol',
        'definition': 'Any 3-hydroxy steroid whose skeleton is closely related to cholestan-3-ol (additional carbon atoms may be present in the side chain).'
    }
}