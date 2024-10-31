from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from collections import defaultdict

def is_oxa_steroid(smiles: str):
    """
    Determines if a molecule is an oxa-steroid (steroid where a carbon atom is replaced by oxygen)
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an oxa-steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
        
    # Check for presence of oxygen
    if not any(atom.GetAtomicNum() == 8 for atom in mol.GetAtoms()):
        return False, "No oxygen atoms present"
    
    # Get ring information
    ring_info = mol.GetRingInfo()
    if not ring_info.NumRings() >= 3:
        return False, "Less than 3 rings present"
        
    # Get all rings and their atoms
    rings = ring_info.AtomRings()
    
    # Create graph of ring connectivity
    ring_graph = defaultdict(set)
    for i, ring1 in enumerate(rings):
        ring1_atoms = set(ring1)
        for j, ring2 in enumerate(rings):
            if i != j:
                ring2_atoms = set(ring2)
                if ring1_atoms & ring2_atoms:  # If rings share atoms
                    ring_graph[i].add(j)
                    ring_graph[j].add(i)
    
    # Find largest connected ring system
    def get_connected_rings(start, visited=None):
        if visited is None:
            visited = set()
        visited.add(start)
        for neighbor in ring_graph[start]:
            if neighbor not in visited:
                get_connected_rings(neighbor, visited)
        return visited
        
    largest_system = max(len(get_connected_rings(i)) for i in range(len(rings)))
    if largest_system < 3:
        return False, "No steroid-like ring system found"
        
    # Check for oxygen atoms in rings
    ring_atoms = set()
    for ring in rings:
        ring_atoms.update(ring)
        
    ring_oxygens = sum(1 for atom_idx in ring_atoms 
                      if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 8)
    
    if ring_oxygens == 0:
        return False, "No oxygen atoms in ring system"
        
    # Check for tetracyclic or tricyclic system with appropriate ring sizes
    ring_sizes = sorted([len(ring) for ring in rings])
    valid_sizes = [5,6]
    if not any(size in valid_sizes for size in ring_sizes):
        return False, "Does not contain appropriate ring sizes"

    # Additional check for steroid-like carbon skeleton
    carbon_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
    if carbon_count < 15:  # Steroids typically have 17+ carbons
        return False, "Too few carbons for steroid skeleton"

    return True, f"Contains steroid-like ring system with {ring_oxygens} oxygen atoms in rings"
# Pr=0.782608695652174
# Recall=1.0