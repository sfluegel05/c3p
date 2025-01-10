"""
Classifies: CHEBI:50126 tetrasaccharide
"""
"""
Classifies: tetrasaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetrasaccharide(smiles: str):
    """
    Determines if a molecule is a tetrasaccharide based on its SMILES string.
    A tetrasaccharide is an oligosaccharide comprising four monomeric monosaccharide units.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a tetrasaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information
    ring_info = mol.GetRingInfo()
    
    # List to store monosaccharide rings
    mono_rings = []
    
    # Loop over rings
    for ring_atoms in ring_info.AtomRings():
        ring = [mol.GetAtomWithIdx(idx) for idx in ring_atoms]
        ring_size = len(ring)
        
        # Check if ring is 5 or 6 membered
        if ring_size not in [5,6]:
            continue
        
        # Count number of oxygens in ring
        num_oxygens = sum(1 for atom in ring if atom.GetAtomicNum() == 8)
        # Count number of carbons in ring
        num_carbons = sum(1 for atom in ring if atom.GetAtomicNum() == 6)
        
        # Check for one oxygen and rest carbons
        if num_oxygens == 1 and num_carbons == ring_size - 1:
            mono_rings.append(set(ring_atoms))
    
    num_mono_units = len(mono_rings)
    
    if num_mono_units != 4:
        return False, f"Found {num_mono_units} monosaccharide units, expected 4"
    
    # Check connectivity between monosaccharide units (glycosidic linkages)
    # Build a graph of rings connected via oxygen atoms (glycosidic oxygen bridges)
    edges = []
    for i, ring1 in enumerate(mono_rings):
        for j, ring2 in enumerate(mono_rings):
            if i >= j:
                continue
            # Check for bonds between ring1 and ring2 via oxygen atoms
            for atom_idx1 in ring1:
                atom1 = mol.GetAtomWithIdx(atom_idx1)
                for bond in atom1.GetBonds():
                    neighbor = bond.GetOtherAtom(atom1)
                    neighbor_idx = neighbor.GetIdx()
                    # If neighbor is in ring2 and atom is oxygen
                    if neighbor_idx in ring2 and atom1.GetAtomicNum() == 8:
                        edges.append((i,j))
    # Build connectivity graph
    from collections import defaultdict, deque
    
    graph = defaultdict(list)
    for i,j in edges:
        graph[i].append(j)
        graph[j].append(i)
    
    # Check if the monosaccharide units are connected
    # For tetrasaccharide, we expect a connected graph of 4 nodes
    visited = [False]*num_mono_units
    queue = deque()
    queue.append(0)
    visited[0] = True
    while queue:
        current = queue.popleft()
        for neighbor in graph[current]:
            if not visited[neighbor]:
                visited[neighbor] = True
                queue.append(neighbor)
    if not all(visited):
        return False, "Monosaccharide units are not properly connected via glycosidic linkages"
    
    # Additional checks (optional)
    # Check molecular weight range for typical tetrasaccharides
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500 or mol_wt > 1000:
        return False, f"Molecular weight {mol_wt:.2f} not typical for a tetrasaccharide"
    
    return True, "Molecule contains 4 monosaccharide units connected via glycosidic linkages"