"""
Classifies: CHEBI:33848 polycyclic arene
"""
"""
Classifies: polycyclic arene (polycyclic aromatic hydrocarbon)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polycyclic_arene(smiles: str):
    """
    Determines if a molecule is a polycyclic arene based on its SMILES string.
    A polycyclic arene is a polycyclic aromatic hydrocarbon (PAH).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polycyclic arene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Get ring information
    ri = mol.GetRingInfo()
    if ri.NumRings() < 2:
        return False, "Must contain at least 2 rings"
    
    # Check if molecule has aromatic atoms
    if not any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "No aromatic rings found"
    
    # Get all aromatic rings (both 5 and 6 membered)
    aromatic_rings = []
    for ring in ri.AtomRings():
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_rings.append(ring)
    
    if len(aromatic_rings) < 2:
        return False, "Must contain at least 2 aromatic rings"
    
    # Count heteroatoms and check if they're acceptable for a PAH
    heteroatom_count = 0
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in {1, 6}:  # Only H and C allowed in core structure
            # Allow some substituents (OH, NO2, etc) but count them
            if atomic_num not in {7, 8, 9, 17}:  # N, O, F, Cl allowed as substituents
                return False, "Contains disallowed elements"
            heteroatom_count += 1
    
    # Check carbon proportion - should be primarily hydrocarbon
    total_heavy_atoms = mol.GetNumHeavyAtoms()
    carbon_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    if carbon_atoms / total_heavy_atoms < 0.75:
        return False, "Too few carbon atoms for a polycyclic arene"
    
    if heteroatom_count / total_heavy_atoms > 0.2:
        return False, "Too many heteroatoms for a polycyclic arene"
    
    # Check ring connectivity
    ring_sets = [set(ring) for ring in aromatic_rings]
    
    def are_rings_connected(ring1, ring2):
        """Check if two rings share atoms or are connected by a bond"""
        shared_atoms = ring1.intersection(ring2)
        if len(shared_atoms) >= 1:
            return True
        # Check for single bond connections between rings
        for atom1 in ring1:
            for atom2 in ring2:
                if mol.GetBondBetweenAtoms(atom1, atom2) is not None:
                    return True
        return False
    
    # Build connectivity graph
    n_rings = len(ring_sets)
    connected = [[False] * n_rings for _ in range(n_rings)]
    for i in range(n_rings):
        for j in range(i+1, n_rings):
            if are_rings_connected(ring_sets[i], ring_sets[j]):
                connected[i][j] = connected[j][i] = True
    
    # Check if rings form a connected system using DFS
    visited = [False] * n_rings
    def dfs(v):
        visited[v] = True
        for u in range(n_rings):
            if connected[v][u] and not visited[u]:
                dfs(u)
    
    dfs(0)
    if not all(visited):
        return False, "Aromatic rings must form a connected system"
    
    # Check if at least one pair of rings is fused
    has_fused_rings = False
    for i in range(len(ring_sets)):
        for j in range(i+1, len(ring_sets)):
            shared_atoms = ring_sets[i].intersection(ring_sets[j])
            if len(shared_atoms) >= 2:
                # Verify shared atoms are adjacent
                shared_list = list(shared_atoms)
                if mol.GetBondBetweenAtoms(shared_list[0], shared_list[1]) is not None:
                    has_fused_rings = True
                    break
        if has_fused_rings:
            break
    
    if not has_fused_rings:
        # Check if it's a multi-ring system connected by single bonds
        aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
        if aromatic_atoms / total_heavy_atoms < 0.8:
            return False, "Non-fused aromatic system must be primarily composed of aromatic rings"
    
    return True, "Contains multiple aromatic rings forming a polycyclic system"