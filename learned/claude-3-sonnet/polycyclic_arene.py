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
    
    # Get all rings
    rings = ri.AtomRings()
    
    # Check each ring for aromaticity and carbon-only composition
    aromatic_carbocyclic_rings = []
    for ring in rings:
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        
        # Check if all atoms in ring are carbon
        if not all(atom.GetAtomicNum() == 6 for atom in ring_atoms):
            continue
            
        # Check if all atoms in ring are aromatic
        if not all(atom.GetIsAromatic() for atom in ring_atoms):
            continue
            
        # Check if all carbons are sp2 hybridized
        if not all(atom.GetHybridization() == Chem.HybridizationType.SP2 for atom in ring_atoms):
            continue
            
        aromatic_carbocyclic_rings.append(ring)
    
    if len(aromatic_carbocyclic_rings) < 2:
        return False, "Must contain at least 2 aromatic carbocyclic rings"
    
    # Create set of all aromatic ring atoms
    aromatic_ring_atoms = set()
    for ring in aromatic_carbocyclic_rings:
        aromatic_ring_atoms.update(ring)
    
    # Check if the aromatic system is the main part of the molecule
    total_heavy_atoms = mol.GetNumHeavyAtoms()
    if len(aromatic_ring_atoms) / total_heavy_atoms < 0.6:
        return False, "Aromatic system is not the main part of the molecule"
    
    # Check substituents - only allow H, OH, NO2, Cl
    allowed_substituents = {1, 6, 8, 7, 17}  # H, C, O, N, Cl
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_substituents:
            return False, f"Contains disallowed element: {atom.GetSymbol()}"
        
        # If it's not part of the aromatic system, it must be a substituent
        if atom.GetIdx() not in aromatic_ring_atoms:
            if atom.GetAtomicNum() not in {1, 8, 7, 17}:  # Only allow H, O, N, Cl as substituents
                return False, "Contains non-aromatic carbon framework"
    
    # Check ring connectivity using aromatic bonds
    ring_sets = [set(ring) for ring in aromatic_carbocyclic_rings]
    
    def are_rings_connected(ring1, ring2):
        """Check if two rings share atoms or are connected by aromatic bonds"""
        shared_atoms = ring1.intersection(ring2)
        if len(shared_atoms) >= 1:
            return True
        # Check for aromatic bond connections between rings
        for atom1 in ring1:
            for atom2 in ring2:
                bond = mol.GetBondBetweenAtoms(atom1, atom2)
                if bond is not None and bond.GetIsAromatic():
                    return True
        return False
    
    # Check if rings form a connected aromatic system
    visited_rings = set()
    def dfs_rings(ring_idx):
        visited_rings.add(ring_idx)
        for other_idx in range(len(ring_sets)):
            if other_idx not in visited_rings and are_rings_connected(ring_sets[ring_idx], ring_sets[other_idx]):
                dfs_rings(other_idx)
    
    dfs_rings(0)
    if len(visited_rings) != len(ring_sets):
        return False, "Aromatic rings must form a connected system"
    
    # Success - we have a polycyclic arene
    return True, "Contains multiple connected aromatic carbocyclic rings forming a polycyclic system"