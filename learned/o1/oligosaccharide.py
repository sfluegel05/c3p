"""
Classifies: CHEBI:50699 oligosaccharide
"""
from rdkit import Chem

def is_sugar_ring(mol, ring):
    """
    Check if a ring is a sugar ring.
    A sugar ring is defined as a 5 or 6 membered ring containing exactly one ring oxygen atom,
    and all ring carbons have at least one exocyclic oxygen atom (e.g., hydroxyl or other oxygen-containing groups).
    """
    ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
    ring_size = len(ring_atoms)
    if ring_size not in [5, 6]:
        return False
    num_ring_oxygens = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8)
    if num_ring_oxygens != 1:
        return False
    carbons_in_ring = [atom for atom in ring_atoms if atom.GetAtomicNum() == 6]
    if len(carbons_in_ring) != ring_size - 1:
        return False
    # Check that each ring carbon is connected to at least one oxygen not in the ring (exocyclic oxygen)
    for carbon in carbons_in_ring:
        exocyclic_oxygens = [nbr for nbr in carbon.GetNeighbors()
                             if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in ring]
        if len(exocyclic_oxygens) == 0:
            return False
    return True

def is_oligosaccharide(smiles: str):
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string.
    An oligosaccharide is a compound in which monosaccharide units are joined by glycosidic linkages.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an oligosaccharide, False otherwise
        str: Reason for classification
    """
    from rdkit.Chem import rdMolDescriptors
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Kekulize the molecule to ensure all valencies are handled correctly
    try:
        Chem.Kekulize(mol)
    except:
        pass  # Some molecules cannot be kekulized; ignore for this context
    
    # Identify monosaccharide units
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    monosaccharide_rings = []
    ring_atom_sets = []  # Keep track of ring atom sets
    for ring in rings:
        if is_sugar_ring(mol, ring):
            monosaccharide_rings.append(set(ring))
            ring_atom_sets.append(set(ring))
    
    num_monosaccharides = len(monosaccharide_rings)
    if num_monosaccharides == 0:
        return False, "No monosaccharide units found"
    
    # Build atom to ring index mapping
    atom_ring_indices = {}
    for idx, ring_set in enumerate(ring_atom_sets):
        for atom_idx in ring_set:
            atom_ring_indices.setdefault(atom_idx, []).append(idx)
    
    # Check for glycosidic linkages between monosaccharide units
    glycosidic_bonds = []
    for bond in mol.GetBonds():
        if bond.IsInRing():
            continue  # Skip ring bonds
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        # Look for ether linkages (C-O-C)
        if set([atom1.GetAtomicNum(), atom2.GetAtomicNum()]) != {6, 8}:
            continue  # Not a C-O bond
        oxygen_atom = atom1 if atom1.GetAtomicNum() == 8 else atom2
        carbon_atom = atom1 if atom1.GetAtomicNum() == 6 else atom2
        # Check if oxygen is connected to two carbons
        if oxygen_atom.GetDegree() != 2:
            continue
        neighbors = oxygen_atom.GetNeighbors()
        carbon_neighbors = [n for n in neighbors if n.GetAtomicNum() == 6]
        if len(carbon_neighbors) != 2:
            continue
        # Check if carbons are part of different monosaccharide rings
        ring_indices = []
        for carbon_neighbor in carbon_neighbors:
            atom_idx = carbon_neighbor.GetIdx()
            rings_for_atom = atom_ring_indices.get(atom_idx, [])
            ring_indices.extend(rings_for_atom)
        ring_indices = set(ring_indices)
        if len(ring_indices) >= 2:
            glycosidic_bonds.append(bond)
    
    num_glycosidic_bonds = len(glycosidic_bonds)
    if num_glycosidic_bonds < num_monosaccharides - 1:
        return False, "Not all monosaccharide units are connected via glycosidic linkages"
    
    # Ensure the number of monosaccharide units is within the typical range for oligosaccharides
    if num_monosaccharides < 2:
        return False, "Less than 2 monosaccharide units found"
    elif num_monosaccharides > 10:
        return False, "More than 10 monosaccharide units found; may be a polysaccharide"
    
    # Optional: Check if all substituents are acceptable (e.g., acetamido groups)
    # This can be enhanced based on the specific definition and exclusions
    
    return True, f"Contains {num_monosaccharides} monosaccharide units connected via glycosidic linkages"