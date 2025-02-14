"""
Classifies: CHEBI:50699 oligosaccharide
"""
from rdkit import Chem

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
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify monosaccharide units
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    monosaccharide_rings = []
    for ring in rings:
        ring_size = len(ring)
        if ring_size not in [5, 6]:
            continue  # Skip rings that are not 5 or 6-membered
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Count oxygen atoms in the ring
        oxygens_in_ring = [atom for atom in ring_atoms if atom.GetAtomicNum() == 8]
        if len(oxygens_in_ring) != 1:
            continue  # Monosaccharide rings have exactly one ring oxygen
        # Check for unsaturation in the ring (no double bonds)
        has_double_bond = False
        for bond in mol.GetBonds():
            if bond.IsInRing():
                if bond.GetBondTypeAsDouble() > 1.1:
                    if bond.GetBeginAtomIdx() in ring and bond.GetEndAtomIdx() in ring:
                        has_double_bond = True
                        break
        if has_double_bond:
            continue  # Skip rings with double bonds
        # Check if ring carbons have hydroxyl groups
        sufficient_hydroxyls = True
        for atom in ring_atoms:
            if atom.GetAtomicNum() == 6:  # Carbon atom
                hydroxyl_found = False
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:
                        hydroxyl_found = True  # Found hydroxyl group
                        break
                if not hydroxyl_found:
                    sufficient_hydroxyls = False
                    break
        if not sufficient_hydroxyls:
            continue  # Skip if not all carbons have hydroxyls
        # If all checks passed, it's a monosaccharide ring
        monosaccharide_rings.append(set(ring))
    
    num_monosaccharides = len(monosaccharide_rings)
    if num_monosaccharides == 0:
        return False, "No monosaccharide units found"
    
    # Check for glycosidic linkages between monosaccharide units
    glycosidic_bonds = []
    for bond in mol.GetBonds():
        if bond.IsInRing():
            continue  # Skip ring bonds
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        # Look for ether linkages (C-O-C)
        if atom1.GetAtomicNum() == 8 and atom2.GetAtomicNum() == 6:
            oxygen_atom = atom1
            carbon_atom = atom2
        elif atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 8:
            carbon_atom = atom1
            oxygen_atom = atom2
        else:
            continue  # Not an ether linkage
        # Check if oxygen is connected to two carbons
        if oxygen_atom.GetDegree() != 2:
            continue
        neighbors = oxygen_atom.GetNeighbors()
        if len(neighbors) != 2:
            continue
        carbon_neighbors = [n for n in neighbors if n.GetAtomicNum() == 6]
        if len(carbon_neighbors) != 2:
            continue
        # Check if carbons are part of different monosaccharide rings
        rings_with_carbons = []
        for ring_set in monosaccharide_rings:
            if carbon_neighbors[0].GetIdx() in ring_set or carbon_neighbors[1].GetIdx() in ring_set:
                rings_with_carbons.append(ring_set)
        if len(set(map(tuple, rings_with_carbons))) >= 2:
            glycosidic_bonds.append(bond)
    
    num_glycosidic_bonds = len(glycosidic_bonds)
    if num_glycosidic_bonds < num_monosaccharides - 1:
        return False, "Not all monosaccharide units are connected via glycosidic linkages"
    
    # Ensure the number of monosaccharide units is within the typical range for oligosaccharides
    if num_monosaccharides < 2:
        return False, "Less than 2 monosaccharide units found"
    elif num_monosaccharides > 10:
        return False, "More than 10 monosaccharide units found; may be a polysaccharide"
    
    return True, f"Contains {num_monosaccharides} monosaccharide units connected via glycosidic linkages"