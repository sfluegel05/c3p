"""
Classifies: CHEBI:35915 sterol ester
"""
"""
Classifies: sterol ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    A sterol ester is a steroid ester obtained by formal condensation of the carboxy group of any carboxylic acid with the 3-hydroxy group of a sterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Find rings of size 5 and 6
    ring_sizes = [len(r) for r in atom_rings]
    six_membered_rings = [r for r in atom_rings if len(r) == 6]
    five_membered_rings = [r for r in atom_rings if len(r) == 5]

    if len(six_membered_rings) < 3 or len(five_membered_rings) < 1:
        return False, "Does not contain required number of rings for sterol nucleus"

    # Build ring adjacency map
    # The ring adjacency can be determined by checking the bonds shared between rings
    ring_bond_map = {}
    for idx1, ring1 in enumerate(atom_rings):
        bonds1 = set()
        for i in range(len(ring1)):
            a1 = ring1[i]
            a2 = ring1[(i+1)%len(ring1)]
            bond = mol.GetBondBetweenAtoms(a1, a2).GetIdx()
            bonds1.add(bond)
        for idx2, ring2 in enumerate(atom_rings):
            if idx1 >= idx2:
                continue
            bonds2 = set()
            for i in range(len(ring2)):
                a1 = ring2[i]
                a2 = ring2[(i+1)%len(ring2)]
                bond = mol.GetBondBetweenAtoms(a1, a2).GetIdx()
                bonds2.add(bond)
            shared_bonds = bonds1.intersection(bonds2)
            if shared_bonds:
                ring_bond_map.setdefault(idx1, set()).add(idx2)
                ring_bond_map.setdefault(idx2, set()).add(idx1)

    # Identify fused ring system of size 6-6-6-5
    # Try to find a set of four rings that are connected and have sizes 6-6-6-5
    found_steroid_nucleus = False
    for idx1 in range(len(atom_rings)):
        for idx2 in ring_bond_map.get(idx1, []):
            for idx3 in ring_bond_map.get(idx2, []):
                for idx4 in ring_bond_map.get(idx3, []):
                    ring_sizes_sequence = [len(atom_rings[idx]) for idx in [idx1, idx2, idx3, idx4]]
                    if sorted(ring_sizes_sequence) == [5,6,6,6]:
                        unique_indices = set([idx1, idx2, idx3, idx4])
                        if len(unique_indices) == 4:
                            # Check if they are all connected
                            connections = [
                                idx2 in ring_bond_map[idx1],
                                idx3 in ring_bond_map[idx2],
                                idx4 in ring_bond_map[idx3],
                            ]
                            if all(connections):
                                found_steroid_nucleus = True
                                steroid_ring_indices = [idx1, idx2, idx3, idx4]
                                break
                if found_steroid_nucleus:
                    break
            if found_steroid_nucleus:
                break
        if found_steroid_nucleus:
            break

    if not found_steroid_nucleus:
        return False, "No steroid nucleus found"

    # Get atoms in steroid nucleus
    steroid_atom_indices = set()
    for idx in steroid_ring_indices:
        steroid_atom_indices.update(atom_rings[idx])

    # Find oxygen atoms attached to steroid nucleus
    steroid_atoms = [mol.GetAtomWithIdx(idx) for idx in steroid_atom_indices]
    steroid_atom_indices_set = set(steroid_atom_indices)
    oxygen_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in steroid_atom_indices_set:
                    oxygen_atoms.append(atom)
                    break

    if not oxygen_atoms:
        return False, "No oxygen atoms attached to steroid nucleus"

    # Check if any of these oxygen atoms are part of an ester group
    ester_pattern = Chem.MolFromSmarts('[CX3](=O)[OX2H0]')  # Ester group
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    ester_oxygen_indices = [match[2] for match in ester_matches]  # index of oxygen atom in ester

    for o_atom in oxygen_atoms:
        if o_atom.GetIdx() in ester_oxygen_indices:
            return True, "Contains steroid nucleus with esterified oxygen attached"

    return False, "No esterified oxygen attached to steroid nucleus found"