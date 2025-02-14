"""
Classifies: CHEBI:36141 quinone
"""
from rdkit import Chem

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    A quinone is defined as compounds having a fully conjugated cyclic dione structure,
    derived from aromatic compounds by conversion of an even number of -CH= groups into
    -C(=O)- groups with any necessary rearrangement of double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    bond_rings = ring_info.BondRings()

    if not atom_rings:
        return False, "Molecule does not contain any rings"

    # Build ring systems (fused rings)
    # Each ring system is a set of ring indices that are connected via shared bonds
    def get_ring_systems(bond_rings):
        ring_count = len(bond_rings)
        adjacency = {i: set() for i in range(ring_count)}
        for i in range(ring_count):
            for j in range(i + 1, ring_count):
                # Check if rings i and j share any bonds
                if set(bond_rings[i]).intersection(bond_rings[j]):
                    adjacency[i].add(j)
                    adjacency[j].add(i)
        # Find connected components
        visited = set()
        ring_systems = []
        for i in range(ring_count):
            if i not in visited:
                stack = [i]
                component = set()
                while stack:
                    idx = stack.pop()
                    if idx not in visited:
                        visited.add(idx)
                        component.add(idx)
                        stack.extend(adjacency[idx] - visited)
                ring_systems.append(component)
        return ring_systems

    ring_systems_indices = get_ring_systems(bond_rings)

    # For each ring system, analyze if it satisfies quinone criteria
    for ring_sys_idxs in ring_systems_indices:
        # Collect atoms in the ring system
        ring_atoms = set()
        for ring_idx in ring_sys_idxs:
            ring_atoms.update(atom_rings[ring_idx])

        # Identify carbonyl carbons in the ring system
        carbonyl_carbons = []
        for idx in ring_atoms:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:
                # Check for C=O group
                is_carbonyl = False
                for nbr in atom.GetNeighbors():
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    if nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE:
                        is_carbonyl = True
                        break
                if is_carbonyl:
                    carbonyl_carbons.append(idx)
        if len(carbonyl_carbons) >= 2 and len(carbonyl_carbons) % 2 == 0:
            # Check if all atoms in the ring system are sp2 hybridized
            is_conjugated = True
            for idx in ring_atoms:
                atom = mol.GetAtomWithIdx(idx)
                # Consider only heavy atoms (exclude hydrogens)
                if atom.GetAtomicNum() > 1:
                    hyb = atom.GetHybridization()
                    if hyb not in (Chem.rdchem.HybridizationType.SP2, Chem.rdchem.HybridizationType.AROMATIC):
                        is_conjugated = False
                        break
            if not is_conjugated:
                continue  # Skip this ring system

            # Optionally, check if the ring system is planar (not enforced here)
            # Determine if the ring system was derived from an aromatic compound
            # For simplicity, we can accept that it is derived from an aromatic system

            return True, "Molecule contains a fully conjugated cyclic diketone ring system (quinone)"

    return False, "Molecule does not contain the characteristic quinone structure"