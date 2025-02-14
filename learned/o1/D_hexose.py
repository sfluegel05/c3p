"""
Classifies: CHEBI:4194 D-hexose
"""
"""
Classifies: CHEBI:16551 D-hexose
"""
from rdkit import Chem

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    A D-hexose is a monosaccharide with six carbons and D-configuration at position 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a D-hexose, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add hydrogens
    mol = Chem.AddHs(mol)

    # Assign stereochemistry
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)

    # Find all carbon atoms
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbons) < 6:
        return False, f"Number of carbon atoms is {len(carbons)}, expected at least 6 for hexose"

    # Find aldehyde group (C=O) to identify open-chain form
    aldehyde_smarts = '[CX3H1](=O)'
    aldehyde = Chem.MolFromSmarts(aldehyde_smarts)
    aldehyde_matches = mol.GetSubstructMatches(aldehyde)
    if aldehyde_matches:
        # Open-chain form detected
        c1_idx = aldehyde_matches[0][0]
        # Traverse the chain to number the carbons from C1 to C6
        carbon_chain = [c1_idx]
        visited = set(carbon_chain)
        current_idx = c1_idx

        while len(carbon_chain) < 6:
            atom = mol.GetAtomWithIdx(current_idx)
            neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors()
                         if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited]
            if neighbors:
                current_idx = neighbors[0]
                carbon_chain.append(current_idx)
                visited.add(current_idx)
            else:
                break  # No further carbons found

        if len(carbon_chain) != 6:
            return False, "Could not find a chain of 6 connected carbons starting from the aldehyde carbon"

        # Get the C5 atom (5th carbon in the chain)
        c5_idx = carbon_chain[4]
        c5_atom = mol.GetAtomWithIdx(c5_idx)

        # Check if C5 is chiral and get its configuration
        if not c5_atom.HasProp('_CIPCode'):
            return False, "C5 is not a chiral center"
        c5_chirality = c5_atom.GetProp('_CIPCode')
        if c5_chirality == 'R':
            return True, "Molecule is a D-hexose; C5 has R configuration"
        elif c5_chirality == 'S':
            return False, "Molecule is an L-hexose; C5 has S configuration"
        else:
            return False, "Chirality at C5 is undefined"
    else:
        # Cyclic form detected

        # Find rings containing oxygen
        rings = mol.GetRingInfo().AtomRings()
        o_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8]
        potential_rings = []
        for ring in rings:
            if any(atom_idx in o_atoms for atom_idx in ring):
                # Ring contains oxygen
                ring_size = len(ring)
                if ring_size == 5 or ring_size == 6:
                    # Potential furanose or pyranose
                    potential_rings.append(ring)

        # For each potential ring, attempt to find the anomeric carbon
        for ring in potential_rings:
            # Find oxygen atom in the ring
            o_in_ring = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
            if len(o_in_ring) !=1:
                continue  # Skip if multiple oxygens in ring
            o_idx = o_in_ring[0]
            # Find the anomeric carbon (connected to ring oxygen and exocyclic oxygen)
            for atom_idx in ring:
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetAtomicNum() == 6:
                    neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors()]
                    if o_idx in neighbors:
                        # Atom is connected to the ring oxygen
                        # Check if it's also connected to an exocyclic oxygen (OH group)
                        exocyclic_oxygens = [nbr.GetIdx() for nbr in atom.GetNeighbors()
                                             if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in ring]
                        if exocyclic_oxygens:
                            # Found anomeric carbon (C1)
                            c1_idx = atom_idx
                            # Now assign numbering around the ring
                            # Number carbons from C1 to C5
                            ring_atoms = ring.copy()
                            ring_atoms.remove(o_idx)
                            # Create an ordered list starting from c1_idx
                            ring_order = []
                            current_idx = c1_idx
                            previous_idx = o_idx
                            for _ in range(len(ring_atoms)):
                                ring_order.append(current_idx)
                                current_atom = mol.GetAtomWithIdx(current_idx)
                                # Find next carbon in ring
                                neighbors = [nbr.GetIdx() for nbr in current_atom.GetNeighbors()
                                             if nbr.GetIdx() in ring_atoms and nbr.GetIdx() != previous_idx]
                                if neighbors:
                                    previous_idx = current_idx
                                    current_idx = neighbors[0]
                                else:
                                    break
                            # Now ring_order should have the carbons in order from C1
                            if len(ring_order) >= 5:
                                c5_idx = ring_order[4]  # C5 is the 5th carbon in ring starting from C1
                                c5_atom = mol.GetAtomWithIdx(c5_idx)
                                # Check if C5 is chiral and get its configuration
                                if c5_atom.HasProp('_CIPCode'):
                                    c5_chirality = c5_atom.GetProp('_CIPCode')
                                    if c5_chirality == 'R':
                                        return True, "Molecule is a D-hexose; C5 has R configuration"
                                    elif c5_chirality == 'S':
                                        return False, "Molecule is an L-hexose; C5 has S configuration"
                                    else:
                                        return False, "Chirality at C5 is undefined"
                                else:
                                    return False, "C5 is not a chiral center"
        # If we've checked all potential rings without returning, unable to classify
        return False, "Could not find an appropriate cyclic structure or cannot determine chirality at C5"