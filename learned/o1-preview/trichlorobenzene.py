"""
Classifies: CHEBI:27096 trichlorobenzene
"""
from rdkit import Chem

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule is a trichlorobenzene based on its SMILES string.
    A trichlorobenzene is any chlorobenzene carrying three chloro substituents at unspecified positions.
    
    The molecule must contain a single benzene ring with exactly three chlorine substituents,
    and the ring should not be fused to or directly connected to any other ring systems.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a trichlorobenzene, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Iterate over all rings in the molecule
    for ring in atom_rings:
        # Check if the ring is aromatic and has 6 members
        if len(ring) == 6 and all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            # Check if all atoms in the ring are carbons
            if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
                # Check if the ring is not fused (atoms are only part of one ring)
                if all(ring_info.NumAtomRings(idx) == 1 for idx in ring):
                    # Check if the ring is not connected to other rings
                    ring_atoms = set(ring)
                    is_connected_to_other_rings = False
                    for idx in ring:
                        atom = mol.GetAtomWithIdx(idx)
                        for neighbor in atom.GetNeighbors():
                            neighbor_idx = neighbor.GetIdx()
                            # If neighbor is in a ring and not in the current ring, it's connected to another ring
                            if neighbor_idx not in ring_atoms and ring_info.NumAtomRings(neighbor_idx) > 0:
                                is_connected_to_other_rings = True
                                break
                        if is_connected_to_other_rings:
                            break
                    if is_connected_to_other_rings:
                        continue  # Skip this ring
                    # Count the number of chlorine atoms attached directly to the ring
                    num_cl = 0
                    for idx in ring:
                        atom = mol.GetAtomWithIdx(idx)
                        for neighbor in atom.GetNeighbors():
                            # Neighbor not in ring and is chlorine
                            if neighbor.GetAtomicNum() == 17 and neighbor.GetIdx() not in ring:
                                num_cl += 1
                    if num_cl == 3:
                        return True, "Contains benzene ring with exactly three chlorine substituents, not connected to other rings"
    
    # No suitable benzene ring with exactly three chlorine substituents found
    return False, "Does not contain appropriate trichlorobenzene ring"