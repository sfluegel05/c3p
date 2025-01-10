"""
Classifies: CHEBI:27096 trichlorobenzene
"""
from rdkit import Chem

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule is a trichlorobenzene based on its SMILES string.
    A trichlorobenzene is any chlorobenzene carrying three chloro substituents at unspecified positions.
    
    The molecule must contain at least one benzene ring with exactly three chlorine substituents
    attached directly to the ring carbons.
    
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
    
    # Initialize a flag to indicate if a suitable ring is found
    found_trichlorobenzene_ring = False
    
    # Iterate over all rings in the molecule
    for ring in atom_rings:
        # Check if the ring is aromatic and has 6 members (benzene ring)
        if len(ring) == 6 and all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            # Check if all atoms in the ring are carbons
            if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
                # Count the number of chlorine atoms attached directly to the ring carbons
                num_cl = 0
                for idx in ring:
                    atom = mol.GetAtomWithIdx(idx)
                    for neighbor in atom.GetNeighbors():
                        # Neighbor not in ring and is chlorine
                        if neighbor.GetAtomicNum() == 17 and neighbor.GetIdx() not in ring:
                            num_cl += 1
                if num_cl == 3:
                    found_trichlorobenzene_ring = True
                    break  # No need to check other rings
    
    if found_trichlorobenzene_ring:
        return True, "Contains benzene ring with exactly three chlorine substituents"
    else:
        return False, "Does not contain benzene ring with exactly three chlorine substituents"