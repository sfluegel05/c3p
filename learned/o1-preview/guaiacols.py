"""
Classifies: CHEBI:134251 guaiacols
"""
from rdkit import Chem

def is_guaiacols(smiles: str):
    """
    Determines if a molecule is a guaiacol based on its SMILES string.
    A guaiacol is any phenol carrying an additional methoxy substituent at the ortho-position.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a guaiacol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information
    rings = mol.GetRingInfo().AtomRings()
    for ring in rings:
        # Check for aromatic ring of size 6 (benzene ring)
        if len(ring) == 6 and all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            # Iterate over atoms in the ring
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            num_atoms = len(ring_atoms)
            for i in range(num_atoms):
                atom1 = ring_atoms[i]
                atom2 = ring_atoms[(i + 1) % num_atoms]  # Next atom in ring (adjacent)
                
                # Check if atom1 has -OH attached
                has_OH1 = False
                for nbr in atom1.GetNeighbors():
                    if nbr.GetAtomicNum() == 8 and nbr.GetDegree() == 1 and nbr.GetTotalNumHs() == 1:
                        has_OH1 = True
                        break
                # Check if atom2 has -OCH3 attached
                has_OCH3_2 = False
                for nbr in atom2.GetNeighbors():
                    if nbr.GetAtomicNum() == 8 and nbr.GetDegree() == 2:
                        # Oxygen with degree 2 may be methoxy
                        attached_methyl = False
                        oxy_neighbors = nbr.GetNeighbors()
                        if len(oxy_neighbors) == 2:
                            for nbr2 in oxy_neighbors:
                                if nbr2.GetIdx() != atom2.GetIdx():
                                    if nbr2.GetAtomicNum() == 6 and nbr2.GetDegree() == 1 and nbr2.GetTotalNumHs() == 3:
                                        # Found methyl group attached to oxygen
                                        attached_methyl = True
                                        break
                        if attached_methyl:
                            has_OCH3_2 = True
                            break
                # Check if adjacent atoms have -OH and -OCH3
                if has_OH1 and has_OCH3_2:
                    return True, "Contains phenol with methoxy group at ortho-position"
                
                # Also check the reverse (atom1 has -OCH3, atom2 has -OH)
                # Check if atom1 has -OCH3 attached
                has_OCH3_1 = False
                for nbr in atom1.GetNeighbors():
                    if nbr.GetAtomicNum() == 8 and nbr.GetDegree() == 2:
                        # Oxygen with degree 2 may be methoxy
                        attached_methyl = False
                        oxy_neighbors = nbr.GetNeighbors()
                        if len(oxy_neighbors) == 2:
                            for nbr2 in oxy_neighbors:
                                if nbr2.GetIdx() != atom1.GetIdx():
                                    if nbr2.GetAtomicNum() == 6 and nbr2.GetDegree() == 1 and nbr2.GetTotalNumHs() == 3:
                                        # Found methyl group attached to oxygen
                                        attached_methyl = True
                                        break
                        if attached_methyl:
                            has_OCH3_1 = True
                            break
                # Check if atom2 has -OH attached
                has_OH2 = False
                for nbr in atom2.GetNeighbors():
                    if nbr.GetAtomicNum() == 8 and nbr.GetDegree() == 1 and nbr.GetTotalNumHs() == 1:
                        has_OH2 = True
                        break
                if has_OCH3_1 and has_OH2:
                    return True, "Contains phenol with methoxy group at ortho-position"
    return False, "Does not contain phenol with methoxy group at ortho-position"