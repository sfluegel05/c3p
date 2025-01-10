"""
Classifies: CHEBI:4194 D-hexose
"""
"""
Classifies: D-hexose
"""

from rdkit import Chem

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    A D-hexose is a hexose (six-carbon monosaccharide) that has D-configuration at position 5.

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
    
    # Assign stereochemistry
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Check for 6 carbons
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons != 6:
        return False, "Molecule does not have 6 carbons"
    
    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Initialize variables
    c5_atom_idx = None
    c5_conf = None
    
    # Identify rings containing oxygen atoms
    for ring in atom_rings:
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        ring_oxygen_atoms = [atom for atom in ring_atoms if atom.GetAtomicNum() == 8]
        
        # Process rings with exactly one ring oxygen (sugar rings)
        if len(ring_oxygen_atoms) == 1:
            ring_oxygen = ring_oxygen_atoms[0]
            # Find carbons connected to the ring oxygen
            oxygen_neighbors = ring_oxygen.GetNeighbors()
            for carbon in oxygen_neighbors:
                if carbon.GetAtomicNum() != 6:
                    continue
                # Check if this carbon (candidate C-5) is connected to:
                # - the ring oxygen
                # - another ring carbon
                # - a CH2OH group (carbon with two hydrogens and one oxygen)
                neighbors = carbon.GetNeighbors()
                ring_connections = 0
                has_ch2oh = False
                for neighbor in neighbors:
                    if neighbor.GetIdx() == ring_oxygen.GetIdx():
                        continue  # Skip the ring oxygen
                    elif neighbor.GetIdx() in ring:
                        ring_connections += 1
                    elif neighbor.GetAtomicNum() == 6:
                        # Potential CH2OH group
                        o_count = sum(1 for nbr in neighbor.GetNeighbors() if nbr.GetAtomicNum() == 8)
                        h_count = neighbor.GetTotalNumHs()
                        if o_count == 1 and h_count == 2:
                            has_ch2oh = True
                # Candidate C-5 must have:
                # - connected to one other ring carbon (ring_connections == 1)
                # - connected to CH2OH group (has_ch2oh == True)
                if ring_connections == 1 and has_ch2oh:
                    c5_atom_idx = carbon.GetIdx()
                    break
            if c5_atom_idx is not None:
                break
    
    if c5_atom_idx is None:
        return False, "C-5 chiral center not found"
    
    # Get chiral tags and configuration
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False, useLegacyImplementation=False)
    chiral_dict = {idx: conf for idx, conf in chiral_centers}
    
    if c5_atom_idx not in chiral_dict:
        return False, "C-5 is not a chiral center"
    
    c5_conf = chiral_dict[c5_atom_idx]
    
    # Check configuration at C-5
    if c5_conf == 'R':
        return True, "Molecule is a D-hexose with R configuration at C-5"
    elif c5_conf == 'S':
        return False, "Molecule is an L-hexose with S configuration at C-5"
    else:
        return False, f"Cannot determine configuration at C-5: {c5_conf}"