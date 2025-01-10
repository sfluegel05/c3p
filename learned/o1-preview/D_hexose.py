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
    
    # Identify chiral centers
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False, useLegacyImplementation=False)
    if not chiral_centers:
        return False, "No chiral centers found"
    
    # Initialize variables
    found_c5 = False
    c5_conf = None
    
    # Iterate over chiral centers to find C-5
    for idx, conf in chiral_centers:
        atom = mol.GetAtomWithIdx(idx)
        
        # Check if atom is a carbon
        if atom.GetAtomicNum() != 6:
            continue
        
        # Check for connections to oxygen (part of ring) and a carbon (possible CH2OH group)
        neighbor_atomic_nums = [nbr.GetAtomicNum() for nbr in atom.GetNeighbors()]
        if neighbor_atomic_nums.count(8) < 1 or neighbor_atomic_nums.count(6) < 1:
            continue
        
        # Look for connected CH2OH group
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                # Count hydrogens and oxygens attached to this neighbor
                h_count = sum(1 for nbr in neighbor.GetNeighbors() if nbr.GetAtomicNum() == 1)
                o_count = sum(1 for nbr in neighbor.GetNeighbors() if nbr.GetAtomicNum() == 8)
                if h_count == 2 and o_count == 1:
                    # Found C-5
                    c5_conf = conf
                    found_c5 = True
                    break
        if found_c5:
            break
    
    if not found_c5:
        return False, "C-5 chiral center not found"
    
    # Check configuration at C-5
    if c5_conf == 'R':
        return True, "Molecule is a D-hexose with R configuration at C-5"
    elif c5_conf == 'S':
        return False, "Molecule is an L-hexose with S configuration at C-5"
    else:
        return False, f"Cannot determine configuration at C-5: {c5_conf}"