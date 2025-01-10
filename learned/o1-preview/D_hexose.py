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
    This function identifies the C-5 chiral center connected to a CH2OH group and checks its configuration.

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
    
    # Find all chiral carbons
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False, useLegacyImplementation=False)
    
    if not chiral_centers:
        return False, "No chiral centers found"
    
    c5_atom_idx = None
    c5_conf = None
    # Iterate over chiral centers to find C-5
    for idx, conf in chiral_centers:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() != 6:
            continue  # Skip if not carbon
        # Check neighbors to find CH2OH group
        neighbors = atom.GetNeighbors()
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() != 6:
                continue  # Skip if not carbon
            # Check if neighbor is CH2OH group
            h_count = neighbor.GetTotalNumHs()
            o_neighbors = [nbr for nbr in neighbor.GetNeighbors() if nbr.GetAtomicNum() == 8]
            if h_count == 2 and len(o_neighbors) == 1:
                # Found CH2OH group connected to chiral carbon atom
                c5_atom_idx = idx
                c5_conf = conf
                break
        if c5_atom_idx is not None:
            break  # Stop if C-5 is found
    
    if c5_atom_idx is None:
        return False, "C-5 chiral center not found"
    
    # Check configuration at C-5
    if c5_conf == 'R':
        return True, "Molecule is a D-hexose with R configuration at C-5"
    elif c5_conf == 'S':
        return False, "Molecule is an L-hexose with S configuration at C-5"
    else:
        return False, f"Unknown configuration at C-5: {c5_conf}"