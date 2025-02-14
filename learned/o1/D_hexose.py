"""
Classifies: CHEBI:4194 D-hexose
"""
"""
Classifies: CHEBI:16551 D-hexose
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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

    # Add hydrogens to correctly identify chiral centers
    mol = Chem.AddHs(mol)

    # Check for six carbons (hexose)
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    num_carbons = len(c_atoms)
    if num_carbons != 6:
        return False, f"Number of carbons is {num_carbons}, not 6"

    # Check for monosaccharide (single ring or open chain)
    ri = mol.GetRingInfo()
    num_rings = ri.NumRings()
    if num_rings > 1:
        return False, "More than one ring detected, molecule is not a monosaccharide"

    # Check for appropriate functional groups (aldehyde or ketone group and hydroxyl groups)
    # Count hydroxyl groups
    oh_group = Chem.MolFromSmarts('[OX2H]')
    num_hydroxyls = len(mol.GetSubstructMatches(oh_group))
    if num_hydroxyls < 4:
        return False, f"Number of hydroxyl groups is {num_hydroxyls}, less than expected for a hexose"

    # Identify chiral centers
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if len(chiral_centers) < 4:
        return False, f"Found {len(chiral_centers)} chiral centers, expected at least 4 for a hexose"

    # Map atom indices to chiral centers
    chiral_dict = dict(chiral_centers)

    # Identify atom corresponding to carbon 5 (C5)
    # For this, we need to find the highest-numbered chiral carbon
    # In open-chain form, carbons are connected linearly; in cyclic forms, we need to consider ring numbering
    # Since numbering is not straightforward from SMILES, we can approximate by looking for chiral carbons connected to oxygen

    # Find all chiral carbons connected to hydroxyl groups
    c_oh_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetIdx() in chiral_dict:
            # Check if atom is connected to hydroxyl group
            neighbors = [nbr.GetAtomicNum() for nbr in atom.GetNeighbors()]
            if 8 in neighbors:
                c_oh_atoms.append(atom)

    if not c_oh_atoms:
        return False, "No chiral carbons connected to hydroxyl groups found"

    # Assume the C5 is the chiral carbon connected to oxygen (hydroxyl group) and further from aldehyde/ketone
    # Sort chiral carbons based on their distance from potential aldehyde/ketone group
    # Identify aldehyde or ketone carbon
    carbonyl = Chem.MolFromSmarts('[C]=[O]')
    carbonyl_matches = mol.GetSubstructMatches(carbonyl)
    if not carbonyl_matches:
        return False, "No aldehyde or ketone group found"

    # Get the carbon atom index of the carbonyl group
    carbonyl_carbons = [match[0] for match in carbonyl_matches]

    # Compute shortest path from carbonyl carbon to each chiral carbon
    shortest_paths = []
    for chiral_atom in c_oh_atoms:
        min_path_len = float('inf')
        for carb_coord in carbonyl_carbons:
            path = Chem.rdmolops.GetShortestPath(mol, chiral_atom.GetIdx(), carb_coord)
            if len(path) < min_path_len:
                min_path_len = len(path)
        shortest_paths.append((min_path_len, chiral_atom))

    # Sort chiral carbons by distance (we expect C5 to be the one with the longest path)
    sorted_chirals = sorted(shortest_paths, key=lambda x: -x[0])

    # Get the chiral atom assumed to be C5
    c5_atom = sorted_chirals[0][1]

    # Get the chirality tag of C5
    c5_idx = c5_atom.GetIdx()
    c5_chirality = c5_atom.GetProp('_CIPCode') if c5_atom.HasProp('_CIPCode') else None

    if c5_chirality == 'R':
        return True, "Molecule is a D-hexose; C5 has R configuration"
    elif c5_chirality == 'S':
        return False, "Molecule is an L-hexose; C5 has S configuration"
    else:
        return False, "Chirality at C5 is undefined"