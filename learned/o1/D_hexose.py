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

    # Exclude molecules with carboxylic acids, esters, phosphates, sulfates, or amines
    forbidden_groups = Chem.MolFromSmarts('C(=O)[OH,OR,N,S]')
    if mol.HasSubstructMatch(forbidden_groups):
        return False, "Molecule contains forbidden functional groups (carboxylic acid, ester, amine, etc.)"

    # Identify anomeric carbon (C1): carbon connected to two oxygens, one of which is in a ring
    anomeric_carbons = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            oxygen_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
            if len(oxygen_neighbors) == 2:
                # Check if one oxygen is in a ring
                in_ring = False
                for o_atom in oxygen_neighbors:
                    if o_atom.IsInRing():
                        in_ring = True
                if in_ring and atom.IsInRing():
                    anomeric_carbons.append(atom)
    if not anomeric_carbons:
        # Try open-chain form: aldehyde group at C1
        aldehyde = Chem.MolFromSmarts('[C;H1](=O)')
        aldehyde_matches = mol.GetSubstructMatches(aldehyde)
        if not aldehyde_matches:
            return False, "No anomeric carbon or aldehyde group found"
        else:
            # Open-chain form
            # Number carbons starting from aldehyde carbon
            c1_atom = mol.GetAtomWithIdx(aldehyde_matches[0][0])
            chain = Chem.rdmolops.GetShortestPath(mol, c1_atom.GetIdx(), c_atoms[-1].GetIdx())
            if len(chain) != 6:
                return False, "Chain length is not 6 carbons"
            carbon_atoms = [mol.GetAtomWithIdx(idx) for idx in chain]
    else:
        # Cyclic form
        # Use the first anomeric carbon found
        c1_atom = anomeric_carbons[0]
        # Find the ring containing the anomeric carbon
        ring_info = mol.GetRingInfo()
        rings = ring_info.AtomRings()
        relevant_ring = None
        for ring in rings:
            if c1_atom.GetIdx() in ring and any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring):
                relevant_ring = ring
                break
        if not relevant_ring:
            return False, "No suitable ring containing the anomeric carbon found"
        # Traverse the ring to assign carbon numbers
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in relevant_ring]
        # Start numbering from the anomeric carbon
        starting_idx = relevant_ring.index(c1_atom.GetIdx())
        reordered_ring = ring_atoms[starting_idx:] + ring_atoms[:starting_idx]
        # Number of carbons in the ring (excluding oxygen)
        ring_carbon_atoms = [atom for atom in reordered_ring if atom.GetAtomicNum() == 6]
        if len(ring_carbon_atoms) not in [4, 5]:  # Furanose or Pyranose ring
            return False, f"Ring size {len(ring_carbon_atoms)} not typical for hexose"
        # Number the carbons
        carbon_atoms = []
        for atom in ring_carbon_atoms:
            carbon_atoms.append(atom)
        # For pyranose, C6 is connected to C5 outside the ring
        if len(ring_carbon_atoms) == 5:
            # Find C6 connected to C5
            c5_atom = carbon_atoms[4]
            for neighbor in c5_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor not in carbon_atoms:
                    carbon_atoms.append(neighbor)  # C6
        elif len(ring_carbon_atoms) == 4:
            # Furanose: C5 is outside the ring, connected to C4
            c4_atom = carbon_atoms[3]
            for neighbor in c4_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor not in carbon_atoms:
                    carbon_atoms.append(neighbor)  # C5
                    # Find C6 connected to C5
                    c5_atom = carbon_atoms[4]
                    for nbr in c5_atom.GetNeighbors():
                        if nbr.GetAtomicNum() == 6 and nbr not in carbon_atoms:
                            carbon_atoms.append(nbr)  # C6
                    break
        if len(carbon_atoms) != 6:
            return False, "Could not identify all six carbon atoms"
    # Identify chiral centers
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    chiral_dict = dict(chiral_centers)
    # Get the C5 atom
    c5_atom = carbon_atoms[4]
    c5_idx = c5_atom.GetIdx()
    # Check if C5 is chiral
    if c5_idx not in chiral_dict:
        return False, "C5 is not a chiral center"
    c5_chirality = c5_atom.GetProp('_CIPCode') if c5_atom.HasProp('_CIPCode') else None
    if c5_chirality == 'R':
        return True, "Molecule is a D-hexose; C5 has R configuration"
    elif c5_chirality == 'S':
        return False, "Molecule is an L-hexose; C5 has S configuration"
    else:
        return False, "Chirality at C5 is undefined"