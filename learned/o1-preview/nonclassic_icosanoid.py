"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.
    A nonclassic icosanoid is any biologically active signalling molecule made by oxygenation of C20 fatty acids
    other than the classic icosanoids (the leukotrienes and the prostanoids).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nonclassic icosanoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons != 20:
        return False, f"Molecule has {num_carbons} carbons, expected exactly 20 carbons"

    # Check for terminal carboxylic acid group (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[O;H1]')
    matches = mol.GetSubstructMatches(carboxylic_acid_pattern)

    terminal_carboxylic_acid = False
    for match in matches:
        # Check if carboxylic acid is at the end of a chain (terminal)
        idx = match[0]  # carbonyl carbon
        atom = mol.GetAtomWithIdx(idx)
        neighbors = atom.GetNeighbors()
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() == 6:  # carbon
                # Check if neighbor carbon is terminal (has only one heavy neighbor)
                heavy_neighbors = [n for n in neighbor.GetNeighbors() if n.GetAtomicNum() > 1]
                if len(heavy_neighbors) <= 2:
                    terminal_carboxylic_acid = True
                    break
        if terminal_carboxylic_acid:
            break

    if not terminal_carboxylic_acid:
        return False, "No terminal carboxylic acid group found"

    # Check for oxygenated functional groups (hydroxyl, epoxide)
    hydroxyl_pattern = Chem.MolFromSmarts('[OX2H]')
    epoxide_pattern = Chem.MolFromSmarts('[C;!R][C;R1]1[O;R1][C;R1]1[C;!R]')  # Epoxide attached to chain
    ketone_pattern = Chem.MolFromSmarts('C(=O)[C;!$(C(=O))]')  # Ketone

    has_hydroxyl = mol.HasSubstructMatch(hydroxyl_pattern)
    has_epoxide = mol.HasSubstructMatch(epoxide_pattern)
    has_ketone = mol.HasSubstructMatch(ketone_pattern)

    if not (has_hydroxyl or has_epoxide or has_ketone):
        return False, "No hydroxyl, epoxide, or ketone functional groups found"

    # Ensure there are at least two oxygen atoms in addition to carboxylic acid
    num_oxygen_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    num_oxygen_in_carboxyl = 2  # From the terminal carboxylic acid
    if num_oxygen_atoms - num_oxygen_in_carboxyl < 1:
        return False, "Not enough additional oxygen atoms in functional groups"

    # Check that molecule is mostly linear (no rings except epoxides)
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings > 0:
        # Check if all rings are epoxides (three-membered rings with one oxygen)
        for ring in ring_info.AtomRings():
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            if len(ring) != 3:
                return False, "Molecule contains non-epoxide rings"
            num_oxygen_in_ring = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8)
            if num_oxygen_in_ring != 1:
                return False, "Molecule contains non-epoxide rings"

    # Exclude prostanoids (must not contain cyclopentane ring with adjacent oxygenated groups)
    prostanoid_pattern = Chem.MolFromSmarts('C1C(CCC1)C(=O)')  # Cyclopentane ring with ketone
    if mol.HasSubstructMatch(prostanoid_pattern):
        return False, "Molecule matches prostanoid pattern (possible prostanoid)"

    # Exclude leukotrienes (conjugated triene system from C5 to C12)
    leukotriene_pattern = Chem.MolFromSmarts('C=CC=CC=CC=CC(=O)O')  # Simplified pattern
    if mol.HasSubstructMatch(leukotriene_pattern):
        return False, "Molecule matches leukotriene pattern (possible leukotriene)"

    return True, "Molecule is a nonclassic icosanoid"