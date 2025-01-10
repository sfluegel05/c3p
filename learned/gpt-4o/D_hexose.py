"""
Classifies: CHEBI:4194 D-hexose
"""
from rdkit import Chem

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    A D-hexose is defined as a hexose that has D-configuration at position 5.

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

    # Ensure the molecule contains 6 carbons typical of hexoses
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 6:
        return False, "Molecule does not have 6 carbon atoms typical of hexoses"

    # Attempt to identify the ring or linear structure
    pyranose_pattern = Chem.MolFromSmarts("C1OC(O)C(O)C(O)C1") # typical pyranose pattern
    furanose_pattern = Chem.MolFromSmarts("O1C(O)C(O)C(O)C1") # typical furanose pattern

    found_pyranose = mol.HasSubstructMatch(pyranose_pattern)
    found_furanose = mol.HasSubstructMatch(furanose_pattern)
    ring_found = found_pyranose or found_furanose

    # Identify the chiral centers in the molecule
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)

    # Tracking whether we find a valid D-configuration
    valid_d_configuration = False

    # Look specifically for D-configuration at position 5
    for center, config in chiral_centers:
        # Checking for manual or structured D-configuration recognition
        atom = mol.GetAtomWithIdx(center)
        if atom.GetIdx() == 4:  # Check the 5th carbon position, accounting for 0-based index
            if config == 'R':
                # Assuming 'R' as indicative for D-configuration requires domain knowledge
                valid_d_configuration = True

    # Check the stereochemistry of the positions outside typical ring paradigms if necessary
    if not valid_d_configuration and not ring_found:
        # For non-ring structures (aldehyde forms), further stereochemical insight required
        for center, config in chiral_centers:
            if config == 'R':
                valid_d_configuration = True
                break

    if valid_d_configuration:
        return True, "Molecule has D-configuration at position 5"

    if ring_found:
        return False, "Molecule does not have D-configuration at the expected position"

    return False, "Molecule does not conform to D-hexose classification due to stereochemistry"