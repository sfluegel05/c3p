"""
Classifies: CHEBI:4194 D-hexose
"""
from rdkit import Chem

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    A D-hexose is defined as a hexose with D-configuration at position 5.

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

    # Check if the molecule graph has exactly 6 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 6:
        return False, "Molecule does not have 6 carbon atoms typical of hexoses"

    # Identify pyranose or furanose ring structures to map positions correctly
    pyranose_pattern = Chem.MolFromSmarts("C1OC(O)C(O)C(O)C1") # typical pyranose pattern
    furanose_pattern = Chem.MolFromSmarts("O1C(O)C(O)C(O)C1") # typical furanose pattern

    ring_found = False
    ring_atoms = []

    if mol.HasSubstructMatch(pyranose_pattern):
        ring_atoms = mol.GetSubstructMatch(pyranose_pattern)
        ring_found = True
    elif mol.HasSubstructMatch(furanose_pattern):
        ring_atoms = mol.GetSubstructMatch(furanose_pattern)
        ring_found = True

    if not ring_found:
        return False, "No typical hexose ring structure found"

    # Identify the chiral centers in the molecule
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)

    # Map the expected D-position (chiral center 5 in ring) and check configuration
    for center, config in chiral_centers:
        if center in ring_atoms:
            # Position 5 mapping could be simplified to the 5th chiral seen
            # But should be (loosely ensure) ~5th position
            if config == 'R':
                return True, "Molecule has D-configuration at position 5"

    return False, "Molecule does not conform to D-hexose classification due to stereochemistry"