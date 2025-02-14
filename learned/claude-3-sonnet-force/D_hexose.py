"""
Classifies: CHEBI:4194 D-hexose
"""
"""
Classifies: CHEBI:16646 D-hexose
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    A D-hexose is a hexose that has D-configuration at position 5.

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

    # Check if molecule contains exactly 6 carbon atoms and 6 oxygen atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count != 6 or o_count != 6:
        return False, "Not a hexose (does not contain 6 carbons and 6 oxygens)"

    # Check for the presence of a ring structure
    ring_info = mol.GetRingInfo()
    if not ring_info.AtomRings():
        return False, "Not a hexose (no ring structure found)"

    # Get the chiral centers of the molecule
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnspecified=True)

    # Check if there are at least 4 chiral centers (minimum for a hexose)
    if len(chiral_centers) < 4:
        return False, "Not a hexose (less than 4 chiral centers)"

    # Find the second-to-last chiral center (position 5 in a hexose)
    second_to_last_chiral = chiral_centers[-2]

    # Check if the second-to-last chiral center has D-configuration
    if second_to_last_chiral == "R":
        return True, "Molecule has D-configuration at position 5 (second-to-last chiral center)"
    else:
        return False, "Molecule does not have D-configuration at position 5 (second-to-last chiral center)"