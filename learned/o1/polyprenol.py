"""
Classifies: CHEBI:26199 polyprenol
"""
"""
Classifies: polyprenol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    A polyprenol is a long-chain isoprenoid alcohol composed of more than one isoprene unit.
    They have the general formula H-[CH2C(Me)=CHCH2]nOH, with n > 1.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for terminal primary alcohol group (O-H connected to a primary carbon)
    alcohol_pattern = Chem.MolFromSmarts('[C;H2][O;H1]')  # Primary alcohol
    if not mol.HasSubstructMatch(alcohol_pattern):
        return False, "No terminal primary alcohol group found"

    # Define isoprene unit pattern (captures the repeating isoprene unit in polyprenols)
    isoprene_pattern = Chem.MolFromSmarts('[CH2]-[C](=C([CH3]))-[CH2]')  # Isoprene unit
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    num_isoprene_units = len(isoprene_matches)

    if num_isoprene_units < 2:
        return False, f"Found {num_isoprene_units} isoprene units, need more than one"

    return True, f"Contains terminal alcohol and {num_isoprene_units} isoprene units"