"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
"""
Classifies: tricarboxylic acid
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid based on its SMILES string.
    A tricarboxylic acid is an oxoacid containing exactly three carboxy groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tricarboxylic acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the carboxy group pattern (both protonated and deprotonated forms)
    carboxy_pattern = Chem.MolFromSmarts('[CX3](=O)[O;H1,-]')
    if carboxy_pattern is None:
        return False, "Failed to create carboxy group pattern"

    # Find all carboxy groups in the molecule
    carboxy_matches = mol.GetSubstructMatches(carboxy_pattern)
    num_carboxy_groups = len(carboxy_matches)

    if num_carboxy_groups != 3:
        return False, f"Found {num_carboxy_groups} carboxy groups, need exactly 3"

    # Check if the molecule is an oxoacid (contains oxygen and hydrogen in acidic groups)
    # Presence of carboxy groups ensures it's an oxoacid
    is_oxoacid = True  # Since carboxy groups are oxoacid functional groups

    if not is_oxoacid:
        return False, "Molecule is not an oxoacid"

    return True, "Contains exactly three carboxy groups (tricarboxylic acid)"