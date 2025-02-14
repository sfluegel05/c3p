"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
"""
Classifies: tricarboxylic acid
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdMolStandardize import MetalDisconnector
from rdkit.Chem.rdMolStandardize import LargestFragmentChooser

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

    # Standardize molecule: disconnect metals, keep largest fragment
    # Disconnect metals
    md = MetalDisconnector()
    mol = md.Disconnect(mol)

    # Keep the largest organic fragment
    lfc = LargestFragmentChooser()
    mol = lfc.choose(mol)

    # Remove salts and small inorganic fragments
    if mol.GetNumAtoms() == 0:
        return False, "No atoms in molecule after removing salts and metals"

    # Check for metal ions
    metals = [3,4,11,12,13,19,20,26,29,30]
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in metals:
            return False, "Contains metal ions, likely a salt or coordination complex"

    # Identify carboxy groups (both protonated and deprotonated forms)
    carboxy_pattern = Chem.MolFromSmarts('[CX3](=O)[O,H1,-1]')
    carboxy_matches = mol.GetSubstructMatches(carboxy_pattern)
    num_carboxy_groups = len(carboxy_matches)

    if num_carboxy_groups != 3:
        return False, f"Found {num_carboxy_groups} carboxy groups, need exactly 3"

    # Check for other acidic groups (e.g., sulfonic acids)
    other_acid_pattern = Chem.MolFromSmarts('[S](=O)(=O)[O,H1,-1]')
    other_acid_matches = mol.GetSubstructMatches(other_acid_pattern)
    if other_acid_matches:
        return False, "Contains other acidic groups (e.g., sulfonic acids)"

    # Optionally, ensure molecule is organic (contains carbon)
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Molecule does not contain carbon, not an organic acid"

    return True, "Contains exactly three carboxy groups (tricarboxylic acid)"