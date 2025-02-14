"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_anthocyanidin_cation(smiles: str):
    """
    Determines if a molecule is an anthocyanidin cation based on its SMILES string.
    An anthocyanidin cation is an aglycon of anthocyanin cation, characterized by a flavylium core.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthocyanidin cation, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the flavylium core SMARTS pattern
    # A more general pattern that matches the core structure
    # Specifically: an O+ connected to an aromatic carbon, part of a 6-member ring,
    # and attached to another aromatic ring.
    flavylium_core_pattern = Chem.MolFromSmarts("[O+1]1[c]~[c]~[c]~[c]~[c]1-[c]~[c]:[c]:[c]:[c]:[c]")


    # Check for the flavylium core
    if not mol.HasSubstructMatch(flavylium_core_pattern):
        return False, "Flavylium core not found"

    return True, "Contains the flavylium cation core (aglycone)"